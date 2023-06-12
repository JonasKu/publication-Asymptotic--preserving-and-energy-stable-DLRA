__precompile__

using ProgressMeter
using LinearAlgebra, FastGaussQuadrature
using LegendrePolynomials
using QuadGK
using SparseArrays

struct Solver
    # spatial grid of cell interfaces
    x::Array{Float64};

    # Solver settings
    settings::Settings;
    
    # squared L2 norms of Legendre coeffs
    gamma::Array{Float64,1};
    # flux matrix PN system
    A::Array{Float64,2};
    ABar::Array{Float64,1};
    # Roe matrix
    AbsA::Array{Float64,2};
    M::Array{Float64,2};

    # stencil matrices
    Dx::Tridiagonal{Float64, Vector{Float64}};
    Dxx::Tridiagonal{Float64, Vector{Float64}};
    DxS::SparseMatrixCSC{Float64, Int64};
    DxST::SparseMatrixCSC{Float64, Int64};

    # physical parameters
    sigmaA::Float64;
    sigmaS::Float64;

    # constructor
    function Solver(settings)
        x = settings.x;
        nx = settings.NCells;
        dx = settings.dx;

        # setup flux matrix
        gamma = ones(settings.nPN+1);

        # setup gamma vector
        gamma = zeros(settings.nPN+1);
        for i = 1:settings.nPN+1
            n = i-1;
            gamma[i] = 2/(2*n+1);
        end
        
        # setup flux matrix
        AFull = zeros(settings.nPN+1,settings.nPN+1)

        for i = 1:settings.nPN
            n = i-1;
            AFull[i,i+1] = (n+1)/(2*n+1)*sqrt(gamma[i+1])/sqrt(gamma[i]);
        end

        for i = 2:settings.nPN+1
            n = i-1;
            AFull[i,i-1] = n/(2*n+1)*sqrt(gamma[i-1])/sqrt(gamma[i]);
        end
        ABar = AFull[1,2:end];
        A = AFull[2:end,2:end];

        # setup Roe matrix according to quadrature
        nq = 2*settings.nPN+1
        TFull = zeros(settings.nPN+1,settings.nPN+1) # allocate transformation matrix
        M = zeros(settings.nPN+1,nq) # allocate transformation matrix
        mu, w = gausslegendre(settings.nPN+1)
        for k = 1:settings.nPN+1
            P = collectPl(mu[k],lmax = settings.nPN);
            for i = 1:settings.nPN+1
                TFull[i,k] = P[i-1]*sqrt(w[k])/sqrt(gamma[i]);
                
            end
        end

        muF, wF = gausslegendre(nq)
        for k = 1:nq
            P = collectPl(muF[k],lmax = settings.nPN);
            for i = 1:settings.nPN+1
                M[i,k] = P[i-1]*wF[k]/sqrt(gamma[i]);
            end
        end

        
        T = TFull[2:end,:];
        #AbsA = T*abs.(diagm([0.0; mu[2:end]]))*T';
        AbsA = T*abs.(diagm(mu))*T';
        A1 = T*diagm([0.0; mu[2:end]])*T';

        # setup Roe matrix
        #S = eigvals(A)
        #V = eigvecs(A)
        #AbsA = V*abs.(diagm(S))*inv(V)

        # set up spatial stencil matrices
        Dx = Tridiagonal(-ones(nx-1)./dx/2.0,zeros(nx),ones(nx-1)./dx/2.0) # central difference matrix
        Dxx = Tridiagonal(ones(nx-1)./dx/2.0,-ones(nx)./dx,ones(nx-1)./dx/2.0) # stabilization matrix
        DxS = Tridiagonal(zeros(nx),-ones(nx+1)./dx,ones(nx)./dx) # difference matrix for staggered grid
        DxS = DxS[1:nx,:]
        DxST = Tridiagonal([-ones(nx-1)./dx;0.0],[0.0;ones(nx-1)./dx; 0.0],zeros(nx)) # difference matrix for staggered grid
        DxST = DxST[:,1:(end-1)];

        new(x,settings,gamma,A,ABar,AbsA,M,Dx,Dxx,DxS,DxST,settings.sigmaA,settings.sigmaS);
    end
end

function SetupIC(obj::Solver)
    u = zeros(obj.settings.NCells,obj.settings.nPN); # u has NCells cells, rho has NCells+1 cells
    rho = 2.0/sqrt(obj.gamma[1])*IC(obj.settings,obj.settings.x);

    nq = 2*obj.settings.nPN+1
    mu, w = gausslegendre(nq)
    v = zeros(nq)

    s1 = 0.03
    s2 = s1^2
    x0 = 0;
    x = obj.settings.xMid
    mu0 = 1.0;
    for k = 1:nq
        v[k] = 1.0/(sqrt(2*pi)*s1) *exp(-((mu[k]-mu0)*(mu[k]-mu0))/2.0/s2);
    end
    moments = obj.M*v;
    for j = 1:obj.settings.NCells
        u[j,:] = 1.0/(sqrt(2*pi)*s1) *exp(-((x[j]-x0)*(x[j]-x0))/2.0/s2) * moments[2:end] / obj.settings.epsilon;
    end
    x = obj.settings.x
    for j = 1:length(x)
        rho[j] = 1.0/(sqrt(2*pi)*s1) *exp(-((x[j]-x0)*(x[j]-x0))/2.0/s2) * moments[1];
    end
    return rho,u;
end

function SolveDiffusion(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    ABar = 0.5773502691896257;
    sigmaS = obj.settings.sigmaS;
    dt = (dx^2 * sigmaS)/(2+0.5*ABar^2)
    tEnd = obj.settings.tEnd;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);

    #Compute diagonal of scattering matrix G
    sigS = ones(nx).*obj.settings.sigmaS;
    sigmaS=Diagonal(sigS);
    sigmaA=Diagonal(ones(nx+1)).*obj.settings.sigmaA;
    ABar = obj.ABar;

    energy = zeros(2, nt)

    prog = Progress(nt,1)
    #loop over time
    for n=1:nt
        energy[1,n] = (n-1) * dt
        energy[2,n] = dx * norm(rho)^2
        u = -Diagonal(1.0./sigS)*obj.DxS*rho*ABar';
        rho = rho .- dt * obj.DxST*u*ABar .- dt * sigmaA*rho;

        rho[1] = rho[end-1];
        rho[end] = rho[2];
        
        

        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho, energy;

end

function Solve(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt
    dx = obj.settings.dx;

    N = obj.settings.nPN;
    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);

    #Compute diagonal of scattering matrix G
    G = Diagonal(ones(N));
    sigmaS=Diagonal(ones(nx)).*obj.settings.sigmaS;
    sigmaA=Diagonal(ones(nx)).*obj.settings.sigmaA;
    sigmaAR=Diagonal(ones(nx+1)).*obj.settings.sigmaA;
    A = obj.A;
    ABar = obj.ABar;
    AbsA = obj.AbsA;
    epsilon = obj.settings.epsilon;

    energy = zeros(2, nt)

    prog = Progress(nt,1)
    #loop over time
    for n=1:nt
        energy[1,n] = (n-1) * dt
        energy[2,n] = dx * norm(rho)^2 + dx * epsilon^2 * norm(u)^2
        u = u .- dt * obj.Dx*u*A'./epsilon .+ dt * obj.Dxx*u*AbsA'./epsilon .- dt*sigmaA*u .- dt*sigmaS*u*G./epsilon^2 .- dt*obj.DxS*rho*ABar'/epsilon^2;
        rho = rho .- dt * obj.DxST*u*ABar .- dt * sigmaAR*rho;
        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho,0.5*sqrt(obj.gamma[1])*u, energy;

end

function SolveIMEX(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt
    dx = obj.settings.dx;

    N = obj.settings.nPN;
    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);

    #Compute diagonal of scattering matrix G
    sigmaSVec = ones(nx).*obj.settings.sigmaS;
    sigmaS=Diagonal(sigmaSVec);
    sigmaA=Diagonal(ones(nx)).*obj.settings.sigmaA;
    sigmaAR=Diagonal(ones(nx+1)).*obj.settings.sigmaA;
    A = obj.A;
    ABar = obj.ABar;
    AbsA = obj.AbsA;
    epsilon = obj.settings.epsilon;
    Scat = Diagonal((1 .+ dt./epsilon^2 .* sigmaSVec).^(-1));

    energy = zeros(2, nt)
    prog = Progress(nt,1);
    #loop over time
    for n=1:nt
        energy[1,n] = (n-1) * dt
        energy[2,n] = dx * norm(rho)^2 + dx * epsilon^2 * norm(u)^2
        u .= u .- dt * obj.Dx*u*A'./epsilon .+ dt * obj.Dxx*u*AbsA'./epsilon .- dt*sigmaA*u .- dt*obj.DxS*rho*ABar'/epsilon^2;

        u = Scat*u;

        rho = rho .- dt * obj.DxST*u*ABar .- dt * sigmaAR*rho;

        rho[1] = rho[end-1]; rho[end] = rho[2];
        u[1,:] = u[end-1,:]; u[end,:] = u[2,:];

        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho,0.5*sqrt(obj.gamma[1])*u, energy;

end

function SolveDLRA(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;
    r = obj.settings.r;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt
    dx = obj.settings.dx;

    N = obj.settings.nPN;
    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);
    sigmaSVec = ones(nx).*obj.settings.sigmaS;
    sigmaS=Diagonal(sigmaSVec);
    epsilon = obj.settings.epsilon;
    sigmaA=Diagonal(ones(nx)).*obj.settings.sigmaA;
    sigmaAR=Diagonal(ones(nx+1)).*obj.settings.sigmaA;
    Scat = Diagonal((1 .+ dt./epsilon^2 .* sigmaSVec).^(-1));
    A = obj.A;
    ABar = obj.ABar;
    AbsA = obj.AbsA;
    

    #if norm(u) < 1e-7
    #    X,S,V = svd(rho*[1.0;zeros(N-1)]'); # make sure initial information in basis has some structure
    #    S = zeros(size(S));
    #else
    X,S,V = svd(u);
    #end
    
    # rank-r truncation:
    X = X[:,1:r];
    V = V[:,1:r];
    S = diagm(S[1:r]);

    # allocate memory for efficiency
    K = zeros(size(X));
    L = zeros(size(V));
    MUp = zeros(r,r)
    NUp = zeros(r,r)

    energy = zeros(2, nt)
    prog = Progress(nt,1);
    #loop over time
    for n=1:nt

        energy[1,n] = (n-1) * dt
        energy[2,n] = dx * norm(rho)^2 + dx * epsilon^2 * norm(S)^2

        ################## K-step ##################

        K .= X*S;

        VAV = V'*A'*V;
        VAbsAV = V'*AbsA'*V;
        ABarV = ABar'*V;

        K .= K .- dt * obj.Dx*K*VAV./epsilon .+ dt * obj.Dxx*K*VAbsAV./epsilon .- dt*sigmaA*K .- dt*obj.DxS*rho*ABarV/epsilon^2;

        # apply implicit scattering
        K = Scat*K;

        XNew,STmp = qr(K);
        XNew = Matrix(XNew)
        XNew = XNew[:,1:r];
        XNew[1,:] .= zeros(r); XNew[end,:] .= zeros(r);

        MUp .= XNew' * X;

        ################## L-step ##################
        L = V*S';
        XDxX = (X'*obj.Dx*X)';
        XDxxX = X'*obj.Dxx*X;
        XsigAX = X'*sigmaA*X;
        XsigSX = X'*sigmaS*X;
        DxX = obj.DxS'*X;

        L .= L .- dt * A*L*XDxX./epsilon .+ dt * AbsA*L*XDxxX./epsilon .- dt*L*XsigAX .- dt*ABar*rho'*DxX/epsilon^2;

        for k = 1:N
            L[k,:] =  (I+dt/epsilon^2 * XsigSX) \ L[k,:];
        end

        VNew,STmp = qr(L);
        VNew = Matrix(VNew)
        VNew = VNew[:,1:r];

        NUp .= VNew' * V;

        ################## S-step ##################
        S .= MUp*S*(NUp')
        V .= VNew;
        X .= XNew;

        VAV = V'*A'*V;
        VAbsAV = V'*AbsA'*V;
        ABarV = ABar'*V;
        XDxX = X'*obj.Dx*X;
        XDxxX = X'*obj.Dxx*X;
        XsigAX = X'*sigmaA*X;
        XsigSX = X'*sigmaS*X;
        XDx = X'*obj.DxS;

        S .= S .- dt * XDxX*S*VAV./epsilon .+ dt * XDxxX*S*VAbsAV./epsilon .- dt*XsigAX*S .- dt*XDx*rho*ABarV/epsilon^2;

        for k = 1:r
            S[:,k] =  (I+dt/epsilon^2 * XsigSX) \ S[:,k];
        end

        ################## rho-step ##################
        DxX = obj.DxST*X;
        rho = rho .- dt * DxX*S*ABarV' .- dt * sigmaAR*rho;

        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho,0.5*sqrt(obj.gamma[1])*X*S*V', energy;

end