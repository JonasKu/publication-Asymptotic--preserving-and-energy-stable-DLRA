__precompile__

using ProgressMeter
using LinearAlgebra
using LegendrePolynomials
using QuadGK

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

    # stencil matrices
    Dx::Tridiagonal{Float64, Vector{Float64}};
    Dxx::Tridiagonal{Float64, Vector{Float64}};

    # physical parameters
    sigmaA::Float64;
    sigmaS::Float64;

    # constructor
    function Solver(settings)
        x = settings.x;
        nx = settings.NCells;
        dx = settings.dx;

        # setup flux matrix
        gamma = ones(settings.nPN);

        # setup gamma vector
        gamma = zeros(settings.nPN);
        for i = 1:settings.nPN
            n = i-1;
            gamma[i] = 2/(2*n+1);
        end
        
        # setup flux matrix
        AFull = zeros(settings.nPN+1,settings.nPN+1)

        for i = 1:(settings.nPN-1)
            n = i-1;
            AFull[i,i+1] = (n+1)/(2*n+1)*sqrt(gamma[i+1])/sqrt(gamma[i]);
        end

        for i = 2:settings.nPN
            n = i-1;
            AFull[i,i-1] = n/(2*n+1)*sqrt(gamma[i-1])/sqrt(gamma[i]);
        end
        ABar = AFull[1,2:end];
        A = AFull[2:end,2:end];

        # setup Roe matrix
        S = eigvals(A)
        V = eigvecs(A)
        AbsA = V*abs.(diagm(S))*inv(V)

        # set up spatial stencil matrices
        Dx = Tridiagonal(-ones(nx-1)./dx/2.0,zeros(nx),ones(nx-1)./dx/2.0) # central difference matrix
        Dxx = Tridiagonal(ones(nx-1)./dx/2.0,-ones(nx)./dx,ones(nx-1)./dx/2.0) # stabilization matrix

        new(x,settings,gamma,A,ABar,AbsA,Dx,Dxx,settings.sigmaA,settings.sigmaS);
    end
end

function SetupIC(obj::Solver)
    u = zeros(obj.settings.NCells,obj.settings.nPN); 
    rho = 2.0/sqrt(obj.gamma[1])*IC(obj.settings,obj.settings.xMid);
    return rho,u;
end

function SolveDiffusion(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    ABar = 0.5773502691896257;
    sigmaS = obj.settings.sigmaS;
    dt = 5*(dx^2 * sigmaS)/(2+0.5*ABar^2)
    tEnd = obj.settings.tEnd;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);

    #Compute diagonal of scattering matrix G
    sigS = ones(nx).*obj.settings.sigmaS;
    sigmaS=Diagonal(sigS);
    sigmaA=Diagonal(ones(nx)).*obj.settings.sigmaA;
    ABar = obj.ABar;

    # create correct diffusion matrix
    Dxx = 2.0*obj.Dxx/dx;

    prog = Progress(nt,1)
    A = ABar'*ABar;
    #loop over time
    for n=1:nt
        rho = rho .+ dt * Dxx*rho*A .- 0.0*dt*obj.Dx*Diagonal(rho)*obj.Dx*(1.0./sigS)*A .- 0.0*dt*Diagonal(rho)*A*Dxx*(1.0./sigS) .- dt * sigmaA*rho;
        rho[1] = rho[end-1];
        rho[end] = rho[2];
        
        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho;

end

function Solve(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    N = obj.settings.nPN;
    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);

    #Compute diagonal of scattering matrix G
    G = Diagonal(ones(N));
    sigmaS=Diagonal(ones(nx)).*obj.settings.sigmaS;
    sigmaA=Diagonal(ones(nx)).*obj.settings.sigmaA;
    A = obj.A;
    ABar = obj.ABar;
    AbsA = obj.AbsA;
    epsilon = obj.settings.epsilon;

    prog = Progress(nt,1)
    #loop over time
    for n=1:nt
        u = u .- dt * obj.Dx*u*A'./epsilon .+ dt * obj.Dxx*u*AbsA'./epsilon .- dt*sigmaA*u .- dt*sigmaS*u*G./epsilon^2 .- dt*obj.Dx*rho*ABar'/epsilon^2;
        rho = rho .- dt * obj.Dx*u*ABar .- dt * sigmaA*rho;
        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho,0.5*sqrt(obj.gamma[1])*u;

end

function SolveIMEX(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    N = obj.settings.nPN;
    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);

    #Compute diagonal of scattering matrix G
    sigmaSVec = ones(nx).*obj.settings.sigmaS;
    sigmaS=Diagonal(sigmaSVec);
    sigmaA=Diagonal(ones(nx)).*obj.settings.sigmaA;
    A = obj.A;
    ABar = obj.ABar;
    AbsA = obj.AbsA;
    epsilon = obj.settings.epsilon;
    Scat = Diagonal((1 .+ dt./epsilon^2 .* sigmaSVec).^(-1));

    prog = Progress(nt,1);
    #loop over time
    for n=1:nt
        u .= u .- dt * obj.Dx*u*A'./epsilon .+ dt * obj.Dxx*u*AbsA'./epsilon .- dt*sigmaA*u .- dt*obj.Dx*rho*ABar'/epsilon^2;

        u = Scat*u;

        rho = rho .- dt * obj.Dx*u*ABar .- dt * sigmaA*rho;

        rho[1] = rho[end-1]; rho[end] = rho[2];
        u[1,:] = u[end-1,:]; u[end,:] = u[2,:];

        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho,0.5*sqrt(obj.gamma[1])*u;

end

function SolveDLRA(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;
    r = obj.settings.r;

    nt = Int(ceil(tEnd/dt));     # number of time steps
    dt = obj.settings.tEnd/nt;           # adjust dt

    N = obj.settings.nPN;
    nx = obj.settings.NCells;

    # Set up initial condition
    rho,u = SetupIC(obj);
    sigmaSVec = ones(nx).*obj.settings.sigmaS;
    sigmaS=Diagonal(sigmaSVec);
    epsilon = obj.settings.epsilon;
    sigmaA=Diagonal(ones(nx)).*obj.settings.sigmaA;
    Scat = Diagonal((1 .+ dt./epsilon^2 .* sigmaSVec).^(-1));
    A = obj.A;
    ABar = obj.ABar;
    AbsA = obj.AbsA;
    

    if norm(u) < 1e-7
        X,S,V = svd(rho*[1.0;zeros(N-1)]'); # make sure initial information in basis has some structure
        S = zeros(size(S));
    else
        X,S,V = svd(u);
    end
    
    # rank-r truncation:
    X = X[:,1:r];
    V = V[:,1:r];
    S = diagm(S[1:r]);

    # allocate memory for efficiency
    K = zeros(size(X));
    L = zeros(size(V));
    MUp = zeros(r,r)
    NUp = zeros(r,r)

    prog = Progress(nt,1);
    #loop over time
    for n=1:nt

        ################## K-step ##################

        K .= X*S;

        VAV = V'*A'*V;
        VAbsAV = V'*AbsA'*V;
        ABarV = ABar'*V;

        K .= K .- dt * obj.Dx*K*VAV./epsilon .+ dt * obj.Dxx*K*VAbsAV./epsilon .- dt*sigmaA*K .- dt*obj.Dx*rho*ABarV/epsilon^2;

        # apply implicit scattering
        K = Scat*K;

        XNew,STmp = qr(K);
        XNew = Matrix(XNew)
        XNew = XNew[:,1:r];
        XNew[1,:] .= zeros(r); XNew[end,:] .= zeros(r);

        MUp .= XNew' * X;

        ################## L-step ##################
        L = V*S';
        XDxX = X'*obj.Dx'*X;
        XDxxX = X'*obj.Dxx'*X;
        XsigAX = X'*sigmaA*X;
        XsigSX = X'*sigmaS*X;
        DxX = obj.Dx'*X;

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
        XDx = X'*obj.Dx;

        S .= S .- dt * XDxX*S*VAV./epsilon .+ dt * XDxxX*S*VAbsAV./epsilon .- dt*XsigAX*S .- dt*XDx*rho*ABarV/epsilon^2;

        for k = 1:r
            S[:,k] =  (I+dt/epsilon^2 * XsigSX) \ S[:,k];
        end

        ################## rho-step ##################
        DxX = obj.Dx*X;
        rho = rho .- dt * DxX*S*ABarV' .- dt * sigmaA*rho;

        next!(prog) # update progress bar
    end
    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*rho,0.5*sqrt(obj.gamma[1])*X*S*V';

end
