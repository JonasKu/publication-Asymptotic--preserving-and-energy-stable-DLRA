__precompile__

using FastGaussQuadrature, LinearAlgebra

mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    # number spatial cells
    NCells::Int64;
    # start and end point
    a::Float64;
    b::Float64;
    # grid cell width
    dx::Float64

    # time settings
    # end time
    tEnd::Float64;
    # time increment
    dt::Float64;
    # CFL number 
    cfl::Float64;
    
    # degree PN
    nPN::Int64;

    # spatial grid
    x
    xMid

    # physical parameters
    sigmaA::Float64;
    sigmaS::Float64;
    epsilon::Float64;

    # low rank parameters
    r::Int;

    function Settings(Nx::Int=502,nPN=100,r=10,epsilon=1.0,problem::String="LineSource")
        # spatial grid setting
        NCells = Nx - 1;
        a = -1.5; # left boundary
        b = 1.5; # right boundary
        
        # time settings
        tEnd = 0.2;
        cfl = 1.0; # CFL condition
        
        x = collect(range(a,stop = b,length = NCells));
        dx = x[2]-x[1];
        x = [x[1]-dx;x]; # add ghost cells so that boundary cell centers lie on a and b
        x = x.+dx/2;
        xMid = x[1:(end-1)].+0.5*dx

        # physical parameters
        sigmaS = 1.0;
        sigmaA = 0.0; 

        # compute CFL condition
        mu, w = gausslegendre(nPN+1)
        dt = 1000;
        for k = 1:(nPN+1)
            dtTmp = 1/(4*(nPN+1)*w[k]) * (epsilon * dx / abs(mu[k]) + sigmaS * dx^2 / 2 / mu[k]^2);
            if dtTmp < dt
                dt = dtTmp;
                #println(w[k], " ", mu[k])
                #println(dt)
            end
        end
        dt = cfl*dt;

        # build class
        new(Nx,NCells,a,b,dx,tEnd,dt,cfl,nPN,x,xMid,sigmaA,sigmaS,epsilon,r);
    end

end

function IC(obj::Settings,x,xi=0.0)
    y = zeros(size(x));
    x0 = 0.0
    s1 = 0.03
    s2 = s1^2
    floor = 1e-4
    x0 = 0.0
    for j = 1:length(y);
        y[j] = max(floor,1.0/(sqrt(2*pi)*s1) *exp(-((x[j]-x0)*(x[j]-x0))/2.0/s2))
    end
    return y;
end