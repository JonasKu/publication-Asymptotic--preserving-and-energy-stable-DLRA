include("settings.jl")
include("Solver.jl")

using DelimitedFiles
using NPZ
using PyPlot

close("all")

s = Settings();

############################
solver = Solver(s);
@time tEnd, rhoDiffusion = SolveDiffusion(solver);
@time tEnd, rho, u = SolveIMEX(solver);
@time tEnd, rhoDLRA, uDLRA = SolveDLRA(solver);

s.tEnd = tEnd;

# read reference solution
v = readdlm("PlaneSourceRaw", ',')
uEx = zeros(length(v));
for i = 1:length(v)
    if v[i] == ""
        uEx[i] = 0.0;
    else
        uEx[i] = Float64(v[i])
    end
end
x = collect(range(-1.5,1.5,length=(2*length(v)-1)));
uEx = [uEx[end:-1:2];uEx]

fig, ax = subplots()
ax.plot(x,uEx, "k-", linewidth=2, label="exact", alpha=0.6)
#ax.plot(s.xMid,rhoDiffusion, "g--", linewidth=2, label="diffusion", alpha=0.6)
ax.plot(s.xMid,rho, "b-.", linewidth=2, label="full, N=$(s.nPN)", alpha=0.6)
ax.plot(s.xMid,rhoDLRA, "r--", linewidth=2, label="DLRA, r=$(s.r)", alpha=0.6)
ax.legend(loc="upper left")
#ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.set_xlabel(L"x",fontsize=13)
ax.set_ylabel(L"\Phi",fontsize=13)
ax.tick_params("both",labelsize=13) 
show()

npzwrite("data/SolutionUnconventionalNx$(s.Nx)N$(s.nPN)tEnd$(s.tEnd)r$(s.r).npy", u)

println("main finished")
