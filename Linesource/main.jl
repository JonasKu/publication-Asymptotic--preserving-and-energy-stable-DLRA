include("settings.jl")
include("Solver.jl")

using DelimitedFiles
using NPZ
using PyPlot

close("all")


N = 500;

############################ epsilon = 1 ############################
epsilon = 1
r1 = 40
nx = 2002;
sE1 = Settings(nx,N,r1,epsilon);
solverE1 = Solver(sE1);
@time tEnd, rho, u, energy = SolveIMEX(solverE1);
@time tEnd, rhoDLRA, uDLRA, energyDLRA = SolveDLRA(solverE1);

############################ epsilon = 1e-5 ############################
epsilon = 1e-5
r = 3
nx = 502;
s = Settings(nx,N,r,epsilon);
s.tEnd = 0.2
solver = Solver(s);
@time tEnd, rhoDiffusion, energyDiffusion = SolveDiffusion(solver);
@time tEnd, rhoD, uD, energyD = SolveIMEX(solver);
@time tEnd, rhoDDLRA, uDDLRA, energyDDLRA = SolveDLRA(solver);


############################ plotting ############################
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

## epsilon = 1 ##
fig, ax = subplots()
ax.plot(x,uEx, "k-", linewidth=2, label="exact", alpha=0.6)
ax.plot(sE1.x,rho, "b-.", linewidth=2, label="full", alpha=0.6)
ax.plot(sE1.x,rhoDLRA, "r--", linewidth=2, label="DLRA", alpha=0.6)
ax.legend(loc="upper left",fontsize=13)
#ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([sE1.a,sE1.b])
ax.set_xlabel(L"x",fontsize=13)
ax.set_ylabel(L"\rho",fontsize=13)
ax.tick_params("both",labelsize=13) 
show()

fig, ax = subplots()
ax.plot(energy[1,:],energy[2,:], "b-.", linewidth=2, label="full", alpha=0.6)
ax.plot(energyDLRA[1,:],energyDLRA[2,:], "r--", linewidth=2, label="DLRA", alpha=0.6)
ax.legend(loc="upper right",fontsize=13)
#ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([energy[1,1],energy[1,end]])
ax.set_xlabel(L"t",fontsize=13)
ax.set_ylabel(L"e",fontsize=13)
ax.tick_params("both",labelsize=13) 
show()

## epsilon = 1e-5 ##
fig, ax = subplots()
ax.plot(s.x,rhoDiffusion, "g-", linewidth=2, label="diffusion", alpha=0.6)
ax.plot(s.x,rhoD, "b-.", linewidth=2, label="full", alpha=0.6)
ax.plot(s.x,rhoDDLRA, "r--", linewidth=2, label="DLRA", alpha=0.6)
ax.legend(loc="upper left",fontsize=13)
#ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.set_xlabel(L"x",fontsize=13)
ax.set_ylabel(L"\rho",fontsize=13)
ax.tick_params("both",labelsize=13) 
show()

fig, ax = subplots()
ax.plot(energyDiffusion[1,:],energyDiffusion[2,:], "g-", linewidth=2, label="diffusion", alpha=0.6)
ax.plot(energyD[1,:],energyD[2,:], "b-.", linewidth=2, label="full", alpha=0.6)
ax.plot(energyDDLRA[1,:],energyDDLRA[2,:], "r--", linewidth=2, label="DLRA", alpha=0.6)
ax.legend(loc="upper right",fontsize=13)
#ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([energyD[1,1],energyD[1,end]])
ax.set_xlabel(L"t",fontsize=13)
ax.set_ylabel(L"e",fontsize=13)
ax.tick_params("both",labelsize=13) 
show()

println("main finished")
