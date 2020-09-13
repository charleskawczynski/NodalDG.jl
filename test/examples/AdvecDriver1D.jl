# Driver script for solving the 1D advection equations
using NodalDG
FT = Float64
using UnPack

# Generate grid
grid = Grid1D(FT;N=8,K=10, xmax=FT(2))

function AdvecRHS1D!(u, grid, t, params)
    a = params.a
    Nfp = get_Nfp(grid)
    Nfaces = get_Nfaces(grid)
    K = n_elems(grid)
    topology = grid.topology
    @unpack vmapM, vmapP, vmapI, mapI, mapO = topology
    @unpack nx, Dr, rx, Fscale, LIFT = grid

    # form field differences at faces
    α=1;
    du = zeros(FT, Nfp*Nfaces,K);
    du[:] .= (u[vmapM] .- u[vmapP]).*(a*nx[:] .- (1-α)*abs.(a*nx[:]))/2;
    # impose boundary condition at x=0
    uin = -sin(a*t);
    du[mapI] = (u[vmapI] .- uin).*(a*nx[mapI] .- (1-α)*abs.(a*nx[mapI]))/2;
    du[mapO] = 0;
    # compute right hand sides of the semi-discrete PDE
    rhsu = -a*rx.*(Dr*u) .+ LIFT*(Fscale .* du);
    return rhsu
end

struct Params{FT}
  "advection speed"
  a::FT
end

FinalTime = 10;
CFL = FT(0.75)
Δt   = CFL/(2*FT(π))*grid.Δxmin
Δt = FT(0.5)*Δt;
Nsteps = convert(Int, ceil(FinalTime/Δt));
Δt = FinalTime/Nsteps;

# advection speed
params = Params(2*FT(π))

# Set initial conditions
u = sin.(grid.x)

# Solve Problem
lserk = LSERK(AdvecRHS1D!, u; Δt=Δt, Nsteps=Nsteps)
solve!(lserk, u, grid, params)

using Plots
plot(grid.x, u, xlabel="x", ylabel="u at t=$FinalTime")
