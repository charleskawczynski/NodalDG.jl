using NodalDG
using LinearAlgebra
using Test
using Plots

@testset "Derivative" begin
    FT = Float64

    K = 10
    N_all = []
    x_all = []
    u_all = []
    u′_all = []
    ∇u_all = []
    for N in 2:6:20

        grid = Grid1D(FT;N=N,K=K, xmax=FT(2))
        push!(x_all, grid.x)

        u = sin.(π .* grid.x)
        push!(u_all, u)
        u′ = cos.(π * grid.x)
        push!(u′_all, u′)
        ∇u = grid.rx .* (grid.Dr*u)
        push!(∇u_all, ∇u)
        push!(N_all, N)
    end

    plot()
    for (i,N) in enumerate(N_all)
        x_a = [x_all[i]...]; y = [u_all[i]...]
        plot!(x_a,y, label="N=$N")
    end
    savefig(joinpath(fig_dir, "funapprox_K=$K.png"))

    plot()
    i = length(N_all)
    x_a = [x_all[i]...]; y = [u_all[i]...];plot!(x_a,y, label="u")
    x_a = [x_all[i]...]; y = [u′_all[i]...];plot!(x_a,y, label="u′")
    for (i,N) in enumerate(N_all)
        x_a = [x_all[i]...]; y = [∇u_all[i]...]; plot!(x_a,y, label="∇u, N=$N")
    end
    savefig(joinpath(fig_dir, "fun_derivative_K=$K.png"))

    @test 1==1
end
