using NodalDG
using LinearAlgebra
using Test
using Plots

@testset "Derivative Error" begin
    FT = Float64

    K = 5
    N_all = []
    Δu′_all = []
    for N in 2:15
        grid = Grid1D(FT;N=N,K=K, xmax=FT(2))
        u = sin.(FT(π) .* grid.x)
        u′ = FT(π)*cos.(FT(π) * grid.x)
        ∇u = grid.rx .* (grid.Dr*u)

        Δu′ = norm(abs.(u′ .- ∇u))
        push!(N_all, N)
        push!(Δu′_all, Δu′)
    end

    plot(N_all,Δu′_all; xlabel="N", ylabel="derivative error", yaxis=:log)
    savefig(joinpath(fig_dir, "derivative_error_K=$K.png"))

    @test 1==1
end
