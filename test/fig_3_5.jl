using NodalDG
using LinearAlgebra
using Test
using Plots

@testset "Figure 3.5" begin
    FT = Float64
    α = FT(0)
    β = FT(0)
    # Derivative operator

    xmin,xmax,K = FT(0), FT(1), 20
    (Nv, VX, K, EToV) = MeshGen1D(xmin,xmax,K)

    N_all = []
    r_all = []
    x_all = []
    u_all = []
    u′_all = []
    ∇u_all = []
    Δu′_all = []
    for N in 2:6:20

        r = JacobiGL(α, β, N)
        push!(r_all, r)

        va = EToV[:,1]';
        vb = EToV[:,2]';
        x = ones(FT,N+1,1) .* VX[va] .+ FT(0.5)*(r .+ 1) .* (VX[vb] .- VX[va])
        push!(x_all, x)

        dP = GradJacobiP(r, α, β, N)
        ∇V = GradVandermonde1D(N,r)
        V = Vandermonde1D(N,r)
        @show size(V)
        Dr = Dmatrix1D(N,r,V)
        (rx,J) = GeometricFactors1D(x,Dr)
        u = sin.(π .* x)
        push!(u_all, u)
        u′ = cos.(π * x)
        push!(u′_all, u′)
        ∇u = rx .* (Dr*u)
        push!(∇u_all, ∇u)

        Δu′ = norm(abs.(u′ .- ∇u))
        push!(N_all, N)
        push!(Δu′_all, Δu′)
    end

    plot()
    for (i,N) in enumerate(N_all)
        x_a = [x_all[i]...]; y = [u_all[i]...]
        plot!(x_a,y, label="N=$N")
    end
    savefig("funapprox_K=$K.png")

    plot(N_all,Δu′_all; xlabel="N", ylabel="derivative error")
    savefig("derivative_error_K=$K.png")

    plot()
    i = length(N_all)
    x_a = [x_all[i]...]; y = [u_all[i]...];plot!(x_a,y, label="u")
    x_a = [x_all[i]...]; y = [u′_all[i]...];plot!(x_a,y, label="u′")
    for (i,N) in enumerate(N_all)
        x_a = [x_all[i]...]; y = [∇u_all[i]...]; plot!(x_a,y, label="∇u, N=$N")
    end
    savefig("fun_derivative_K=$K.png")

    @test 1==1
end
