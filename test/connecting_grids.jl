using NodalDG
using Test

@testset "Connecting grids" begin
    FT = Float64
    K = 20
    N = 5
    α = FT(0)
    β = FT(0)
    NODETOL = 1e-10
    Nfaces = 2
    Nfp = 1
    Np = N+1

    r = JacobiGL(α, β, N)
    V = Vandermonde1D(N,r)
    Dr = Dmatrix1D(N,r,V)

    # Compute masks for edge nodes
    fmask1 = (abs.(r .+ 1) .< NODETOL)'
    fmask2 = (abs.(r .- 1) .< NODETOL)'
    Fmask  = [fmask1;fmask2]'

    (Nv, VX, K, EToV) = MeshGen1D(FT(0),FT(1),K)
    va = EToV[:,1]'
    vb = EToV[:,2]'
    x = ones(FT,N+1,1) .* VX[va] .+ FT(0.5)*(r .+ 1) .* (VX[vb] .- VX[va])
    (rx,J) = GeometricFactors1D(x,Dr)

    @test size(VX) == (1,K+1)
    @test size(EToV) == (K,2)
    @test Nv==K+1
    @test VX == [0.00000;0.05000;0.10000;0.15000;0.20000;
                 0.25000;0.30000;0.35000;0.40000;0.45000;
                 0.50000;0.55000;0.60000;0.65000;0.70000;
                 0.75000;0.80000;0.85000;0.90000;0.95000;
                 1.00000]'
    @test EToV == [1 2; 2 3; 3 4; 4 5; 5 6;
                   6 7; 7 8; 8 9; 9 10; 10 11;
                  11 12; 12 13; 13 14; 14 15; 15 16;
                  16 17; 17 18; 18 19; 19 20; 20 21]

    (EToE, EToF) = Connect1D(EToV)
    @test eltype(EToE) == Int
    @test eltype(EToF) == Int

    @test EToE == [1 2; 1 3; 2 4; 3 5; 4 6;
                   5 7; 6 8; 7 9; 8 10; 9 11;
                  10 12; 11 13; 12 14; 13 15; 14 16;
                  15 17; 16 18; 17 19; 18 20; 19 20]

    @test EToF == [1 1; 2 1; 2 1; 2 1; 2 1;
                   2 1; 2 1; 2 1; 2 1; 2 1;
                   2 1; 2 1; 2 1; 2 1; 2 1;
                   2 1; 2 1; 2 1; 2 1; 2 2]

   (vmapM, vmapP, vmapB, mapB) = BuildMaps1D(K, Np, Nfp, Nfaces, x, Fmask, EToE, EToF)

   # Test that we can instantiate an encapsulating struct:
   grid = Grid1D(FT;N=8,K=10, xmax=FT(2))
   grid = Grid1D(FT;N=1,K=1, xmax=FT(2)); # edge case grid

end