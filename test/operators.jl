using NodalDG
using Test

@testset "Operator correctness" begin
    FT = Float64

    # From Matlab script: https://github.com/tcew/nodal-dg
    P_correct = Dict()
    x_correct = Dict()
    w_correct = Dict()
    V_correct = Dict()
    dP_correct = Dict()
    ∇V_correct = Dict()
    Dr_correct = Dict()
    N=0
    x_correct[N] = [-0]
    w_correct[N] = [2]
    P_correct[N] = [0.86603]
    V_correct[N] = [0.70711, 0.70711, 0.70711, 0.70711, 0.70711, 0.70711]
    dP_correct[N] = [0 0 0 0 0 0]
    ∇V_correct[N] = [0 0 0 0 0 0]
    Dr_correct[N] = [-0   0   0   0   0   0;
                     -0   0   0   0   0   0;
                     -0   0   0   0   0   0;
                     -0   0   0   0   0   0;
                     -0   0   0   0   0   0;
                     -0   0   0   0   0   0]
    N=1
    x_correct[N] = [-0.44721, 0.44721]
    w_correct[N] = [0.66667, 0.66667]
    P_correct[N] = [-0.86603, 0.86603]
    V_correct[N] = [0.70711  -1.22474;
                    0.70711  -0.73485;
                    0.70711  -0.24495;
                    0.70711   0.24495;
                    0.70711   0.73485;
                    0.70711   1.22474]
    dP_correct[N] = [1.9365, 1.9365, 1.9365, 1.9365, 1.9365, 1.9365]
    ∇V_correct[N] = [0.00000   1.22474;
                     0.00000   1.22474;
                     0.00000   1.22474;
                     0.00000   1.22474;
                     0.00000   1.22474;
                     0.00000   1.22474]

    Dr_correct[N] = [-0.357143  -0.214286  -0.071429   0.071429   0.214286   0.357143;
                     -0.357143  -0.214286  -0.071429   0.071429   0.214286   0.357143;
                     -0.357143  -0.214286  -0.071429   0.071429   0.214286   0.357143;
                     -0.357143  -0.214286  -0.071429   0.071429   0.214286   0.357143;
                     -0.357143  -0.214286  -0.071429   0.071429   0.214286   0.357143;
                     -0.357143  -0.214286  -0.071429   0.071429   0.214286   0.357143]
    N=5
    x_correct[N] = [-0.87174,-0.59170,-0.20930, 0.20930, 0.59170, 0.87174]
    w_correct[N] = [0.050584,0.221693,0.394391,0.394391,0.221693,0.050584]
    P_correct[N] = [-1.13374, 0.89104, -0.81033, 0.81033, -0.89104, 1.13374]
    V_correct[N] = [0.707107  -1.224745   1.581139  -1.870829   2.121320  -2.345208;
                    0.707107  -0.734847   0.063246   0.673498  -0.865499   0.357973;
                    0.707107  -0.244949  -0.695701   0.523832   0.492146  -0.721198;
                    0.707107   0.244949  -0.695701  -0.523832   0.492146   0.721198;
                    0.707107   0.734847   0.063246  -0.673498  -0.865499  -0.357973;
                    0.707107   1.224745   1.581139   1.870829   2.121320   2.345208]
    dP_correct[N] = [82.6136,-6.2125, 1.7184, 1.7184,-6.2125, 82.6136]
    ∇V_correct[N] = [0.00000    1.22474   -4.74342   11.22497  -21.21320   35.17812;
                     0.00000    1.22474   -2.84605    2.24499    1.52735   -5.79735;
                     0.00000    1.22474   -0.94868   -2.24499    2.88500    2.08254;
                     0.00000    1.22474    0.94868   -2.24499   -2.88500    2.08254;
                     0.00000    1.22474    2.84605    2.24499   -1.52735   -5.79735;
                     0.00000    1.22474    4.74342   11.22497   21.21320   35.17812]

    Dr_correct[N] = [-5.708333   12.500000  -12.500000    8.333333   -3.125000    0.500000;
                     -0.500000   -2.708333    5.000000   -2.500000    0.833333   -0.125000;
                      0.125000   -1.250000   -0.833333    2.500000   -0.625000    0.083333;
                     -0.083333    0.625000   -2.500000    0.833333    1.250000   -0.125000;
                      0.125000   -0.833333    2.500000   -5.000000    2.708333    0.500000;
                     -0.500000    3.125000   -8.333333   12.500000  -12.500000    5.708333]



    α=FT(1); β=FT(1);
    r = collect(range(FT(-1), FT(1), length=6));

    for N in (0, 1, 5)
        (x,w) = JacobiGQ(α, β, N)
        @test all(isapprox.(x, x_correct[N], rtol=0.1))
        # The order of returned eigen values is
        # different between implementations, so
        # we sort them for the test.
        @test all(isapprox.(sort(w), sort(w_correct[N]), rtol=0.1))

        P = JacobiP(x,α,β,N)
        @test all(isapprox.(P, P_correct[N], rtol=0.1))

        V = Vandermonde1D(N,r)
        @test all(isapprox.(V, V_correct[N], rtol=0.1))

        Dr = Dmatrix1D(N,r,V)
        @test all(isapprox.(Dr, Dr_correct[N], rtol=0.1))

        dP = GradJacobiP(r, α, β, N)
        @test all(isapprox.(dP, dP_correct[N], rtol=0.1))

        ∇V = GradVandermonde1D(N,r)
        @test all(isapprox.(∇V, ∇V_correct[N], rtol=0.1))

    end

end
