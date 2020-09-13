
using SpecialFunctions: gamma
using Test
using LinearAlgebra
export JacobiGQ

"""
    JacobiGQ(α,β,N)

Purpose: Compute the N'th order Gauss quadrature points, x,
         and weights, w, associated with the Jacobi
         polynomial, of type (α,β) > -1 ( <> -0.5).
"""
function JacobiGQ(α::FT,β::FT,N) where {FT}

    if N==0
        x = FT[-(α-β)/(α+β+2)]
        w = FT[2]
        return (x,w)
    end

    # Form symmetric matrix from recurrence.
    J = zeros(FT, N+1)
    h1 = 2*(0:N) .+ α .+ β
    d1 = diagm(0=>-1/2*(α^2-β^2)./(h1 .+ 2)./h1)
    t2 = sqrt.((1:N).*((1:N) .+α .+ β).*
        ((1:N) .+ α).*((1:N) .+ β)./(h1[1:N] .+ 1)./(h1[1:N] .+ 3))
    J = d1 .+ diagm(1=>2 ./ (h1[1:N] .+ 2) .* t2)

    α+β<10*eps(FT) && (J[1,1]=0)

    J = J + J'

    # Compute quadrature by eigenvalue solve
    F = eigen(J)
    D = F.values
    V = F.vectors

    x = D

    w = (V[1,:]).^2*2^(α+β+1)/(α+β+1)*gamma(α+1)*
        gamma(β+1)/gamma(α+β+1);

    return (x,w)
end
