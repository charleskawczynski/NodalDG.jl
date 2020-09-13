export JacobiGL

"""
    JacobiGL(α,β,N)

Purpose: Compute the N'th order Gauss Lobatto quadrature
         points, x, associated with the Jacobi polynomial,
         of type (α,β) > -1 ( <> -0.5).
"""
function JacobiGL(α::FT,β::FT,N) where {FT}

    if N==1
        x = zeros(FT, N+1)
        x[1]= FT(-1)
        x[2]= FT(1)
        return x
    end

    (xint,w) = JacobiGQ(α+1,β+1,N-2)

    return FT[-1, xint..., 1]
end
