
export GradJacobiP

"""
    GradJacobiP(r, α, β, N)

Purpose: Evaluate the derivative of the Jacobi polynomial of type (α,β)>-1,
       at points r for order N and returns dP[1:length(r))]
"""
function GradJacobiP(r, α, β, N)

    FT = eltype(r)
    dP = zeros(FT, length(r), 1)
    if N == 0
        dP .= 0
    else
        dP = sqrt(N*(N+α+β+1)) * JacobiP(r[:],α+1,β+1, N-1)
    end

    return dP

end
