export Vandermonde1D

"""
    Vandermonde1D(N,r)

return [V1D]

Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);
"""
function Vandermonde1D(N,r)
    FT = eltype(r)
    V1D = zeros(FT,length(r),N+1)
    for j=1:N+1
        V1D[:,j] = JacobiP(r[:], FT(0), FT(0), j-1)
    end
    return V1D
end
