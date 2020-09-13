export GradVandermonde1D
"""
    GradVandermonde1D(N,r)

Purpose : Initialize the gradient of the modal basis (i) at (r) at order N
"""
function GradVandermonde1D(N,r)

    FT = eltype(r)
    DVr = zeros(length(r),(N+1))

    # Initialize matrix
    for i=0:N
       DVr[:,i+1] = GradJacobiP(r[:],FT(0),FT(0),i)
    end
    return DVr

end
