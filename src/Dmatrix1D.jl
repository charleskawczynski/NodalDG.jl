export Dmatrix1D

"""
    Dmatrix1D(N,r,V)

Purpose : Initialize the (r) differentiation matrices on the interval,
        evaluated at (r) at order N
"""
function Dmatrix1D(N,r,V)
    Vr = GradVandermonde1D(N, r)
    Dr = Vr/V
    return Dr
end