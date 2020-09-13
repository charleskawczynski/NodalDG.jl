export GeometricFactors1D

"""
    GeometricFactors1D(x,Dr)

Purpose  : Compute the metric elements for the local mappings of the 1D elements
"""
function GeometricFactors1D(x,Dr)
    xr = Dr*x
    J = xr
    rx = 1 ./ J
    return (rx, J)
end
