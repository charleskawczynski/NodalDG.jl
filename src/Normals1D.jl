export Normals1D

"""
    Normals1D(::Type{FT}, Nfp, Nfaces, K) where {FT}

Purpose : Compute outward pointing normals at elements faces
"""
function Normals1D(::Type{FT}, Nfp, Nfaces, K) where {FT}

    nx = zeros(FT, Nfp*Nfaces, K)
    # Define outward normals
    nx[1,:] .= -1
    nx[2,:] .= 1
    return nx
end
