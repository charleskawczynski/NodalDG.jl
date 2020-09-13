# module Grids1D

# using ...NodalDG

abstract type AbstractTopology1D end

struct Topology1D{AI2,AI1,AB1} <: AbstractTopology1D
    "Elem to vertex connectivity"
    EToV::AI2
    "Elem to elem connectivity"
    EToE::AI2
    "Elem to face connectivity"
    EToF::AI2
    "vertex maps"
    vmapM::AI1
    vmapP::AI1
    vmapB::AI1
    mapB::AB1
    "inflow / outflow maps"
    mapI::Int
    mapO::Int
    vmapI::Int
    vmapO::Int
end

abstract type AbstractGrid1D end
export Grid1D
struct Grid1D{N,K,FT,AF1,AF2,TT} <: AbstractGrid1D
    "Polynomial parameter"
    α::FT
    "Polynomial parameter"
    β::FT
    "Node tolerance"
    NODETOL::FT
    "Smallest nodal distance"
    Δxmin::FT
    "Node coordinates"
    VX::AF2
    "Physical coordinates"
    x::AF2
    "Reference coordinates"
    r::AF1
    "Normals"
    nx::AF2
    "Derivative matrix wrt reference elem"
    Dr::AF2
    "Inverse of Vandermonde matrix"
    V⁻¹::AF2
    "Fscale"
    Fscale::AF2
    "LIFT"
    LIFT::AF2
    "Metric terms"
    rx::AF2
    "Jacobian"
    J::AF2
    "Topology (connectivity between elems/vertices)"
    topology::TT
end

n_vertices(::Grid1D{N,K}) where {N,K} = K+1
export n_vertices
n_elems(::Grid1D{N,K}) where {N,K} = K
export n_elems
poly_order(::Grid1D{N}) where {N} = N
export poly_order
n_faces(::Grid1D) = 2
export n_faces
get_Nfp(::Grid1D) = 1
export get_Nfp
get_Nfaces(::Grid1D) = 2
export get_Nfaces

function Grid1D(::Type{FT};
    N::Int,
    K::Int,
    α::FT=FT(0),
    β::FT=FT(0),
    xmin::FT=FT(0),
    xmax::FT=FT(1),
    NODETOL::FT= FT(1e-10)
    ) where {FT<:AbstractFloat}
    (Nv, VX, K, EToV) = MeshGen1D(xmin,xmax,K)
    r = JacobiGL(α, β, N)
    va = EToV[:,1]'
    vb = EToV[:,2]'
    x = ones(FT,N+1,1) .* VX[va] .+ FT(0.5)*(r .+ 1) .* (VX[vb] .- VX[va])

    V = Vandermonde1D(N,r)
    V⁻¹ = inv(V)
    Dr = Dmatrix1D(N,r,V)
    (rx,J) = GeometricFactors1D(x,Dr)
    (EToE, EToF) = Connect1D(EToV)

    fmask1 = findall((abs.(r .+ 1) .< NODETOL))'
    fmask2 = findall((abs.(r .- 1) .< NODETOL))'

    # Fmask  = [fmask1;fmask2]' # open issue, cannot print this for N=1...
    Fmask = [fmask1;fmask2]
    Fmask = permutedims(Fmask)

    Np = N+1
    Nfaces = 2
    Nfp = 1
    Δxmin = min(abs.(x[1,:] .- x[2,:])...)
    nx = Normals1D(FT, Nfp, Nfaces, K)

    (vmapM, vmapP, vmapB, mapB) = BuildMaps1D(
      K, Np, Nfp, Nfaces, x, Fmask, EToE, EToF)
    AI2,AI1,AB1 = typeof(EToV), typeof(vmapM), typeof(mapB)

    # Create specific left (inflow) and right (outflow) maps
    mapI = 1
    mapO = K*Nfaces
    vmapI = 1
    vmapO = K*Np

    topology = Topology1D{AI2,AI1,AB1}(
        EToV, EToE, EToF,
        vmapM, vmapP, vmapB, mapB,
        mapI,mapO,vmapI,vmapO)

    AF1 = typeof(r)
    AF2 = typeof(Dr)
    TT = typeof(topology)

    Fscale = 1 ./ J[Fmask[:],:]
    LIFT = Lift1D(Np, Nfaces, Nfp, V)

    return Grid1D{N,K,FT,AF1,AF2,TT}(α, β, NODETOL,
      Δxmin, VX, x, r, nx, Dr, V⁻¹, Fscale, LIFT, rx, J, topology)
end

# end