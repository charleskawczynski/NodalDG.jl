export Connect1D
using SparseArrays
using LinearAlgebra

subv2indv(s, i) = LinearIndices(s)[CartesianIndex.(i)]

"""
    Connect1D(EToV)

Purpose  : Build global connectivity arrays for 1D grid based on standard
         EToV input array from grid generator
"""
function Connect1D(EToV)

    Nfaces = 2;
    # Find number of elements and vertices
    K = size(EToV,1)
    TotalFaces = Nfaces*K
    Nv = K+1;

    # List of local face to local vertex connections
    vn = [1,2];

    # Build global face to node sparse array
    SpFToV = spzeros(TotalFaces, Nv);
    sk = 1;
    for k=1:K
        for face=1:Nfaces
           SpFToV[sk, EToV[k, vn[face]]] = 1
           sk = sk+1;
        end
    end

    # Build global face to global face sparse array
    SpFToF = SpFToV*SpFToV' - I;

    # Find complete face to face connections
    faces1, faces2, vals = findnz(SpFToF .== 1) # maybe findall?

    # Convert face global number to element and face numbers
    element1 = floor.( (faces1 .- 1) ./ Nfaces ) .+ 1
    face1    =   mod.( (faces1 .- 1)  , Nfaces ) .+ 1
    element2 = floor.( (faces2 .- 1) ./ Nfaces ) .+ 1
    face2    =   mod.( (faces2 .- 1)  , Nfaces ) .+ 1

    element1 = convert.(Int, element1)
    element2 = convert.(Int, element2)
    # # Rearrange into Nelements x Nfaces sized arrays
    ind = subv2indv((K, Nfaces), collect(zip(element1, face1)));
    EToE      = (1:K)*ones(Int,1,Nfaces);
    EToF      = ones(Int,K,1)*(1:Nfaces)';
    EToE[ind] = element2
    EToF[ind] = face2

    return (EToE, EToF)

end
