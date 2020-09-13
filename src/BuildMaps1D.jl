export BuildMaps1D

"""
BuildMaps1D(K, Np, Nfp, Nfaces, Fmask)

Purpose: Connectivity and boundary tables for nodes given in the K # of elements,
         each with N+1 degrees of freedom.
"""
function BuildMaps1D(K, Np, Nfp, Nfaces, x, Fmask, EToE, EToF; NODETOL=1e-10)

    # number volume nodes consecutively
    nodeids = reshape(collect(1:K*Np), Np, K);
    vmapM   = zeros(Int,Nfp, Nfaces, K);
    vmapP   = zeros(Int,Nfp, Nfaces, K);

    for k1=1:K
        for f1=1:Nfaces
            # find index of face nodes with respect to volume node ordering
            vmapM[:,f1,k1] = nodeids[Fmask[:,f1], k1];
        end
    end

    for k1=1:K
        for f1=1:Nfaces
            # find neighbor
            k2 = EToE[k1,f1]
            f2 = EToF[k1,f1]

            # find volume node numbers of left and right nodes
            vidM = vmapM[:,f1,k1]
            vidP = vmapM[:,f2,k2]

            x1  = x[vidM]
            x2  = x[vidP]

            # Compute distance matrix
            D = (x1 - x2).^2
            if all(D .< NODETOL)
                vmapP[:,f1,k1] = vidP
            end
        end
    end

    vmapP = vmapP[:]
    vmapM = vmapM[:]

    # Create list of boundary nodes
    mapB = vmapP .== vmapM
    vmapB = vmapM[mapB]

    # Create specific left (inflow) and right (outflow) maps
    mapI = 1
    mapO = K*Nfaces
    vmapI = 1
    vmapO = K*Np
    return (vmapM, vmapP, vmapB, mapB)
end
