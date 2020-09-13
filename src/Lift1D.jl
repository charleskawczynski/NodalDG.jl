export Lift1D

"""
    Lift1D(Np, Nfaces, Nfp, V)

Purpose  : Compute surface integral term in DG formulation
"""
function Lift1D(Np, Nfaces, Nfp, V)

    FT = eltype(V)
    Emat = zeros(FT, Np,Nfaces*Nfp)

    # Define Emat
    Emat[1,1] = 1
    Emat[Np,2] = 1

    # inv(mass matrix)*\sₙ (Lᵢ,Lⱼ)_{edgeₙ}v0
    LIFT = V*(V'*Emat)
    return LIFT
end
