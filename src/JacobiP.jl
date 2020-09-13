using SpecialFunctions: gamma

export JacobiP

"""
    JacobiP(x,α,β,N)

Purpose: Evaluate Jacobi Polynomial of type (α,β) > -1
         (α+β <> -1) at points x for order N and returns P[1:length(xp))]
Note   : They are normalized to be orthonormal.
"""
function JacobiP(x,α,β,N)

    # Turn points into row if needed.
    xp = x
    size(xp,2)==1 && (xp = xp')

    PL = zeros(N+1,length(xp))

    # Initial values P_0(x) and P_1(x)
    γ₀ = 2^(α+β+1)/(α+β+1)*gamma(α+1)*
        gamma(β+1)/gamma(α+β+1)

    PL[1,:] .= 1/sqrt(γ₀);
    N==0 && return PL

    γ₁ = (α+1)*(β+1)/(α+β+3)*γ₀;
    PL[2,:] .= ((α+β+2) * xp'/2 .+ (α-β)/2) ./ sqrt(γ₁);
    N==1 && return PL[N+1,:]

    # Repeat value in recurrence.
    a_old = 2/(2+α+β)*sqrt((α+1)*(β+1)/(α+β+3));

    # Forward recurrence using the symmetry of the recurrence.
    for i=1:N-1
        h₁ = 2*i+α+β
        a_new = 2/(h₁+2)*sqrt( (i+1)*(i+1+α+β)*(i+1+α)*
            (i+1+β)/(h₁+1)/(h₁+3))
        b_new = - (α^2-β^2)/h₁/(h₁+2)
        rhs = 1/a_new .* ( -a_old .* PL[i,:] .+ (xp .- b_new)' .* PL[i+1,:])
        PL[i+2,:] = rhs
        a_old = a_new
    end

    return PL[N+1,:]

end
