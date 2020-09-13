# module TimeSteppers
export LSERK
export solve!

include("rk_coeffs.jl")

abstract type AbstractTimeStepper end

"""
    LSERK{F<:Function,AFT} <: AbstractTimeStepper

Low Storage Explicit Runge-Kutta
"""
struct LSERK{FT, AFT, F<:Function} <: AbstractTimeStepper
    "RHS eval"
    rhs!::F
    "residual storage"
    res::AFT
    "rhs storage"
    rhs::AFT
    "Time step"
    Δt::FT
    "Number of timesteps"
    Nsteps::Int
end

function LSERK(rhs!::F,
        u::AbstractArray{FT};
        Δt::FT,
        Nsteps::Int
    ) where {F<:Function, FT<:AbstractFloat}
    resu = similar(u); resu .= 0
    rhsu = similar(u); rhsu .= 0
    return LSERK{FT, typeof(resu), F}(rhs!, resu, rhsu, Δt, Nsteps)
end

function step!(
    ::Type{LSERK},
    grid,
    rhs!::RHS!,
    u::AFT,
    resu::AFT,
    rhsu::AFT,
    Δt::FT,
    t::AbstractArray{FT},
    params
    ) where {RHS!, FT, AFT<:AbstractArray{FT}}
    _rk4c = rk4c(FT)
    _rk4a = rk4a(FT)
    _rk4b = rk4b(FT)
    for INTRK = 1:5
        timelocal = t[1] + _rk4c[INTRK]*Δt
        rhsu .= rhs!(u, grid, timelocal, params)
        resu .= _rk4a[INTRK]*resu .+ Δt*rhsu
        u .= u .+ _rk4b[INTRK]*resu
    end
    return nothing
end

"""
    solve!(lserk::LSERK, u::AbstractArray{FT}, grid) where {FT}

"""
function solve!(lserk::LSERK, u::AbstractArray{FT}, grid, params) where {FT}
    t = FT[0];
    Δt = lserk.Δt
    resu = lserk.res
    rhsu = lserk.rhs
    rhs! = lserk.rhs!
    # outer time step loop
    for tstep=1:lserk.Nsteps
        step!(LSERK, grid, rhs!,u,resu,rhsu,lserk.Δt,t,params)
        # Increment time
        t[1] = t[1]+Δt
    end
end

# end