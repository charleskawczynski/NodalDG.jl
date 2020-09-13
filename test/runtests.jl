using NodalDG
using Test
const ndg_dir = dirname(dirname(pathof(NodalDG)))
const out_dir = joinpath(ndg_dir, "output")
const fig_dir = joinpath(out_dir, "figs")
mkpath(fig_dir)

@testset "NodalDG" begin
    include("operators.jl")
    include("connecting_grids.jl")
    include("derivative.jl")
    include("derivative_error.jl")
end
