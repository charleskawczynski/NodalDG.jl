module NodalDG

include("JacobiP.jl")
include("JacobiGQ.jl")
include("JacobiGL.jl")
include("Vandermonde1D.jl")
include("GradJacobiP.jl")
include("Dmatrix1D.jl")
include("GradVandermonde1D.jl")
include("Lift1D.jl")
include("GeometricFactors1D.jl")
include("Normals1D.jl")
include("Connect1D.jl")
include("MeshGen1D.jl")
include("BuildMaps1D.jl")

include("Grids1D.jl") # for encapsulation
include("TimeSteppers.jl")

end # module
