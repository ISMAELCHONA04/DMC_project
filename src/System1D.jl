__precompile__(false)
module System1D

using Random

export Hamiltonian, Walker, DMCParams, DMCSim

include("Domain/hamiltonian.jl")
include("Domain/walker.jl")
include("UseCases/params.jl")
include("UseCases/sim.jl")

end # module System1D
