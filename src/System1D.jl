__precompile__(false)
module System1D

using Random

# ============================================================================
# Domain Layer: Core types and interfaces
# ============================================================================

include("Domain/hamiltonian.jl")
include("Domain/walker.jl")
include("Domain/nodepolicy.jl")
include("Domain/trialwf.jl")
include("Domain/guiding.jl")

# ============================================================================
# UseCases Layer: DMC and VMC implementations
# ============================================================================

# Common utilities shared by DMC and VMC
include("UseCases/common/energy.jl")
include("UseCases/common/record.jl")

# DMC implementation
include("UseCases/dmc/params.jl")
include("UseCases/dmc/sim.jl")
include("UseCases/dmc/moves.jl")
include("UseCases/dmc/branching.jl")
include("UseCases/dmc/et_update.jl")
include("UseCases/dmc/run.jl")
include("UseCases/dmc/plot.jl")

# VMC skeleton (not yet implemented)
include("UseCases/vmc/params.jl")
include("UseCases/vmc/sim.jl")
include("UseCases/vmc/proposals.jl")
include("UseCases/vmc/estimators.jl")
include("UseCases/vmc/run.jl")

# ============================================================================
# Exports: Public API
# ============================================================================

# Domain layer exports
export Hamiltonian, potential, nparticles, diffusion_constant
export Walker, position, clone, initialize!
export AbstractGuiding, NoGuiding, ImportanceGuiding, TrialWF
export AbstractNodePolicy, NoNode, FixedNode
export drift, signpsi, local_energy

# DMC exports
export DMCParams, DMCSim
export step!, run_simulation!, run_dmc!
export propose_move, crosses_node, branching_factor
export update_ET!
export record_state!, record_positions!, estimate_energy, estimate_energy_variance
export plot_snapshot_1d_density, plot_snapshot_1d_points

# VMC exports
export VMCParams, VMCSim
export vmc_step!, run_vmc!
export metropolis_step!, proposal_step!
export estimate_energy_vmc, compute_variance

end # module System1D
