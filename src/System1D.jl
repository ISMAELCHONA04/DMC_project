__precompile__(false)

"""
    System1D

Package module for one-dimensional quantum Monte Carlo experiments.

The public API is organized around a small domain layer (`Hamiltonian`,
boundaries, walkers, trial states) and three method layers (DMC, VMC, and
fixed-population GFMC).
"""
module System1D

using Random

# ============================================================================
# Domain Layer: Core types and interfaces
# ============================================================================

include("Domain/boundaries1D.jl")
include("Domain/hamiltonian.jl")
include("Domain/walker.jl")
include("Domain/nodepolicy.jl")
include("Domain/trialwf.jl")
include("Domain/guiding.jl")
include("Domain/spinless_fermion_ring_1d.jl")

# ============================================================================
# UseCases Layer: DMC and VMC implementations
# ============================================================================

# Common utilities shared by DMC and VMC
include("UseCases/common/walkers.jl")
include("UseCases/common/debug.jl")
include("UseCases/common/energy.jl")
include("UseCases/common/record.jl")

# Fixed-population GFMC implementation
include("UseCases/gfmc/abstract.jl")
include("UseCases/gfmc/params.jl")
include("UseCases/gfmc/reconfiguration.jl")
include("UseCases/gfmc/kernels.jl")
include("UseCases/gfmc/sim.jl")
include("UseCases/gfmc/run.jl")

# DMC implementation
include("UseCases/dmc/params.jl")
include("UseCases/dmc/sim.jl")
include("UseCases/dmc/moves.jl")
include("UseCases/dmc/branching.jl")
include("UseCases/dmc/et_update.jl")
include("UseCases/dmc/run.jl")
include("UseCases/dmc/plot.jl")

# VMC implementation
include("UseCases/vmc/params.jl")
include("UseCases/vmc/sim.jl")
include("UseCases/vmc/proposals.jl")
include("UseCases/vmc/estimators.jl")
include("UseCases/vmc/run.jl")

# Cross-method helpers
include("UseCases/common/warm_start.jl")

# ============================================================================
# Exports: Public API
# ============================================================================

# Domain layer exports
export AbstractBoundary1D, OpenBoundary1D, PeriodicBoundary1D
export isperiodic, cell_length, cell_bounds, in_cell
export wrap_coordinate, wrap_position!, wrap_position, wrap_position_with_phase
export minimum_image, displacement, distance_1d, squared_distance_1d
export Hamiltonian, potential, nparticles, diffusion_constant, boundary
export Walker, position, phase, clone, initialize!
export AbstractGuiding, NoGuiding, ImportanceGuiding, TrialWF
export AbstractNodePolicy, NoNode, FixedNode
export drift, signpsi, local_energy
export SpinlessFermionRing1D, trial_logpsi, trial_gradlogpsi, trial_lapllogpsi
export trial_wavefunction, importance_guiding, hamiltonian, sample_uniform_configurations

# GFMC exports
export AbstractGFMCKernel, GenericGFMCKernel, SpinlessFermionRing1DKernel
export AbstractGFMCBackend, CPUBackend
export AbstractGFMCReconfiguration, MultinomialReconfiguration, SystematicReconfiguration
export GFMCParams, GFMCSim, run_gfmc!

# DMC exports
export DMCParams, DMCSim
export step!, run_simulation!, run_dmc, run_dmc!
export propose_move, crosses_node, branching_factor
export update_ET!
export record_state!, record_positions!, estimate_energy, estimate_energy_variance
export plot_snapshot_1d_density, plot_snapshot_1d_points, plot_snapshot_1d_smoothed_density
export mean_walker_position_1d, sliding_window_average, animate_dmc_walker_distribution_1d

# VMC exports
export VMCParams, VMCSim
export AbstractVMCProposal, DriftGaussianProposal, GaussianProposal
export vmc_step!, run_vmc, run_vmc!
export metropolis_step!, proposal_step!
export estimate_energy_vmc, compute_variance

# Cross-method exports
export sample_vmc_positions, run_gfmc_with_vmc_init

end # module System1D
