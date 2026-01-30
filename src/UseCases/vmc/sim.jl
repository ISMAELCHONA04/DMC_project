# UseCases/vmc: VMC simulation container with guided importance sampling

mutable struct VMCSim{H,W,RNG,TW}
    H::H
    params::VMCParams
    walkers::Vector{W}
    trial_wf::TW  # Trial wavefunction (which is also the guiding)
    guiding::ImportanceGuiding{TW,H}  # Derived from trial_wf and H
    rng::RNG
    step::Int
    energy_history::Vector{Float64}
    energy_variance_history::Vector{Float64}
    walker_positions_history::Vector{Vector{Vector{Float64}}}
    acceptance_rate::Float64
    acceptance_count::Int
end

"""
    VMCSim(h::H, params::VMCParams, walkers::Vector{W}, trial_wf::TW, rng::RNG)

Create a VMC simulation with guided importance sampling.
The trial wavefunction serves as both the variational ansatz and the guiding wavefunction.

# Arguments
- `h::H`: Hamiltonian
- `params::VMCParams`: VMC parameters
- `walkers::Vector{W}`: Initial walker positions
- `trial_wf::TW`: Trial wavefunction (TrialWF type)
- `rng::RNG`: Random number generator
"""
function VMCSim(h::H, params::VMCParams, walkers::Vector{W}, trial_wf::TW, rng::RNG) where {H,W,RNG,TW}
    guiding = ImportanceGuiding(trial_wf, h)
    return VMCSim{H,W,RNG,TW}(h, params, walkers, trial_wf, guiding, rng, 0, Float64[], Float64[], Vector{Vector{Vector{Float64}}}(), 0.0, 0)
end

"""
    initialize!(sim::VMCSim, positions::Vector{<:AbstractVector})

Initialize VMC simulation with given walker positions.
"""
function initialize!(sim::VMCSim, positions::Vector{<:AbstractVector})
    sim.walkers = [Walker(pos) for pos in positions]
    sim.step = 0
    sim.energy_history = Float64[]
    sim.energy_variance_history = Float64[]
    sim.walker_positions_history = Vector{Vector{Vector{Float64}}}()
    sim.acceptance_rate = 0.0
    sim.acceptance_count = 0
    return sim
end
