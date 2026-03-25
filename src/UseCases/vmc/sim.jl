# UseCases/vmc: VMC simulation container with guided importance sampling

abstract type AbstractVMCProposal end

struct DriftGaussianProposal <: AbstractVMCProposal end

struct GaussianProposal <: AbstractVMCProposal end

_resolve_vmc_proposal(p::AbstractVMCProposal) = p

function _resolve_vmc_proposal(p::Symbol)
    if p === :drift_gaussian || p === :drifted_gaussian || p === :drift
        return DriftGaussianProposal()
    elseif p === :gaussian || p === :position_gaussian
        return GaussianProposal()
    end
    throw(ArgumentError("Unknown VMC proposal mode: $(p). Use :drift_gaussian (default) or :gaussian."))
end

mutable struct VMCSim{H,W,RNG,TW,P<:AbstractVMCProposal}
    H::H
    params::VMCParams
    walkers::Vector{W}
    trial_wf::TW  # Trial wavefunction (which is also the guiding)
    guiding::ImportanceGuiding{TW,H}  # Derived from trial_wf and H
    proposal::P   # proposal kernel policy
    rng::RNG
    step::Int
    energy_history::Vector{Float64}
    energy_variance_history::Vector{Float64}
    acceptance_history::Vector{Float64}
    walker_positions_history::Vector{Vector{Vector{Float64}}}
    acceptance_rate::Float64
    acceptance_count::Int
end

"""
    VMCSim(h::H, params::VMCParams, configs, trial_wf::TW, rng::RNG; proposal=DriftGaussianProposal(), step=0)

Create a VMC simulation with guided importance sampling.
The trial wavefunction serves as both the variational ansatz and the guiding wavefunction.

# Arguments
- `h::H`: Hamiltonian
- `params::VMCParams`: VMC parameters
- `configs`: Initial walkers or raw configurations
- `trial_wf::TW`: Trial wavefunction (TrialWF type)
- `rng::RNG`: Random number generator
- `proposal`: Proposal kernel policy (`DriftGaussianProposal()` default, `GaussianProposal()`, or symbols `:drift_gaussian` / `:gaussian`)
"""
function VMCSim(h::H, params::VMCParams, configs, trial_wf::TW, rng::RNG;
    proposal::Union{AbstractVMCProposal,Symbol}=DriftGaussianProposal(),
    step::Integer=0) where {H,RNG,TW}
    walkers = _wrapped_walkers(h, configs)
    guiding = ImportanceGuiding(trial_wf, h)
    proposal_policy = _resolve_vmc_proposal(proposal)
    WW = eltype(walkers)
    return VMCSim{H,WW,RNG,TW,typeof(proposal_policy)}(
        h,
        params,
        walkers,
        trial_wf,
        guiding,
        proposal_policy,
        rng,
        Int(step),
        Float64[],
        Float64[],
        Float64[],
        Vector{Vector{Vector{Float64}}}(),
        0.0,
        0
    )
end

"""
    initialize!(sim::VMCSim, configs)

Initialize VMC simulation with walkers or raw configurations.
"""
function initialize!(sim::VMCSim, configs)
    sim.walkers = _wrapped_walkers(sim.H, configs)
    sim.step = 0
    sim.energy_history = Float64[]
    sim.energy_variance_history = Float64[]
    sim.acceptance_history = Float64[]
    sim.walker_positions_history = Vector{Vector{Vector{Float64}}}()
    sim.acceptance_rate = 0.0
    sim.acceptance_count = 0
    return sim
end
