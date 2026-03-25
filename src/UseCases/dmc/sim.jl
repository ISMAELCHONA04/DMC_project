# UseCases/dmc: DMC simulation container and state

mutable struct DMCSim{H,W,RNG,G<:AbstractGuiding,NP<:AbstractNodePolicy}
    H::H
    params::DMCParams
    walkers::Vector{W}
    guiding::G # guiding wavefunction / policy
    nodepolicy::NP # nodal surface policy
    rng::RNG
    step::Int
    ET::Float64 # Reference energy
    ET_history::Vector{Float64}
    population_history::Vector{Int}
    energy_mean_history::Vector{Float64}
    energy_variance_history::Vector{Float64}
    walker_positions_history::Vector{Vector{Vector{Float64}}}
end

# Outer constructor with sensible defaults for histories & ET.
function DMCSim(h::H, params::DMCParams, configs, rng::RNG;
    guiding::G=NoGuiding(),
    nodepolicy::NP=NoNode(),
    step::Integer=0) where {H,RNG,G<:AbstractGuiding,NP<:AbstractNodePolicy}
    wrapped_walkers = _wrapped_walkers(h, configs)
    WW = eltype(wrapped_walkers)
    return DMCSim{H,WW,RNG,G,NP}(h, params, wrapped_walkers, guiding, nodepolicy, rng, Int(step), params.ET0, Float64[], Int[], Float64[], Float64[], Vector{Vector{Vector{Float64}}}())
end

function DMCSim(h::H, params::DMCParams, configs, guiding::G, rng::RNG, step::Integer=0) where {H,G<:AbstractGuiding,RNG}
    return DMCSim(h, params, configs, rng; guiding=guiding, nodepolicy=NoNode(), step=step)
end

function DMCSim(h::H, params::DMCParams, configs, rng::RNG, step::Integer) where {H,RNG}
    return DMCSim(h, params, configs, rng; guiding=NoGuiding(), nodepolicy=NoNode(), step=step)
end

# Initialize simulation with given walkers or raw positions.
function initialize!(sim::DMCSim, configs)
    sim.walkers = _wrapped_walkers(sim.H, configs)
    sim.step = 0
    sim.ET = sim.params.ET0
    sim.ET_history = Float64[]
    sim.population_history = Int[]
    sim.energy_mean_history = Float64[]
    sim.energy_variance_history = Float64[]
    sim.walker_positions_history = Vector{Vector{Vector{Float64}}}()
    return sim
end
