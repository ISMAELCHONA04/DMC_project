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
function DMCSim(h::H, params::DMCParams, walkers::Vector{W}, rng::RNG;
    guiding::G=NoGuiding(),
    nodepolicy::NP=NoNode(),
    step::Int=0) where {H,W,RNG,G<:AbstractGuiding,NP<:AbstractNodePolicy}
    return DMCSim{H,W,RNG,G,NP}(h, params, walkers, guiding, nodepolicy, rng, step, params.ET0, Float64[], Int[], Float64[], Float64[], Vector{Vector{Vector{Float64}}}())
end

function DMCSim(h::H, params::DMCParams, walkers::Vector{W}, guiding::G, rng::RNG, step::Int=0) where {H,W,G<:AbstractGuiding,RNG}
    return DMCSim(h, params, walkers, rng; guiding=guiding, nodepolicy=NoNode(), step=step)
end

function DMCSim(h::H, params::DMCParams, walkers::Vector{W}, rng::RNG, step::Int=0) where {H,W,RNG}
    return DMCSim(h, params, walkers, rng; guiding=NoGuiding(), nodepolicy=NoNode(), step=step)
end

# Initialize simulation with given walkers.
function initialize!(sim::DMCSim, positions::Vector{<:AbstractVector})
    sim.walkers = [Walker(pos) for pos in positions]
    sim.step = 0
    sim.ET = sim.params.ET0
    sim.ET_history = Float64[]
    sim.population_history = Int[]
    sim.energy_mean_history = Float64[]
    sim.energy_variance_history = Float64[]
    sim.walker_positions_history = Vector{Vector{Vector{Float64}}}()
    return sim
end
