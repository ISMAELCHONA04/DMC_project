# UseCases/vmc: VMC simulation container (skeleton)

mutable struct VMCSim{H,W,RNG,G<:AbstractGuiding}
    H::H
    params::VMCParams
    walkers::Vector{W}
    guiding::G
    rng::RNG
    step::Int
    energy_history::Vector{Float64}
    walker_positions_history::Vector{Vector{Vector{Float64}}}
end

function VMCSim(h::H, params::VMCParams, walkers::Vector{W}, rng::RNG;
    guiding::G=NoGuiding()) where {H,W,RNG,G<:AbstractGuiding}
    return VMCSim{H,W,RNG,G}(h, params, walkers, guiding, rng, 0, Float64[], Vector{Vector{Vector{Float64}}}())
end

"""
    initialize!(sim::VMCSim, positions::Vector{<:AbstractVector})

Initialize VMC simulation with given walker positions.
"""
function initialize!(sim::VMCSim, positions::Vector{<:AbstractVector})
    sim.walkers = [Walker(pos) for pos in positions]
    sim.step = 0
    sim.energy_history = Float64[]
    sim.walker_positions_history = Vector{Vector{Vector{Float64}}}()
    return sim
end
