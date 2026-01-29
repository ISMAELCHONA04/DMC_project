# UseCases/common: Snapshot recording and state tracking

"""
    record_positions!(sim)

Record current walker positions as a deep copy to freeze snapshot.
"""
function record_positions!(sim)
    snap = [copy(position(w)) for w in sim.walkers]   # deep copy each R
    push!(sim.walker_positions_history, snap)
    return sim
end

"""
    record_state!(sim, record_positions::Bool=true)

Record current state of the simulation for history tracking.
Includes population, energy mean, reference energy, and optionally walker positions.
"""
function record_state!(sim, record_positions::Bool=true)
    push!(sim.population_history, length(sim.walkers))
    push!(sim.energy_mean_history, estimate_energy(sim))
    push!(sim.ET_history, sim.ET)
    if record_positions
        record_positions!(sim)
    end
    return sim
end
