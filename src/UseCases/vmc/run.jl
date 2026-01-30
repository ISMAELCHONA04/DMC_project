# UseCases/vmc: Main VMC simulation loop with guided importance sampling

"""
    vmc_step!(sim::VMCSim)

Perform a single VMC step: apply Metropolis accept/reject moves to all walkers.
"""
function vmc_step!(sim::VMCSim)
    new_walkers = similar(sim.walkers, 0)  # same element type, empty
    
    for w in sim.walkers
        w_new = metropolis_step!(sim, w)
        push!(new_walkers, w_new)
    end
    
    sim.walkers = new_walkers
    sim.step += 1
    
    # Update acceptance rate
    sim.acceptance_rate = sim.acceptance_count / (sim.step * length(sim.walkers))
    
    return sim
end

"""
    run_vmc!(sim::VMCSim; snapshot_steps=Int[])

Run the VMC simulation for `sim.params.nsteps` steps using guided importance sampling.

# Arguments
- `sim::VMCSim`: The VMC simulation object
- `snapshot_steps::AbstractVector{<:Integer}`: Steps at which to record walker positions

# Returns
- `sim::VMCSim`: The simulation object (modified in-place)
"""
function run_vmc!(sim::VMCSim; snapshot_steps::AbstractVector{<:Integer}=Int[])
    snapshot_set = isempty(snapshot_steps) ? nothing : Set(snapshot_steps)
    
    # Record initial state
    push!(sim.energy_history, estimate_energy(sim))
    push!(sim.energy_variance_history, estimate_energy_variance(sim))
    if snapshot_set === nothing || (0 in snapshot_set)
        snap = [copy(position(w)) for w in sim.walkers]
        push!(sim.walker_positions_history, snap)
    end
    
    # Run VMC steps
    for s in 1:sim.params.nsteps
        vmc_step!(sim)
        
        push!(sim.energy_history, estimate_energy(sim))
        push!(sim.energy_variance_history, estimate_energy_variance(sim))
        if snapshot_set === nothing || (s in snapshot_set)
            snap = [copy(position(w)) for w in sim.walkers]
            push!(sim.walker_positions_history, snap)
        end
    end
    
    return sim
end
