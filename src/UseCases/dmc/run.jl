# UseCases/dmc: Main DMC simulation loop

"""
    step!(sim::DMCSim)

Perform a single DMC step: propose moves for each walker, check for node crossings,
calculate branching factors, and update walker population.
"""
function step!(sim::DMCSim)
    new_walkers = similar(sim.walkers, 0)  # same element type, empty

    for w in sim.walkers
        Rold = position(w)
        Rnew = propose_move(sim, w, sim.guiding)

        if crosses_node(sim.nodepolicy, sim.guiding, Rold, Rnew)
            continue
        end

        P = branching_factor(sim, Rold, Rnew, sim.guiding)
        m = floor(Int, P + rand(sim.rng))

        for _ in 1:m
            push!(new_walkers, Walker(copy(Rnew)))  # copy is crucial for Vector R
        end
    end

    isempty(new_walkers) && error("No walkers left (extinction).")

    sim.walkers = new_walkers
    sim.step += 1
    return sim
end

"""
    run_simulation!(sim::DMCSim; snapshot_steps=Int[])

Run the DMC simulation for `sim.params.nsteps` steps, including an equilibration
phase of `sim.params.nequil` steps.

When `snapshot_steps` is non-empty, walker positions are recorded only at those
integer steps (including step 0). Otherwise, positions are recorded at every step.

# Arguments
- `sim::DMCSim`: The simulation object to run
- `snapshot_steps::AbstractVector{<:Integer}`: Steps at which to record walker positions

# Returns
- `sim::DMCSim`: The simulation object (modified in-place)
"""
function run_simulation!(sim::DMCSim; snapshot_steps::AbstractVector{<:Integer}=Int[])
    snapshot_set = isempty(snapshot_steps) ? nothing : Set(snapshot_steps)
    record_state!(sim, snapshot_set === nothing ? true : (0 in snapshot_set))

    for s in 1:sim.params.nsteps
        step!(sim)

        if s > sim.params.nequil
            update_ET!(sim)   # update reference energy after equilibration
        end
        record_state!(sim, snapshot_set === nothing ? true : (s in snapshot_set))
    end
    return sim
end

# Alias for compatibility
const run_dmc! = run_simulation!
