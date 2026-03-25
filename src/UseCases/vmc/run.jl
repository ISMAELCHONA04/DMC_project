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
- `show_progress::Bool`: Enable or disable progress printing
- `progress_every::Int`: Progress print interval in steps (`<= 0` means auto)
- `progress_label::AbstractString`: Optional label appended to progress lines
- `debug::Bool`: Emit step-level diagnostics every `debug_every` steps
- `debug_every::Int`: Debug print cadence in steps
- `debug_io::IO`: Output stream used for debug lines

# Returns
- `sim::VMCSim`: The simulation object (modified in-place)
"""
function run_vmc!(sim::VMCSim;
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout)
    snapshot_set = isempty(snapshot_steps) ? nothing : Set(snapshot_steps)

    nsteps = sim.params.nsteps
    nsteps >= 1 || throw(ArgumentError("sim.params.nsteps must be >= 1, got $nsteps"))
    debug_every_int = Int(debug_every)
    debug_every_int >= 1 || throw(ArgumentError("debug_every must be >= 1, got $debug_every"))

    stride = progress_every > 0 ? Int(progress_every) : max(1, fld(nsteps, 20))
    next_mark = stride
    t0 = time()
    label_prefix = isempty(progress_label) ? "[progress]" : "[progress][$progress_label]"

    # Record initial state
    push!(sim.energy_history, estimate_energy(sim))
    push!(sim.energy_variance_history, estimate_energy_variance(sim))
    push!(sim.acceptance_history, sim.acceptance_rate)
    if snapshot_set === nothing || (0 in snapshot_set)
        snap = [copy(position(w)) for w in sim.walkers]
        push!(sim.walker_positions_history, snap)
    end

    if show_progress
        println(label_prefix, " step 0/", nsteps, " (0.0%, elapsed 0.0s)")
    end
    if debug
        _emit_step_debug(
            debug_io,
            progress_label,
            sim.step,
            nsteps,
            0.0,
            0.0,
            "population" => length(sim.walkers),
            "acceptance_rate" => sim.acceptance_rate,
            "energy_mean" => sim.energy_history[end],
            "energy_var" => sim.energy_variance_history[end],
        )
    end

    # Run VMC steps
    for s in 1:nsteps
        step_t0 = time()
        vmc_step!(sim)

        push!(sim.energy_history, estimate_energy(sim))
        push!(sim.energy_variance_history, estimate_energy_variance(sim))
        push!(sim.acceptance_history, sim.acceptance_rate)
        if snapshot_set === nothing || (s in snapshot_set)
            snap = [copy(position(w)) for w in sim.walkers]
            push!(sim.walker_positions_history, snap)
        end
        step_elapsed = time() - step_t0

        if show_progress && (s >= next_mark || s == nsteps)
            frac = round((100.0 * s) / nsteps; digits=1)
            elapsed = round(time() - t0; digits=1)
            println(label_prefix, " step ", s, "/", nsteps, " (", frac, "%, elapsed ", elapsed, "s)")
            next_mark += stride
        end

        if debug && (s % debug_every_int == 0 || s == nsteps)
            _emit_step_debug(
                debug_io,
                progress_label,
                sim.step,
                nsteps,
                step_elapsed,
                time() - t0,
                "population" => length(sim.walkers),
                "acceptance_rate" => sim.acceptance_rate,
                "accepted_total" => sim.acceptance_count,
                "energy_mean" => sim.energy_history[end],
                "energy_var" => sim.energy_variance_history[end],
            )
        end
    end

    return sim
end

"""
    run_vmc(h, params, configs, trial_wf; kwargs...) -> VMCSim

Construct a `VMCSim` from walkers or raw positions, run it, and return the
simulation object.
"""
function run_vmc(h, params::VMCParams, configs, trial_wf;
    rng::AbstractRNG=Random.default_rng(),
    proposal::Union{AbstractVMCProposal,Symbol}=DriftGaussianProposal(),
    step::Integer=0,
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout)
    sim = VMCSim(h, params, configs, trial_wf, rng; proposal=proposal, step=step)
    return run_vmc!(sim;
        snapshot_steps=snapshot_steps,
        show_progress=show_progress,
        progress_every=progress_every,
        progress_label=progress_label,
        debug=debug,
        debug_every=debug_every,
        debug_io=debug_io,
    )
end
