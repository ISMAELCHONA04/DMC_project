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
        proposal = _propose_move(sim, w, sim.guiding)
        Rnew = proposal.position

        if crosses_node(sim.nodepolicy, sim.guiding, Rold, Rnew)
            continue
        end

        P = branching_factor(sim, Rold, Rnew, sim.guiding)
        m = floor(Int, P + rand(sim.rng))
        propagated_phase = phase(w) * proposal.phase_factor

        for _ in 1:m
            push!(new_walkers, Walker(copy(Rnew), propagated_phase))  # copy is crucial for Vector R
        end
    end

    isempty(new_walkers) && error("No walkers left (extinction).")

    sim.walkers = new_walkers
    sim.step += 1
    return sim
end

"""
    run_simulation!(sim::DMCSim; snapshot_steps=Int[], show_progress=false, progress_every=0, progress_label="")

Run the DMC simulation for `sim.params.nsteps` steps, including an equilibration
phase of `sim.params.nequil` steps.

When `snapshot_steps` is non-empty, walker positions are recorded only at those
integer steps (including step 0). Otherwise, positions are recorded at every step.

When `show_progress=true`, prints progress updates in terms of `step/nsteps`.
If `progress_every <= 0`, progress is printed every ~5% of `nsteps`.

# Arguments
- `sim::DMCSim`: The simulation object to run
- `snapshot_steps::AbstractVector{<:Integer}`: Steps at which to record walker positions
- `show_progress::Bool`: Enable/disable textual progress output
- `progress_every::Int`: Progress print interval in steps (`<= 0` means auto)
- `progress_label::AbstractString`: Optional label included in progress lines
- `debug::Bool`: Emit step-level diagnostics every `debug_every` steps
- `debug_every::Int`: Debug print cadence in steps
- `debug_io::IO`: Output stream used for debug lines

# Returns
- `sim::DMCSim`: The simulation object (modified in-place)
"""
function run_simulation!(sim::DMCSim;
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout)

    snapshot_set = isempty(snapshot_steps) ? nothing : Set(snapshot_steps)
    record_state!(sim, snapshot_set === nothing ? true : (0 in snapshot_set))

    nsteps = sim.params.nsteps
    nsteps >= 1 || throw(ArgumentError("sim.params.nsteps must be >= 1, got $nsteps"))
    debug_every_int = Int(debug_every)
    debug_every_int >= 1 || throw(ArgumentError("debug_every must be >= 1, got $debug_every"))

    stride = progress_every > 0 ? Int(progress_every) : max(1, fld(nsteps, 20))  # ~5% cadence
    next_mark = stride
    t0 = time()
    label_prefix = isempty(progress_label) ? "[progress]" : "[progress][$progress_label]"
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
            "population" => sim.population_history[end],
            "ET" => sim.ET_history[end],
            "energy_mean" => sim.energy_mean_history[end],
            "energy_var" => sim.energy_variance_history[end],
        )
    end

    for s in 1:sim.params.nsteps
        step_t0 = time()
        step!(sim)

        if s > sim.params.nequil
            update_ET!(sim)   # update reference energy after equilibration
        end
        record_state!(sim, snapshot_set === nothing ? true : (s in snapshot_set))
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
                "population" => sim.population_history[end],
                "ET" => sim.ET_history[end],
                "energy_mean" => sim.energy_mean_history[end],
                "energy_var" => sim.energy_variance_history[end],
            )
        end
    end
    return sim
end

"""
    run_dmc(h, params, configs; kwargs...) -> DMCSim

Construct a `DMCSim` from walkers or raw positions, run it, and return the
simulation object.
"""
function run_dmc(h, params::DMCParams, configs;
    rng::AbstractRNG=Random.default_rng(),
    guiding::AbstractGuiding=NoGuiding(),
    nodepolicy::AbstractNodePolicy=NoNode(),
    step::Integer=0,
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout)
    sim = DMCSim(h, params, configs, rng; guiding=guiding, nodepolicy=nodepolicy, step=step)
    return run_simulation!(sim;
        snapshot_steps=snapshot_steps,
        show_progress=show_progress,
        progress_every=progress_every,
        progress_label=progress_label,
        debug=debug,
        debug_every=debug_every,
        debug_io=debug_io,
    )
end

# Alias for compatibility
const run_dmc! = run_simulation!
