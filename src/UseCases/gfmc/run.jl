# UseCases/gfmc: Main fixed-population GFMC loop

function _branch_weight(sim::GFMCSim, Eold::Real, Enew::Real)
    logw = -0.5 * sim.params.dt * ((Float64(Eold) + Float64(Enew)) - 2.0 * sim.ET)
    w = exp(clamp(logw, -700.0, 700.0))
    return min(w, sim.params.branch_cap)
end

function _update_reference_energy!(sim::GFMCSim, mean_weight::Real)
    window = sim.params.energy_window
    if isempty(sim.energy_mean_history)
        Eblock = estimate_energy(sim)
    else
        hi = lastindex(sim.energy_mean_history)
        lo = max(1, hi - window + 1)
        block = @view sim.energy_mean_history[lo:hi]
        Eblock = sum(block) / length(block)
    end

    dt_block = sim.params.dt * sim.params.reconfiguration_interval
    correction = sim.params.feedback == 0 ? 0.0 : (sim.params.feedback / dt_block) * log(max(Float64(mean_weight), eps(Float64)))
    sim.ET = Eblock - correction
    return sim
end

function _reconfigure!(sim::GFMCSim)
    indices, mean_weight = resample_indices(sim.reconfiguration, sim.rng, sim.weights)
    copy_columns!(sim.backend, sim.proposal_positions, sim.positions, indices)
    copy_batch_data!(sim.backend, sim.proposal_data, sim.state_data, indices)

    tmp_positions = sim.positions
    sim.positions = sim.proposal_positions
    sim.proposal_positions = tmp_positions

    tmp_data = sim.state_data
    sim.state_data = sim.proposal_data
    sim.proposal_data = tmp_data

    fill_ones!(sim.backend, sim.weights)
    return mean_weight
end

function step!(sim::GFMCSim)
    D = diffusion_constant(sim.kernel)
    dt = sim.params.dt
    sigma = sqrt(2.0 * D * dt)
    nwalkers = size(sim.positions, 2)
    bc = boundary(sim.kernel)

    if uses_importance_sampling(sim.kernel)
        randn!(sim.backend, sim.rng, sim.proposal_positions)
        @inbounds for j in 1:nwalkers
            for i in axes(sim.positions, 1)
                sim.proposal_positions[i, j] = sim.positions[i, j] + D * dt * sim.state_data.drift[i, j] + sigma * sim.proposal_positions[i, j]
            end
        end
    else
        randn!(sim.backend, sim.rng, sim.proposal_positions)
        @inbounds for j in 1:nwalkers
            for i in axes(sim.positions, 1)
                sim.proposal_positions[i, j] = sim.positions[i, j] + sigma * sim.proposal_positions[i, j]
            end
        end
    end

    wrap_configurations!(sim.kernel, sim.proposal_positions)
    evaluate_configuration_data!(sim.kernel, sim.proposal_data, sim.proposal_positions)

    accept_count = 0
    @inbounds for j in 1:nwalkers
        killed = false

        if uses_fixed_node(sim.kernel)
            sold = sim.state_data.sign[j]
            snew = sim.proposal_data.sign[j]
            killed = (sold == 0.0) || (snew == 0.0) || (sold != snew)
        end

        if killed
            for i in axes(sim.positions, 1)
                sim.proposal_positions[i, j] = sim.positions[i, j]
                sim.proposal_data.drift[i, j] = sim.state_data.drift[i, j]
            end
            sim.proposal_data.logpsi[j] = sim.state_data.logpsi[j]
            sim.proposal_data.branch_energy[j] = sim.state_data.branch_energy[j]
            sim.proposal_data.sign[j] = sim.state_data.sign[j]
            sim.weights[j] = 0.0
            continue
        end

        if uses_importance_sampling(sim.kernel)
            denom = 4.0 * D * dt
            sq_f = 0.0
            sq_b = 0.0

            for i in axes(sim.positions, 1)
                mu_f = sim.positions[i, j] + D * dt * sim.state_data.drift[i, j]
                mu_b = sim.proposal_positions[i, j] + D * dt * sim.proposal_data.drift[i, j]
                df = displacement(bc, mu_f, sim.proposal_positions[i, j])
                db = displacement(bc, mu_b, sim.positions[i, j])
                sq_f += df * df
                sq_b += db * db
            end

            log_ratio = -(sq_b - sq_f) / denom + 2.0 * (sim.proposal_data.logpsi[j] - sim.state_data.logpsi[j])
            if log(rand(sim.rng)) >= min(0.0, log_ratio)
                for i in axes(sim.positions, 1)
                    sim.proposal_positions[i, j] = sim.positions[i, j]
                    sim.proposal_data.drift[i, j] = sim.state_data.drift[i, j]
                end
                sim.proposal_data.logpsi[j] = sim.state_data.logpsi[j]
                sim.proposal_data.branch_energy[j] = sim.state_data.branch_energy[j]
                sim.proposal_data.sign[j] = sim.state_data.sign[j]
            else
                accept_count += 1
            end
        else
            accept_count += 1
        end

        sim.weights[j] *= _branch_weight(sim, sim.state_data.branch_energy[j], sim.proposal_data.branch_energy[j])
    end

    tmp_positions = sim.positions
    sim.positions = sim.proposal_positions
    sim.proposal_positions = tmp_positions

    tmp_data = sim.state_data
    sim.state_data = sim.proposal_data
    sim.proposal_data = tmp_data

    push!(sim.mean_weight_history, sum(sim.weights) / length(sim.weights))
    push!(sim.effective_population_history, _effective_population(sim.weights))
    push!(sim.acceptance_history, nwalkers > 0 ? (accept_count / nwalkers) : 0.0)

    sim.step += 1

    if sim.step % sim.params.reconfiguration_interval == 0
        _reconfigure!(sim)
    end

    return sim
end

function run_simulation!(sim::GFMCSim;
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout)

    snapshot_set = isempty(snapshot_steps) ? nothing : Set(snapshot_steps)
    debug_every_int = Int(debug_every)
    debug_every_int >= 1 || throw(ArgumentError("debug_every must be >= 1, got $debug_every"))

    push!(sim.mean_weight_history, 1.0)
    push!(sim.effective_population_history, Float64(size(sim.positions, 2)))
    push!(sim.acceptance_history, 1.0)
    record_state!(sim, snapshot_set === nothing ? true : (0 in snapshot_set))

    nsteps = sim.params.nsteps
    stride = progress_every > 0 ? Int(progress_every) : max(1, fld(nsteps, 20))
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
            "mean_weight" => sim.mean_weight_history[end],
            "effective_population" => sim.effective_population_history[end],
            "acceptance_rate" => sim.acceptance_history[end],
            "energy_mean" => sim.energy_mean_history[end],
            "energy_var" => sim.energy_variance_history[end],
        )
    end

    for s in 1:nsteps
        step_t0 = time()
        step!(sim)

        if s > sim.params.nequil && s % sim.params.reconfiguration_interval == 0
            _update_reference_energy!(sim, sim.mean_weight_history[end])
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
                "mean_weight" => sim.mean_weight_history[end],
                "effective_population" => sim.effective_population_history[end],
                "acceptance_rate" => sim.acceptance_history[end],
                "reconfigured" => (sim.step % sim.params.reconfiguration_interval == 0),
                "energy_mean" => sim.energy_mean_history[end],
                "energy_var" => sim.energy_variance_history[end],
            )
        end
    end

    return sim
end

const run_gfmc! = run_simulation!
