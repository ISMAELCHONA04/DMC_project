function nb_find_project_root(start::AbstractString=pwd())
    dir = abspath(start)
    while true
        has_project = isfile(joinpath(dir, "Project.toml"))
        has_src = isfile(joinpath(dir, "src", "System1D.jl"))
        if has_project && has_src
            return dir
        end
        parent = dirname(dir)
        parent == dir && error("Could not locate project root from $start")
        dir = parent
    end
end

function nb_paths(project_root::AbstractString, notebook_rel_dir::AbstractString)
    notebook_dir = normpath(joinpath(project_root, notebook_rel_dir))
    output_dir = normpath(joinpath(notebook_dir, "..", "outputs"))
    figures_dir = joinpath(output_dir, "figures")
    tables_dir = joinpath(output_dir, "tables")

    mkpath(output_dir)
    mkpath(figures_dir)
    mkpath(tables_dir)

    return (
        project_root=project_root,
        notebook_dir=notebook_dir,
        output_dir=output_dir,
        figures_dir=figures_dir,
        tables_dir=tables_dir,
    )
end

function nb_include_formatting(notebook_dir::AbstractString)
    formatting_path = joinpath(notebook_dir, "formatting.jl")
    if isfile(formatting_path)
        include(formatting_path)
        return formatting_path
    end
    return nothing
end

function nb_default_snapshot_steps(nsteps::Integer)
    nsteps_int = Int(nsteps)
    nsteps_int >= 1 || throw(ArgumentError("nsteps must be >= 1"))
    return sort(unique([
        0,
        max(1, fld(nsteps_int, 4)),
        max(1, fld(nsteps_int, 2)),
        max(1, fld(3 * nsteps_int, 4)),
        nsteps_int,
    ]))
end

function nb_mean_sem(values::AbstractVector{<:Real})
    finite_vals = Float64[v for v in values if isfinite(v)]
    isempty(finite_vals) && return (NaN, NaN)
    length(finite_vals) == 1 && return (finite_vals[1], NaN)
    return mean(finite_vals), std(finite_vals) / sqrt(length(finite_vals))
end

function nb_moving_average(values::AbstractVector{<:Real}, window::Integer)
    n = length(values)
    n == 0 && return Float64[]

    w = max(1, Int(window))
    half = fld(w, 2)
    out = Vector{Float64}(undef, n)

    for i in 1:n
        ilo = max(1, i - half)
        ihi = min(n, i + half)
        total = 0.0
        count = 0
        for j in ilo:ihi
            v = Float64(values[j])
            if isfinite(v)
                total += v
                count += 1
            end
        end
        out[i] = count == 0 ? NaN : total / count
    end

    return out
end

function nb_padded_limits(values::AbstractVector{<:Real}; pad_frac::Float64=0.08)
    finite_vals = Float64[v for v in values if isfinite(v)]
    isempty(finite_vals) && return (-1.0, 1.0)

    lo = minimum(finite_vals)
    hi = maximum(finite_vals)
    span = hi - lo
    if span <= eps(Float64)
        center = 0.5 * (lo + hi)
        pad = max(1.0e-6, 0.05 * (abs(center) + 1.0))
        return (center - pad, center + pad)
    end

    pad = pad_frac * span
    return (lo - pad, hi + pad)
end

function nb_snapshot_coordinate(snapshot, coord::Integer=1)
    idx = Int(coord)
    return Float64[R[idx] for R in snapshot]
end

function _nb_snapshot_from_positions(X::AbstractMatrix{<:Real})
    snap = Vector{Vector{Float64}}(undef, size(X, 2))
    @inbounds for j in axes(X, 2)
        snap[j] = collect(Float64, @view X[:, j])
    end
    return snap
end

function nb_last_snapshot(sim)
    if !isempty(sim.walker_positions_history)
        return sim.walker_positions_history[end]
    end
    if hasproperty(sim, :walkers)
        return [copy(position(w)) for w in sim.walkers]
    elseif hasproperty(sim, :positions)
        return _nb_snapshot_from_positions(sim.positions)
    end
    error("Could not determine how to extract a snapshot from $(typeof(sim)).")
end

function nb_all_coordinates(sim; coord::Integer=1)
    coords = Float64[]
    for snapshot in sim.walker_positions_history
        append!(coords, nb_snapshot_coordinate(snapshot, coord))
    end
    if isempty(coords)
        append!(coords, nb_snapshot_coordinate(nb_last_snapshot(sim), coord))
    end
    return coords
end

function nb_density_curve(
    xs::AbstractVector{<:Real};
    nbins::Integer=120,
    xmin::Union{Nothing,Real}=nothing,
    xmax::Union{Nothing,Real}=nothing,
    smoothing_window::Integer=1,
)
    isempty(xs) && error("Cannot build a density curve from an empty coordinate set.")
    nbins_int = max(2, Int(nbins))
    smooth_int = max(1, Int(smoothing_window))

    xvals = Float64[x for x in xs if isfinite(x)]
    isempty(xvals) && error("No finite coordinates available for density curve.")

    xmin_f = xmin === nothing ? minimum(xvals) : Float64(xmin)
    xmax_f = xmax === nothing ? maximum(xvals) : Float64(xmax)

    if xmax_f <= xmin_f
        center = 0.5 * (xmax_f + xmin_f)
        xmin_f = center - 1.0
        xmax_f = center + 1.0
    end

    width = (xmax_f - xmin_f) / nbins_int
    counts = zeros(Float64, nbins_int)

    for x in xvals
        idx = floor(Int, (x - xmin_f) / width) + 1
        idx = clamp(idx, 1, nbins_int)
        counts[idx] += 1.0
    end

    density = counts ./ (length(xvals) * width)
    smoothed = nb_moving_average(density, smooth_int)
    centers = Float64[xmin_f + (i - 0.5) * width for i in 1:nbins_int]
    return centers, smoothed
end

function nb_density_curve_from_snapshot(
    snapshot;
    coord::Integer=1,
    nbins::Integer=120,
    xmin::Union{Nothing,Real}=nothing,
    xmax::Union{Nothing,Real}=nothing,
    smoothing_window::Integer=1,
)
    xs = nb_snapshot_coordinate(snapshot, coord)
    return nb_density_curve(
        xs;
        nbins=nbins,
        xmin=xmin,
        xmax=xmax,
        smoothing_window=smoothing_window,
    )
end

function nb_periodic_kde_curve(
    xs::AbstractVector{<:Real};
    xmin::Real,
    xmax::Real,
    grid_points::Integer=320,
    bandwidth::Union{Nothing,Real}=nothing,
)
    isempty(xs) && error("Cannot build a periodic KDE from an empty coordinate set.")

    xlo = Float64(xmin)
    xhi = Float64(xmax)
    L = xhi - xlo
    L > 0 || throw(ArgumentError("xmax must be greater than xmin"))

    raw_vals = Float64[x for x in xs if isfinite(x)]
    isempty(raw_vals) && error("No finite coordinates available for periodic KDE.")

    xvals = Vector{Float64}(undef, length(raw_vals))
    @inbounds for i in eachindex(raw_vals)
        y = raw_vals[i] - xlo
        xvals[i] = xlo + (y - L * floor(y / L))
    end

    n = length(xvals)
    grid_n = max(32, Int(grid_points))

    h = if bandwidth === nothing
        σ = std(xvals)
        candidate = σ > 0 ? (0.9 * σ * n^(-1 / 5)) : (L / 40)
        clamp(candidate, L / 200, L / 8)
    else
        Float64(bandwidth)
    end
    h > 0 || throw(ArgumentError("bandwidth must be > 0"))

    dx = L / grid_n
    centers = Vector{Float64}(undef, grid_n)
    density = zeros(Float64, grid_n)
    norm = 1.0 / (n * h * sqrt(2π))

    @inbounds for i in 1:grid_n
        xi = xlo + (i - 0.5) * dx
        centers[i] = xi
        accum = 0.0
        for x in xvals
            d = x - xi
            d -= L * floor((d + 0.5 * L) / L)
            z = d / h
            accum += exp(-0.5 * z * z)
        end
        density[i] = norm * accum
    end

    return centers, density
end

function nb_csv_escape(value)
    str = replace(string(value), '"' => "\"\"")
    if occursin(',', str) || occursin('"', str) || occursin('\n', str)
        return "\"$str\""
    end
    return str
end

function nb_write_csv(path::AbstractString, rows::AbstractVector{<:NamedTuple})
    isempty(rows) && error("Cannot write CSV with no rows.")
    keys0 = collect(keys(first(rows)))

    open(path, "w") do io
        println(io, join(string.(keys0), ","))
        for row in rows
            vals = [nb_csv_escape(getproperty(row, key)) for key in keys0]
            println(io, join(vals, ","))
        end
    end

    return path
end

function nb_save_figure(fig, figures_dir::AbstractString, stem::AbstractString, suffix::AbstractString; enabled::Bool=false)
    enabled || return nothing
    path = joinpath(figures_dir, "$(stem)_$(suffix).png")
    savefig(fig, path)
    println("Saved figure: ", abspath(path))
    return path
end

function nb_dmc_rows(label::AbstractString, sim)
    n = length(sim.energy_mean_history)
    taus = (0:(n - 1)) .* sim.params.dt
    rows = NamedTuple[]

    for i in 1:n
        push!(rows, (
            label=label,
            step=i - 1,
            tau=taus[i],
            population=sim.population_history[i],
            ET=sim.ET_history[i],
            energy_mean=sim.energy_mean_history[i],
            energy_var=sim.energy_variance_history[i],
        ))
    end

    return rows
end

function nb_vmc_rows(label::AbstractString, sim)
    n = length(sim.energy_history)
    taus = (0:(n - 1)) .* sim.params.dt
    rows = NamedTuple[]

    for i in 1:n
        push!(rows, (
            label=label,
            step=i - 1,
            tau=taus[i],
            population=sim.params.targetN,
            acceptance_rate=sim.acceptance_history[i],
            energy_mean=sim.energy_history[i],
            energy_var=sim.energy_variance_history[i],
        ))
    end

    return rows
end

function nb_gfmc_rows(label::AbstractString, sim)
    n = length(sim.energy_mean_history)
    taus = (0:(n - 1)) .* sim.params.dt
    rows = NamedTuple[]

    for i in 1:n
        push!(rows, (
            label=label,
            step=i - 1,
            tau=taus[i],
            population=sim.population_history[i],
            ET=sim.ET_history[i],
            energy_mean=sim.energy_mean_history[i],
            energy_var=sim.energy_variance_history[i],
            mean_weight=sim.mean_weight_history[i],
            effective_population=sim.effective_population_history[i],
            acceptance_rate=sim.acceptance_history[i],
        ))
    end

    return rows
end

function nb_plot_dmc_history(
    sims::AbstractVector;
    labels::AbstractVector{<:AbstractString},
    colors::AbstractVector,
    title_prefix::AbstractString="DMC",
)
    length(sims) == length(labels) == length(colors) || throw(ArgumentError("sims, labels, and colors must have the same length"))

    p_et = plot(xlabel="imaginary time", ylabel="E_T", title="$title_prefix: reference energy")
    p_pop = plot(xlabel="imaginary time", ylabel="walkers", title="$title_prefix: population")
    p_mean = plot(xlabel="imaginary time", ylabel="mean local energy", title="$title_prefix: mean energy")
    p_var = plot(xlabel="imaginary time", ylabel="variance", title="$title_prefix: energy variance")

    for (sim, label, color) in zip(sims, labels, colors)
        taus = (0:(length(sim.energy_mean_history) - 1)) .* sim.params.dt
        plot!(p_et, taus, sim.ET_history; label=label, color=color, linewidth=2.2)
        plot!(p_pop, taus, sim.population_history; label=label, color=color, linewidth=2.2)
        plot!(p_mean, taus, sim.energy_mean_history; label=label, color=color, linewidth=2.2)
        plot!(p_var, taus, sim.energy_variance_history; label=label, color=color, linewidth=2.2)
    end

    return plot(p_et, p_pop, p_mean, p_var; layout=(2, 2), size=(1100, 800))
end

function nb_plot_vmc_history(
    sims::AbstractVector;
    labels::AbstractVector{<:AbstractString},
    colors::AbstractVector,
    title_prefix::AbstractString="VMC",
)
    length(sims) == length(labels) == length(colors) || throw(ArgumentError("sims, labels, and colors must have the same length"))

    p_mean = plot(xlabel="imaginary time", ylabel="mean local energy", title="$title_prefix: mean energy")
    p_var = plot(xlabel="imaginary time", ylabel="variance", title="$title_prefix: energy variance")
    p_acc = plot(xlabel="imaginary time", ylabel="acceptance rate", title="$title_prefix: acceptance")
    p_pop = plot(xlabel="imaginary time", ylabel="walkers", title="$title_prefix: population")

    for (sim, label, color) in zip(sims, labels, colors)
        taus = (0:(length(sim.energy_history) - 1)) .* sim.params.dt
        population = fill(sim.params.targetN, length(taus))
        plot!(p_mean, taus, sim.energy_history; label=label, color=color, linewidth=2.2)
        plot!(p_var, taus, sim.energy_variance_history; label=label, color=color, linewidth=2.2)
        plot!(p_acc, taus, sim.acceptance_history; label=label, color=color, linewidth=2.2)
        plot!(p_pop, taus, population; label=label, color=color, linewidth=2.2)
    end

    return plot(p_mean, p_var, p_acc, p_pop; layout=(2, 2), size=(1100, 800))
end

function nb_plot_gfmc_history(
    sims::AbstractVector;
    labels::AbstractVector{<:AbstractString},
    colors::AbstractVector,
    title_prefix::AbstractString="GFMC",
)
    length(sims) == length(labels) == length(colors) || throw(ArgumentError("sims, labels, and colors must have the same length"))

    p_et = plot(xlabel="imaginary time", ylabel="E_T", title="$title_prefix: reference energy")
    p_mean = plot(xlabel="imaginary time", ylabel="mean branch energy", title="$title_prefix: mean energy")
    p_var = plot(xlabel="imaginary time", ylabel="variance", title="$title_prefix: energy variance")
    p_weight = plot(xlabel="imaginary time", ylabel="mean weight", title="$title_prefix: mean weight")
    p_eff = plot(xlabel="imaginary time", ylabel="effective walkers", title="$title_prefix: effective population")
    p_acc = plot(xlabel="imaginary time", ylabel="acceptance rate", title="$title_prefix: acceptance")

    for (sim, label, color) in zip(sims, labels, colors)
        taus = (0:(length(sim.energy_mean_history) - 1)) .* sim.params.dt
        plot!(p_et, taus, sim.ET_history; label=label, color=color, linewidth=2.2)
        plot!(p_mean, taus, sim.energy_mean_history; label=label, color=color, linewidth=2.2)
        plot!(p_var, taus, sim.energy_variance_history; label=label, color=color, linewidth=2.2)
        plot!(p_weight, taus, sim.mean_weight_history; label=label, color=color, linewidth=2.2)
        plot!(p_eff, taus, sim.effective_population_history; label=label, color=color, linewidth=2.2)
        plot!(p_acc, taus, sim.acceptance_history; label=label, color=color, linewidth=2.2)
    end

    return plot(p_et, p_mean, p_var, p_weight, p_eff, p_acc; layout=(2, 3), size=(1400, 800))
end
