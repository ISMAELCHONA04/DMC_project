# UseCases/dmc: Plotting utilities for DMC visualization

import Plots
import Plots: histogram, scatter, plot, plot!, vline!, hline!, @animate, gif

# Extract first coordinate for 1D plotting utilities.
_snapshot_xs_1d(snap::AbstractVector{<:AbstractVector}) = Float64[R[1] for R in snap]

"""
    plot_snapshot_1d_density(snap::AbstractVector{<:AbstractVector}; kwargs...)

Plot walker density as a histogram for a 1D snapshot.

# Arguments
- `snap`: Vector of walker positions (each position is a vector with first element x)
- `nbins::Int=200`: Number of histogram bins
- `xmin`, `xmax`: Custom axis limits (auto-detected if not provided)
- `ylims=nothing`: Optional fixed y-axis limits
- `normalize::Bool=true`: Whether to normalize as probability density
- `title::AbstractString`: Plot title
- `curve_label=false`: Legend label (or `false` to hide)
"""
function plot_snapshot_1d_density(snap::AbstractVector{<:AbstractVector};
    nbins::Int=200,
    xmin=nothing,
    xmax=nothing,
    ylims=nothing,
    normalize::Bool=true,
    title::AbstractString="Walker density",
    curve_label=false)
    xs = _snapshot_xs_1d(snap)
    isempty(xs) && error("Cannot plot density for an empty snapshot.")

    if xmin === nothing || xmax === nothing
        lo, hi = minimum(xs), maximum(xs)
        pad = 0.05 * (hi - lo + eps())   # small padding; eps avoids zero width
        xmin = (xmin === nothing) ? lo - pad : xmin
        xmax = (xmax === nothing) ? hi + pad : xmax
    end

    # `normalize=:pdf` makes it a density
    norm = normalize ? :pdf : :none
    return histogram(xs;
        bins=nbins,
        normalize=norm,
        xlabel="x",
        ylabel=normalize ? "density" : "count",
        title=title,
        label=curve_label,
        xlims=(xmin, xmax),
        ylims=ylims
    )
end

"""
    plot_snapshot_1d_points(snap::AbstractVector{<:AbstractVector}; kwargs...)

Plot walker positions as a scatter plot for a 1D snapshot.

# Arguments
- `snap`: Vector of walker positions (each position is a vector with first element x)
- `title::AbstractString`: Plot title
"""
function plot_snapshot_1d_points(snap::AbstractVector{<:AbstractVector};
    title::AbstractString="Walker positions")
    xs = _snapshot_xs_1d(snap)
    return scatter(xs, zeros(length(xs));
        xlabel="x",
        yticks=false,
        ylabel="",
        title=title,
        label=false
    )
end

function _centered_sliding_average(values::AbstractVector{<:Real}, window::Int)::Vector{Float64}
    window >= 1 || throw(ArgumentError("window must be >= 1"))
    n = length(values)
    out = Vector{Float64}(undef, n)
    half = fld(window, 2)

    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        total = 0.0
        for j in lo:hi
            total += Float64(values[j])
        end
        out[i] = total / (hi - lo + 1)
    end
    return out
end

function _snapshot_density_bins_1d(snap::AbstractVector{<:AbstractVector};
    nbins::Int=120,
    xmin=nothing,
    xmax=nothing,
    normalize::Bool=true)
    nbins >= 2 || throw(ArgumentError("nbins must be >= 2"))
    xs = _snapshot_xs_1d(snap)
    isempty(xs) && error("Cannot compute density for an empty snapshot.")

    if xmin === nothing || xmax === nothing
        lo, hi = minimum(xs), maximum(xs)
        pad = 0.05 * (hi - lo + eps())
        xmin = (xmin === nothing) ? lo - pad : xmin
        xmax = (xmax === nothing) ? hi + pad : xmax
    end
    xmax > xmin || throw(ArgumentError("xmax must be greater than xmin"))

    width = (xmax - xmin) / nbins
    counts = zeros(Float64, nbins)
    for x in xs
        idx = floor(Int, (x - xmin) / width) + 1
        idx = clamp(idx, 1, nbins)
        counts[idx] += 1.0
    end

    vals = if normalize
        counts ./ (length(xs) * width)
    else
        counts
    end
    centers = Float64[xmin + (i - 0.5) * width for i in 1:nbins]
    return centers, vals
end

"""
    harmonic_oscillator_ground_density_1d(x::Real; omega::Float64=1.0) -> Float64

Exact 1D harmonic-oscillator ground-state probability density `|psi(x)|^2`
in units where `m = 1` and `hbar = 1`.
"""
function harmonic_oscillator_ground_density_1d(x::Real; omega::Float64=1.0)::Float64
    omega > 0 || throw(ArgumentError("omega must be > 0"))
    return sqrt(omega / pi) * exp(-omega * x^2)
end

function _snapshot_local_energy_stats(sim::DMCSim, snap::AbstractVector{<:AbstractVector})::Tuple{Float64,Float64}
    vals = if sim.guiding isa NoGuiding
        Float64[potential(sim.H, R) for R in snap]
    else
        Float64[local_energy(sim.guiding, R) for R in snap]
    end
    isempty(vals) && return 0.0, 0.0

    mean_E = sum(vals) / length(vals)
    var_E = sum((E - mean_E)^2 for E in vals) / length(vals)
    return mean_E, var_E
end

"""
    plot_snapshot_1d_smoothed_density(snap::AbstractVector{<:AbstractVector}; kwargs...)

Plot a smoothed 1D walker probability density using a sliding-window average over binned density values.

# Arguments
- `snap`: Vector of walker positions (each position is a vector with first element x)
- `nbins::Int=120`: Number of bins used to estimate the density
- `smoothing_window::Int=9`: Sliding-window size in bin units for smoothing
- `xmin`, `xmax`: Optional fixed x-axis limits
- `ylims=nothing`: Optional fixed y-axis limits
- `normalize::Bool=true`: Plot probability density when true, otherwise counts
- `title::AbstractString`: Plot title
- `curve_label::AbstractString="Smoothed density"`: Legend label for the sampled curve
"""
function plot_snapshot_1d_smoothed_density(snap::AbstractVector{<:AbstractVector};
    nbins::Int=120,
    smoothing_window::Int=9,
    xmin=nothing,
    xmax=nothing,
    ylims=nothing,
    normalize::Bool=true,
    title::AbstractString="Walker density (smoothed)",
    curve_label::AbstractString="Smoothed density")
    centers, vals = _snapshot_density_bins_1d(
        snap;
        nbins=nbins,
        xmin=xmin,
        xmax=xmax,
        normalize=normalize
    )
    smooth = _centered_sliding_average(vals, smoothing_window)

    return plot(
        centers,
        smooth;
        linewidth=2.2,
        color=:navy,
        xlabel="x",
        ylabel=normalize ? "density" : "count",
        title=title,
        label=curve_label,
        legend=:topright,
        xlims=(minimum(centers), maximum(centers)),
        ylims=ylims
    )
end

"""
    mean_walker_position_1d(snap::AbstractVector{<:AbstractVector}) -> Float64

Compute the mean walker position for a 1D snapshot.
"""
function mean_walker_position_1d(snap::AbstractVector{<:AbstractVector})::Float64
    xs = _snapshot_xs_1d(snap)
    isempty(xs) && return 0.0
    return sum(xs) / length(xs)
end

"""
    sliding_window_average(values::AbstractVector{<:Real}, window::Int) -> Vector{Float64}

Compute a causal sliding-window average.
For the first `window-1` points, uses the available prefix.
"""
function sliding_window_average(values::AbstractVector{<:Real}, window::Int)::Vector{Float64}
    window >= 1 || throw(ArgumentError("window must be >= 1"))
    n = length(values)
    out = Vector{Float64}(undef, n)
    running = 0.0

    for i in 1:n
        running += Float64(values[i])
        if i > window
            running -= Float64(values[i - window])
            out[i] = running / window
        else
            out[i] = running / i
        end
    end
    return out
end

function _x_limits_1d(snapshots::AbstractVector{<:AbstractVector{<:AbstractVector}})
    lo = Inf
    hi = -Inf
    for snap in snapshots
        xs = _snapshot_xs_1d(snap)
        isempty(xs) && continue
        xlo = minimum(xs)
        xhi = maximum(xs)
        lo = min(lo, xlo)
        hi = max(hi, xhi)
    end
    (isfinite(lo) && isfinite(hi)) || error("No walker positions available for x-limits.")

    if lo == hi
        pad = 1.0
    else
        pad = 0.05 * (hi - lo)
    end
    return lo - pad, hi + pad
end

function _validate_density_plot_style(style::Symbol)::Symbol
    if style in (:smoothed, :histogram)
        return style
    end
    throw(ArgumentError("density_plot_style must be :smoothed or :histogram"))
end

function _validated_axis_limits(lims, name::AbstractString)
    lims === nothing && return nothing
    try
        lo = Float64(first(lims))
        hi = Float64(last(lims))
        hi > lo || throw(ArgumentError("$name upper limit must be greater than lower limit"))
        return (lo, hi)
    catch
        throw(ArgumentError("$name must be `nothing` or a 2-element tuple/vector `(ymin, ymax)`"))
    end
end

function _padded_limits(lo::Real, hi::Real; pad_frac::Float64=0.08)::Tuple{Float64,Float64}
    lo_f = Float64(lo)
    hi_f = Float64(hi)
    if !isfinite(lo_f) || !isfinite(hi_f)
        return (-1.0, 1.0)
    end
    if hi_f < lo_f
        lo_f, hi_f = hi_f, lo_f
    end

    span = hi_f - lo_f
    if span <= eps(Float64)
        center = 0.5 * (lo_f + hi_f)
        pad = max(1e-6, 0.05 * (abs(center) + 1.0))
        return center - pad, center + pad
    end

    pad = pad_frac * span
    return lo_f - pad, hi_f + pad
end

"""
    animate_dmc_walker_distribution_1d(sim::DMCSim; kwargs...) -> AbstractString

Create a GIF animation where each frame is a three-panel plot:
1. Walker probability distribution for the current snapshot (smoothed line or histogram)
2. Mean local energy over snapshots
3. Local-energy variance over snapshots

Returns the output path of the generated GIF.

# Keyword arguments
- `output_path::AbstractString="dmc_walker_distribution.gif"`
- `fps::Int=20`
- `nbins::Int=120`
- `density_smoothing_window::Int=9`
- `density_plot_style::Symbol=:smoothed`: choose `:smoothed` or `:histogram` for the distribution panel
- `fix_axis_scales::Bool=true`: keep y-axis scales fixed across all frames
- `normalize::Bool=true`
- `xmin=nothing`, `xmax=nothing`: Optional fixed x-axis limits for the distribution plot
- `distribution_ylims=nothing`: manual y-axis limits for the top distribution panel
- `mean_energy_ylims=nothing`: manual y-axis limits for the mean local-energy panel
- `variance_ylims=nothing`: manual y-axis limits for the variance panel
- `analytic_density::Union{Nothing,Function}=x -> harmonic_oscillator_ground_density_1d(x)`:
  optional exact density to overlay on the sampled distribution
- `analytic_energy::Union{Nothing,Real}=0.5`: optional exact energy line for the diagnostics panel
- `size::Tuple{Int,Int}=(1000, 900)`
"""
function animate_dmc_walker_distribution_1d(sim::DMCSim;
    output_path::AbstractString="dmc_walker_distribution.gif",
    fps::Int=20,
    nbins::Int=120,
    density_smoothing_window::Int=9,
    density_plot_style::Symbol=:smoothed,
    fix_axis_scales::Bool=true,
    normalize::Bool=true,
    xmin=nothing,
    xmax=nothing,
    distribution_ylims=nothing,
    mean_energy_ylims=nothing,
    variance_ylims=nothing,
    analytic_density::Union{Nothing,Function}=x -> harmonic_oscillator_ground_density_1d(x),
    analytic_energy::Union{Nothing,Real}=0.5,
    size::Tuple{Int,Int}=(1000, 900))

    snapshots = sim.walker_positions_history
    isempty(snapshots) && error("No snapshots found. Run `run_simulation!` first to record walker positions.")
    density_style = _validate_density_plot_style(density_plot_style)

    xlo, xhi = if xmin === nothing || xmax === nothing
        _x_limits_1d(snapshots)
    else
        (xmin, xmax)
    end

    nframes = length(snapshots)
    mean_local_E = Vector{Float64}(undef, nframes)
    var_local_E = Vector{Float64}(undef, nframes)
    for i in 1:nframes
        mean_local_E[i], var_local_E[i] = _snapshot_local_energy_stats(sim, snapshots[i])
    end
    xsnaps = collect(0:(nframes - 1))
    xdiag_hi = max(1, nframes - 1)
    xref = range(xlo, xhi; length=400)
    yref_density = analytic_density === nothing ? nothing : Float64[analytic_density(x) for x in xref]
    dist_panel_ylims = _validated_axis_limits(distribution_ylims, "distribution_ylims")
    mean_panel_ylims = _validated_axis_limits(mean_energy_ylims, "mean_energy_ylims")
    var_panel_ylims = _validated_axis_limits(variance_ylims, "variance_ylims")

    if fix_axis_scales
        if dist_panel_ylims === nothing
            dist_lo = 0.0
            dist_hi = 0.0
            for snap in snapshots
                _, vals = _snapshot_density_bins_1d(
                    snap;
                    nbins=nbins,
                    xmin=xlo,
                    xmax=xhi,
                    normalize=normalize
                )
                ys = density_style == :smoothed ? _centered_sliding_average(vals, density_smoothing_window) : vals
                if !isempty(ys)
                    dist_lo = min(dist_lo, minimum(ys))
                    dist_hi = max(dist_hi, maximum(ys))
                end
            end
            if yref_density !== nothing
                dist_lo = min(dist_lo, minimum(yref_density))
                dist_hi = max(dist_hi, maximum(yref_density))
            end
            dist_lo = min(dist_lo, 0.0)
            dlo, dhi = _padded_limits(dist_lo, dist_hi)
            dist_panel_ylims = dist_lo >= 0.0 ? (0.0, dhi) : (dlo, dhi)
        end

        if mean_panel_ylims === nothing
            mean_lo = minimum(mean_local_E)
            mean_hi = maximum(mean_local_E)
            if analytic_energy !== nothing
                exact_E = Float64(analytic_energy)
                mean_lo = min(mean_lo, exact_E)
                mean_hi = max(mean_hi, exact_E)
            end
            mean_panel_ylims = _padded_limits(mean_lo, mean_hi)
        end

        if var_panel_ylims === nothing
            var_lo = min(0.0, minimum(var_local_E))
            var_hi = max(0.0, maximum(var_local_E))
            vlo, vhi = _padded_limits(var_lo, var_hi)
            var_panel_ylims = var_lo >= 0.0 ? (0.0, vhi) : (vlo, vhi)
        end
    end

    anim = Plots.@animate for i in 1:nframes
        snap = snapshots[i]
        p_dist = if density_style == :smoothed
            plot_snapshot_1d_smoothed_density(
                snap;
                nbins=nbins,
                smoothing_window=density_smoothing_window,
                xmin=xlo,
                xmax=xhi,
                ylims=dist_panel_ylims,
                normalize=normalize,
                title="Smoothed Walker Distribution (snapshot $(i - 1))",
                curve_label="Sampled distribution"
            )
        else
            plot_snapshot_1d_density(
                snap;
                nbins=nbins,
                xmin=xlo,
                xmax=xhi,
                ylims=dist_panel_ylims,
                normalize=normalize,
                title="Histogram Walker Distribution (snapshot $(i - 1))",
                curve_label="Sampled distribution"
            )
        end
        if yref_density !== nothing
            plot!(
                p_dist,
                xref,
                yref_density;
                color=:gray30,
                linestyle=:dash,
                linewidth=2.0,
                label="Exact |psi|^2",
                fillrange=zero(Float64),
                fillcolor=:gray70,
                fillalpha=0.12
            )
        end

        p_avg = plot(
            xsnaps[1:i],
            mean_local_E[1:i];
            label="Mean local E",
            linewidth=2.0,
            color=:steelblue,
            xlabel="Snapshot",
            ylabel="<E_L>",
            title="Mean Local Energy by Snapshot",
            xlims=(0, xdiag_hi),
            ylims=mean_panel_ylims
        )
        if analytic_energy !== nothing
            hline!(
                p_avg,
                [Float64(analytic_energy)];
                label="Exact analytic E0",
                color=:black,
                linestyle=:dash,
                linewidth=1.8
            )
        end
        vline!(p_avg, [xsnaps[i]]; label="Current", color=:gray, linestyle=:dash, alpha=0.65)

        p_var = plot(
            xsnaps[1:i],
            var_local_E[1:i];
            label="Var(local E)",
            linewidth=2.0,
            color=:darkorange,
            xlabel="Snapshot",
            ylabel="Var(E_L)",
            title="Local-Energy Variance by Snapshot",
            xlims=(0, xdiag_hi),
            ylims=var_panel_ylims
        )
        vline!(p_var, [xsnaps[i]]; label="Current", color=:gray, linestyle=:dash, alpha=0.65)

        plot(p_dist, p_avg, p_var; layout=(3, 1), size=size)
    end

    Plots.gif(anim, output_path; fps=fps)
    return output_path
end
