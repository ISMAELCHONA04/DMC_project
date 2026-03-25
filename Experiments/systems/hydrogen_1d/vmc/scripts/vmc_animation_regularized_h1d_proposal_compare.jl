using Random
using Statistics
using Plots
using Logging

include("../../../../../src/System1D.jl")
using .System1D

function main()

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
alpha = parse(Float64, get(ENV, "VMC_H1D_REG_ALPHA", "1.0"))
softening = parse(Float64, get(ENV, "VMC_H1D_REG_SOFTENING", "0.5"))

# Defaults are tuned for the repository's regularized 1D hydrogen VMC comparison.
targetN = parse(Int, get(ENV, "VMC_H1D_REG_TARGETN", "50000"))
dt = parse(Float64, get(ENV, "VMC_H1D_REG_DT", "0.01"))
nsteps = parse(Int, get(ENV, "VMC_H1D_REG_NSTEPS", "160"))

nbins = parse(Int, get(ENV, "VMC_H1D_REG_NBINS", "240"))
smoothing_window = parse(Int, get(ENV, "VMC_H1D_REG_SMOOTH_WINDOW", "9"))
fps = parse(Int, get(ENV, "VMC_H1D_REG_GIF_FPS", "5"))
output_dir = joinpath(dirname(@__DIR__), "outputs", "animations")
mkpath(output_dir)
output_file = joinpath(output_dir, get(ENV, "VMC_H1D_REG_GIF_NAME", "vmc_h1d_reg_proposal_compare.gif"))

# Expand animation x-limits at least to this range for stable overlays.
xmin_hint = parse(Float64, get(ENV, "VMC_H1D_REG_XMIN", "-10.0"))
xmax_hint = parse(Float64, get(ENV, "VMC_H1D_REG_XMAX", "10.0"))

# Snapshot schedule for distribution overlays.
snapshot_candidates = collect(0:nsteps)
snapshot_steps = sort(unique(filter(s -> 0 <= s <= nsteps, snapshot_candidates)))

# -----------------------------------------------------------------------------
# System + trial wavefunction
# -----------------------------------------------------------------------------
V(R) = -1 / sqrt(R[1]^2 + softening^2)
H = Hamiltonian(1, 0.5, V)

logpsi(R) = begin
    x = R[1]
    r = sqrt(x^2 + softening^2)
    return -alpha * r
end

gradlogpsi(R) = begin
    x = R[1]
    r = sqrt(x^2 + softening^2)
    return [-alpha * x / r]
end

lapllogpsi(R) = begin
    x = R[1]
    r = sqrt(x^2 + softening^2)
    return -alpha * softening^2 / (r^3)
end

trial = TrialWF(logpsi, gradlogpsi, lapllogpsi)
params = VMCParams(; dt=dt, nsteps=nsteps, targetN=targetN, ET0=-0.5)

# Shared initial positions for fair proposal comparison.
rng_init = MersenneTwister(1234)
base_positions = [[randn(rng_init)] for _ in 1:targetN]

# Importance-sampling (drifted Gaussian) proposal simulation.
sim_drift = run_vmc(
    H,
    params,
    base_positions,
    trial;
    rng=MersenneTwister(42),
    proposal=DriftGaussianProposal(),
    snapshot_steps=snapshot_steps,
)

# Pure Gaussian proposal simulation.
sim_gauss = run_vmc(
    H,
    params,
    base_positions,
    trial;
    rng=MersenneTwister(43),
    proposal=GaussianProposal(),
    snapshot_steps=snapshot_steps,
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
_snapshot_xs_1d(snap) = Float64[R[1] for R in snap]

function _smoothed_density_curve(
    snap;
    nbins::Int=240,
    smoothing_window::Int=9,
    xmin::Float64=-10.0,
    xmax::Float64=10.0
)
    nbins >= 2 || throw(ArgumentError("nbins must be >= 2"))
    smoothing_window >= 1 || throw(ArgumentError("smoothing_window must be >= 1"))
    xmax > xmin || throw(ArgumentError("xmax must be greater than xmin"))

    xs = _snapshot_xs_1d(snap)
    counts = zeros(Float64, nbins)
    width = (xmax - xmin) / nbins

    for x in xs
        idx = floor(Int, (x - xmin) / width) + 1
        idx = clamp(idx, 1, nbins)
        counts[idx] += 1.0
    end

    density = counts ./ (length(xs) * width)

    smooth = similar(density)
    half = fld(smoothing_window, 2)
    for i in eachindex(density)
        lo = max(1, i - half)
        hi = min(length(density), i + half)
        smooth[i] = sum(@view density[lo:hi]) / (hi - lo + 1)
    end

    centers = Float64[xmin + (i - 0.5) * width for i in 1:nbins]
    return centers, smooth
end

function _combined_x_limits(snapshots_a, snapshots_b)
    lo = Inf
    hi = -Inf

    for snap in snapshots_a
        xs = _snapshot_xs_1d(snap)
        isempty(xs) && continue
        lo = min(lo, minimum(xs))
        hi = max(hi, maximum(xs))
    end
    for snap in snapshots_b
        xs = _snapshot_xs_1d(snap)
        isempty(xs) && continue
        lo = min(lo, minimum(xs))
        hi = max(hi, maximum(xs))
    end

    if !isfinite(lo) || !isfinite(hi)
        return (-10.0, 10.0)
    end

    if lo == hi
        return (lo - 1.0, hi + 1.0)
    end

    pad = 0.08 * (hi - lo)
    return (lo - pad, hi + pad)
end

function _normalize_on_grid(xs, ys)
    length(xs) == length(ys) || throw(ArgumentError("xs and ys must have the same length"))
    length(xs) >= 2 || throw(ArgumentError("need at least 2 points to normalize"))

    area = 0.0
    for i in 1:(length(xs) - 1)
        dx = xs[i + 1] - xs[i]
        area += 0.5 * (ys[i] + ys[i + 1]) * dx
    end
    area > 0 || throw(ArgumentError("normalization area must be positive"))
    return ys ./ area
end

_padded_limits(lo::Real, hi::Real; pad_frac::Float64=0.08) = begin
    lo_f = Float64(lo)
    hi_f = Float64(hi)
    if hi_f < lo_f
        lo_f, hi_f = hi_f, lo_f
    end

    span = hi_f - lo_f
    if span <= eps(Float64)
        center = 0.5 * (lo_f + hi_f)
        pad = max(1e-6, 0.05 * (abs(center) + 1.0))
        return (center - pad, center + pad)
    end

    pad = pad_frac * span
    return (lo_f - pad, hi_f + pad)
end

# -----------------------------------------------------------------------------
# Build frame-wise data with robust snapshot alignment
# -----------------------------------------------------------------------------
available_steps_drift = [s for s in snapshot_steps if s <= sim_drift.step]
available_steps_gauss = [s for s in snapshot_steps if s <= sim_gauss.step]

nframes = min(
    length(available_steps_drift),
    length(available_steps_gauss),
    length(sim_drift.walker_positions_history),
    length(sim_gauss.walker_positions_history)
)
nframes > 0 || error("No snapshots available to animate.")

frame_steps = available_steps_drift[1:nframes]
if frame_steps != available_steps_gauss[1:nframes]
    @warn "Snapshot steps differ between simulations; using drift step labels."
end

snaps_drift = sim_drift.walker_positions_history[1:nframes]
snaps_gauss = sim_gauss.walker_positions_history[1:nframes]

xlo_raw, xhi_raw = _combined_x_limits(snaps_drift, snaps_gauss)
xlo = min(xlo_raw, xmin_hint)
xhi = max(xhi_raw, xmax_hint)

# Target sampled density from trial wavefunction: |psi_T|^2
target_density_unnorm(x) = exp(-2 * alpha * sqrt(x^2 + softening^2))
xref = collect(range(xlo, xhi; length=1200))
yref_unnorm = Float64[target_density_unnorm(x) for x in xref]
yref = _normalize_on_grid(xref, yref_unnorm)

# Precompute smoothed densities for stable y-limits and faster frames.
curves_drift = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, nframes)
curves_gauss = Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef, nframes)

dist_max = maximum(yref)
for i in 1:nframes
    curves_drift[i] = _smoothed_density_curve(
        snaps_drift[i];
        nbins=nbins,
        smoothing_window=smoothing_window,
        xmin=xlo,
        xmax=xhi
    )
    curves_gauss[i] = _smoothed_density_curve(
        snaps_gauss[i];
        nbins=nbins,
        smoothing_window=smoothing_window,
        xmin=xlo,
        xmax=xhi
    )

    dist_max = max(dist_max, maximum(curves_drift[i][2]))
    dist_max = max(dist_max, maximum(curves_gauss[i][2]))
end
dist_ylims = (0.0, 1.08 * dist_max)

# Diagnostics limits shared across both proposal types.
E_drift = sim_drift.energy_history
E_gauss = sim_gauss.energy_history
V_drift = sim_drift.energy_variance_history
V_gauss = sim_gauss.energy_variance_history

energy_ylims = _padded_limits(minimum(vcat(E_drift, E_gauss)), maximum(vcat(E_drift, E_gauss)))
var_lo = min(0.0, minimum(vcat(V_drift, V_gauss)))
var_hi = max(0.0, maximum(vcat(V_drift, V_gauss)))
var_ylims = if var_lo >= 0.0
    (0.0, _padded_limits(var_lo, var_hi)[2])
else
    _padded_limits(var_lo, var_hi)
end

# -----------------------------------------------------------------------------
# Animate: 3 rows x 2 cols (left importance sampling, right gaussian proposal)
# -----------------------------------------------------------------------------
anim = @animate for i in 1:nframes
    step = frame_steps[i]
    t = step * dt

    centers_d, smoothed_d = curves_drift[i]
    centers_g, smoothed_g = curves_gauss[i]

    p_dist_d = plot(
        centers_d,
        smoothed_d;
        color=:navy,
        linewidth=2.2,
        xlabel="x",
        ylabel="density",
        title="Importance Sampling proposal | t=$(round(t; digits=2)) (step $(step))",
        label="sampled",
        xlims=(xlo, xhi),
        ylims=dist_ylims,
        legend=:topright
    )
    plot!(p_dist_d, xref, yref; color=:black, linestyle=:dash, linewidth=2.0, label="target |psi_T|^2")

    p_dist_g = plot(
        centers_g,
        smoothed_g;
        color=:darkorange,
        linewidth=2.2,
        xlabel="x",
        ylabel="density",
        title="Gaussian proposal | t=$(round(t; digits=2)) (step $(step))",
        label="sampled",
        xlims=(xlo, xhi),
        ylims=dist_ylims,
        legend=:topright
    )
    plot!(p_dist_g, xref, yref; color=:black, linestyle=:dash, linewidth=2.0, label="target |psi_T|^2")

    idx = min(step + 1, length(E_drift))
    t_hist = (0:(idx - 1)) .* dt

    p_mean_d = plot(
        t_hist,
        E_drift[1:idx];
        color=:steelblue,
        linewidth=2.0,
        xlabel="time",
        ylabel="mean local E",
        title="Importance Sampling proposal: mean local energy",
        label="mean E",
        xlims=(0, nsteps * dt),
        ylims=energy_ylims,
        legend=:topright
    )
    vline!(p_mean_d, [t]; color=:gray35, linestyle=:dash, alpha=0.7, label="current")

    p_mean_g = plot(
        t_hist,
        E_gauss[1:idx];
        color=:chocolate,
        linewidth=2.0,
        xlabel="time",
        ylabel="mean local E",
        title="Gaussian proposal: mean local energy",
        label="mean E",
        xlims=(0, nsteps * dt),
        ylims=energy_ylims,
        legend=:topright
    )
    vline!(p_mean_g, [t]; color=:gray35, linestyle=:dash, alpha=0.7, label="current")

    p_var_d = plot(
        t_hist,
        V_drift[1:idx];
        color=:purple4,
        linewidth=2.0,
        xlabel="time",
        ylabel="Var(E_L)",
        title="Importance Sampling proposal: local-energy variance",
        label="variance",
        xlims=(0, nsteps * dt),
        ylims=var_ylims,
        legend=:topright
    )
    vline!(p_var_d, [t]; color=:gray35, linestyle=:dash, alpha=0.7, label="current")

    p_var_g = plot(
        t_hist,
        V_gauss[1:idx];
        color=:firebrick3,
        linewidth=2.0,
        xlabel="time",
        ylabel="Var(E_L)",
        title="Gaussian proposal: local-energy variance",
        label="variance",
        xlims=(0, nsteps * dt),
        ylims=var_ylims,
        legend=:topright
    )
    vline!(p_var_g, [t]; color=:gray35, linestyle=:dash, alpha=0.7, label="current")

    plot(
        p_dist_d, p_dist_g,
        p_mean_d, p_mean_g,
        p_var_d, p_var_g;
        layout=(3, 2),
        size=(1800, 1200)
    )
end

Logging.with_logger(Logging.NullLogger()) do
    redirect_stderr(devnull) do
        gif(anim, output_file; fps=fps)
    end
end

println("Saved regularized-H1D VMC proposal-comparison animation: ", output_file)
println("Frames: ", nframes, " | Snapshot steps: ", frame_steps)
println("Final acceptance rates -> importance_sampling proposal: ", round(sim_drift.acceptance_rate; digits=4),
        ", gaussian proposal: ", round(sim_gauss.acceptance_rate; digits=4))

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
