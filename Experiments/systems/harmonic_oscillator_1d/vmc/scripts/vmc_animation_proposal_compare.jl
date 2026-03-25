using Random
using Statistics
using Plots
using Logging
using Printf

include("../../../../../src/System1D.jl")
using .System1D

function main()

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
omega = 1.0
alpha =  parse(Float64, get(ENV, "VMC_ALPHA", "2.0"))

# Defaults are tuned for the repository's harmonic-oscillator VMC comparison.
targetN = parse(Int, get(ENV, "VMC_TARGETN", "50000"))
dt = parse(Float64, get(ENV, "VMC_DT", "0.025"))
nsteps = parse(Int, get(ENV, "VMC_NSTEPS", "120"))

nbins = parse(Int, get(ENV, "VMC_NBINS", "180"))
smoothing_window = parse(Int, get(ENV, "VMC_SMOOTH_WINDOW", "9"))
fps = parse(Int, get(ENV, "VMC_GIF_FPS", "5"))
output_dir = joinpath(dirname(@__DIR__), "outputs", "animations")
mkpath(output_dir)
output_file = joinpath(output_dir, get(ENV, "VMC_GIF_NAME", "vmc_proposal_compare_alpha_2.gif"))
write_gif = lowercase(get(ENV, "VMC_WRITE_GIF", "true")) in ("1", "true", "yes", "y")

figures_dir = joinpath(dirname(@__DIR__), "outputs", "figures")
mkpath(figures_dir)
convergence_energy_file = joinpath(figures_dir, get(ENV, "VMC_ENERGY_FIG_NAME", "vmc_alpha_2_energy_convergence.png"))
convergence_variance_file = joinpath(figures_dir, get(ENV, "VMC_VARIANCE_FIG_NAME", "vmc_alpha_2_variance_convergence.png"))

# Optional dual-panel frame-sequence export (for Beamer/animategraphics).
export_dual_frames = lowercase(get(ENV, "VMC_EXPORT_DUAL_FRAMES", "true")) in ("1", "true", "yes", "y")
frames_root = joinpath(dirname(@__DIR__), "outputs", "frames")
mkpath(frames_root)
dual_frames_dir = get(ENV, "VMC_DUAL_FRAMES_DIR", joinpath(frames_root, "dual_dist_alpha_2"))
dual_prefix = get(ENV, "VMC_DUAL_FRAME_PREFIX", "dual_alpha_2")
frame_pad = parse(Int, get(ENV, "VMC_FRAME_PAD", "4"))

# Snapshot schedule for distribution overlays.
snapshot_candidates = collect(0:nsteps)
snapshot_steps = sort(unique(filter(s -> 0 <= s <= nsteps, snapshot_candidates)))

# -----------------------------------------------------------------------------
# System + trial wavefunction
# -----------------------------------------------------------------------------
V(R) = 0.5 * omega^2 * R[1]^2
H = Hamiltonian(1, 0.5, V)

logpsi(R) = -0.5 * alpha * R[1]^2
gradlogpsi(R) = [-alpha * R[1]]
lapllogpsi(R) = -alpha
trial = TrialWF(logpsi, gradlogpsi, lapllogpsi)

params = VMCParams(; dt=dt, nsteps=nsteps, targetN=targetN, ET0=0.5)

# Shared initial positions for fair comparison.
rng_init = MersenneTwister(1234)
base_positions = [[rand(rng_init)+1] for _ in 1:targetN]

# Importance-sampling (drifted Gaussian) proposal simulation
sim_drift = run_vmc(
    H,
    params,
    base_positions,
    trial;
    rng=MersenneTwister(42),
    proposal=DriftGaussianProposal(),
    snapshot_steps=snapshot_steps,
)

# Pure Gaussian proposal simulation
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
    nbins::Int=180,
    smoothing_window::Int=9,
    xmin::Float64=-4.0,
    xmax::Float64=4.0
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
        return (-4.0, 4.0)
    end

    if lo == hi
        return (lo - 1.0, hi + 1.0)
    end

    pad = 0.08 * (hi - lo)
    return (lo - pad, hi + pad)
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

xlo, xhi = _combined_x_limits(snaps_drift, snaps_gauss)

# Target sampled density from trial wavefunction: |psi_T|^2
target_density(x) = sqrt(alpha / pi) * exp(-alpha * x^2)
xref = range(xlo, xhi; length=800)
yref = Float64[target_density(x) for x in xref]

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
# Theoretical references for the 1D HO (a.u., omega = 1 -> E0 = 0.5).
E_exact = 0.625
var_exact = 0.28125

# Typography / styling.
legend_fs = parse(Int, get(ENV, "VMC_LEGEND_FS", "18"))
guide_fs = parse(Int, get(ENV, "VMC_GUIDE_FS", "24"))
tick_fs = parse(Int, get(ENV, "VMC_TICK_FS", "16"))
title_fs = parse(Int, get(ENV, "VMC_TITLE_FS", "26"))
caption_fs = parse(Int, get(ENV, "VMC_CAPTION_FS", "15"))

# Reserve negative-y space for dynamic per-panel captions.
caption_band = 0.48 * dist_max
dist_ylims_with_caption = (-caption_band, dist_ylims[2])

# -----------------------------------------------------------------------------
# Animate: 1 row x 2 cols (distribution panels only)
# -----------------------------------------------------------------------------
function _build_distribution_dual_frame(i::Int)
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
        ylabel="walker density",
        label="proposal: N(R_i(t)+F_i*dt, dt)",
        xlims=(xlo, xhi),
        ylims=dist_ylims_with_caption,
        legend=:topright,
        legendfontsize=legend_fs,
        guidefontsize=guide_fs,
        tickfontsize=tick_fs,
        titlefontsize=title_fs
    )
    plot!(
        p_dist_d,
        xref,
        yref;
        color=:black,
        linestyle=:dash,
        linewidth=2.0,
        label="target |psi_T|^2"
    )

    p_dist_g = plot(
        centers_g,
        smoothed_g;
        color=:darkorange,
        linewidth=2.2,
        xlabel="x",
        ylabel="walker density",
        label="proposal: N(R_i(t), dt)",
        xlims=(xlo, xhi),
        ylims=dist_ylims_with_caption,
        legend=:topright,
        legendfontsize=legend_fs,
        guidefontsize=guide_fs,
        tickfontsize=tick_fs,
        titlefontsize=title_fs
    )
    plot!(
        p_dist_g,
        xref,
        yref;
        color=:black,
        linestyle=:dash,
        linewidth=2.0,
        label="target |psi_T|^2"
    )

    idx = min(step + 1, length(E_drift))
    e_d = E_drift[idx]
    e_g = E_gauss[idx]
    v_d = V_drift[idx]
    v_g = V_gauss[idx]

    cap_x = xlo + 0.02 * (xhi - xlo)
    cap_y1 = -0.15 * dist_max
    cap_y2 = -0.30 * dist_max

    cap1_d = @sprintf("Energy: <E_L>(t)=%.6f a.u. | Theory: E=%.6f a.u.", e_d, E_exact)
    cap2_d = @sprintf("Variance: Var(E_L)(t)=%.6f (a.u.)^2 | Theory: Var=%.6f (a.u.)^2", v_d, var_exact)
    cap1_g = @sprintf("Energy: <E_L>(t)=%.6f a.u. | Theory: E=%.6f a.u.", e_g, E_exact)
    cap2_g = @sprintf("Variance: Var(E_L)(t)=%.6f (a.u.)^2 | Theory: Var=%.6f (a.u.)^2", v_g, var_exact)

    annotate!(p_dist_d, cap_x, cap_y1, text(cap1_d, caption_fs, :black, :left))
    annotate!(p_dist_d, cap_x, cap_y2, text(cap2_d, caption_fs, :black, :left))

    annotate!(p_dist_g, cap_x, cap_y1, text(cap1_g, caption_fs, :black, :left))
    annotate!(p_dist_g, cap_x, cap_y2, text(cap2_g, caption_fs, :black, :left))

    return plot(p_dist_d, p_dist_g; layout=(1, 2), size=(2200, 950))
end

anim = @animate for i in 1:nframes
    _build_distribution_dual_frame(i)
end

if write_gif
    Logging.with_logger(Logging.NullLogger()) do
        redirect_stderr(devnull) do
            gif(anim, output_file; fps=fps)
        end
    end
    println("Saved VMC distribution-comparison animation: ", output_file)
else
    println("Skipped GIF export (VMC_WRITE_GIF=false).")
end

println("Frames: ", nframes, " | Snapshot steps: ", frame_steps)
println("Final acceptance rates -> importance_sampling proposal: ", round(sim_drift.acceptance_rate; digits=4),
        ", gaussian proposal: ", round(sim_gauss.acceptance_rate; digits=4))

# -----------------------------------------------------------------------------
# Convergence figures (no titles; to be captioned in LaTeX slides)
# -----------------------------------------------------------------------------
t_hist = (0:nsteps) .* dt

pE = plot(
    t_hist,
    E_drift;
    xlabel="time",
    ylabel="mean local energy (a.u.)",
    label="importance-sampled",
    color=:navy,
    linewidth=2.2,
    legend=:topright,
    legendfontsize=11,
    guidefontsize=12,
    tickfontsize=10
)
plot!(
    pE,
    t_hist,
    E_gauss;
    label="gaussian",
    color=:darkorange,
    linewidth=2.2
)
hline!(
    pE,
    [E_exact];
    label=@sprintf("theory = %.6f", E_exact),
    color=:black,
    linestyle=:dash,
    linewidth=2.0
)
savefig(pE, convergence_energy_file)

pV = plot(
    t_hist,
    V_drift;
    xlabel="time",
    ylabel="variance (a.u.^2)",
    label="importance-sampled",
    color=:purple4,
    linewidth=2.2,
    legend=:topright,
    legendfontsize=11,
    guidefontsize=12,
    tickfontsize=10
)
plot!(
    pV,
    t_hist,
    V_gauss;
    label="gaussian",
    color=:firebrick3,
    linewidth=2.2
)
hline!(
    pV,
    [var_exact];
    label=@sprintf("theory = %.6f", var_exact),
    color=:black,
    linestyle=:dash,
    linewidth=2.0
)
savefig(pV, convergence_variance_file)

println("Saved convergence figure: ", abspath(convergence_energy_file))
println("Saved convergence figure: ", abspath(convergence_variance_file))

function _clean_existing_frames!(dir::AbstractString, prefix::AbstractString)
    if !isdir(dir)
        return
    end
    regex = Regex("^" * escape_string(prefix) * raw"\d+\.png$")
    for name in readdir(dir)
        if occursin(regex, name)
            rm(joinpath(dir, name); force=true)
        end
    end
end

if export_dual_frames
    mkpath(dual_frames_dir)
    _clean_existing_frames!(dual_frames_dir, dual_prefix)

    for i in 1:nframes
        idx_str = lpad(string(i - 1), frame_pad, '0')
        dual_panel = _build_distribution_dual_frame(i)
        savefig(dual_panel, joinpath(dual_frames_dir, dual_prefix * idx_str * ".png"))
    end

    println("Saved dual-panel distribution frames with captions to: ", abspath(dual_frames_dir))
    println("Dual frames: ", nframes, " | Prefix: ", dual_prefix, " | Padding: ", frame_pad)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
