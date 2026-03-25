using Random
using Statistics
using Printf
using LinearAlgebra
using Plots
using Logging

include("../notebooks/formatting.jl")

if !isdefined(Main, :System1D)
    include("../../../../../src/System1D.jl")
end
using .System1D: Hamiltonian, TrialWF, DMCParams, ImportanceGuiding, run_dmc, plot_snapshot_1d_density
using .System1D: PeriodicBoundary1D
using .System1D: local_energy

# System: periodic cosine lattice
L = 10.0
V0 = -1.0
k = 2pi / L
lambda_ = -V0 / 2

bc = PeriodicBoundary1D(0.0, L)
cosine_potential(R) = V0 * cos(k * R[1])
H = Hamiltonian(1, 0.5, cosine_potential, bc)

logpsi(R) = lambda_ * cos(k * R[1])
gradlogpsi(R) = [-lambda_ * k * sin(k * R[1])]
lapllogpsi(R) = -lambda_ * k^2 * cos(k * R[1])
trial = TrialWF(logpsi, gradlogpsi, lapllogpsi)
guiding = ImportanceGuiding(trial, H)

local_energy_theory(x) = (V0 + 0.5 * lambda_ * k^2) * cos(k * x) - 0.5 * lambda_^2 * k^2 * sin(k * x)^2

function padded_quantile_limits(values::AbstractVector{<:Real};
    qlo::Float64=0.01,
    qhi::Float64=0.99,
    pad_frac::Float64=0.08)
    0.0 <= qlo < qhi <= 1.0 || throw(ArgumentError("Quantile limits must satisfy 0 <= qlo < qhi <= 1"))
    finite_vals = Float64[v for v in values if isfinite(v)]
    isempty(finite_vals) && return (-1.0, 1.0)

    lo, hi = quantile(finite_vals, (qlo, qhi))
    span = hi - lo
    if span <= eps(Float64)
        center = 0.5 * (lo + hi)
        pad = max(1e-6, 0.05 * (abs(center) + 1.0))
        return center - pad, center + pad
    end
    pad = pad_frac * span
    return lo - pad, hi + pad
end

function mask_outside_range(values::AbstractVector{<:Real}, lo::Real, hi::Real)
    out = Vector{Float64}(undef, length(values))
    nkeep = 0
    for i in eachindex(values)
        v = Float64(values[i])
        if isfinite(v) && lo <= v <= hi
            out[i] = v
            nkeep += 1
        else
            out[i] = NaN
        end
    end
    return out, nkeep
end

function smoothed_density_curve(snap::AbstractVector{<:AbstractVector};
    nbins::Int=120,
    smoothing_window::Int=11,
    xmin::Real=0.0,
    xmax::Real=L)
    nbins >= 2 || throw(ArgumentError("nbins must be >= 2"))
    smoothing_window >= 1 || throw(ArgumentError("smoothing_window must be >= 1"))
    xmax > xmin || throw(ArgumentError("xmax must be greater than xmin"))

    xs = Float64[R[1] for R in snap]
    isempty(xs) && error("Cannot compute smoothed density for an empty snapshot.")

    width = (xmax - xmin) / nbins
    counts = zeros(Float64, nbins)
    for x in xs
        idx = floor(Int, (x - xmin) / width) + 1
        idx = clamp(idx, 1, nbins)
        counts[idx] += 1.0
    end

    dens = counts ./ (length(xs) * width)
    half = fld(smoothing_window, 2)
    smooth = similar(dens)
    for i in eachindex(dens)
        ilo = max(firstindex(dens), i - half)
        ihi = min(lastindex(dens), i + half)
        smooth[i] = sum(@view dens[ilo:ihi]) / (ihi - ilo + 1)
    end

    centers = Float64[xmin + (i - 0.5) * width for i in 1:nbins]
    return centers, smooth
end

function centered_moving_average_ignore_nan(values::AbstractVector{<:Real}; window::Int=31)
    window >= 1 || throw(ArgumentError("window must be >= 1"))
    n = length(values)
    out = Vector{Float64}(undef, n)
    half = fld(window, 2)
    for i in eachindex(values)
        ilo = max(firstindex(values), i - half)
        ihi = min(lastindex(values), i + half)
        total = 0.0
        count = 0
        for j in ilo:ihi
            v = Float64(values[j])
            if isfinite(v)
                total += v
                count += 1
            end
        end
        out[i] = count > 0 ? (total / count) : NaN
    end
    return out
end

# Check periodicity and local-energy algebra
x_test = 0.37 * L
@assert isapprox(cosine_potential([x_test]), cosine_potential([x_test + L]); atol=1e-12)
@assert isapprox(logpsi([x_test]), logpsi([x_test + L]); atol=1e-12)
@assert isapprox(local_energy_theory(x_test), local_energy(guiding, [x_test]); atol=1e-12)

# Finite-difference reference for E0 and |psi0|^2 under PBC
Nref = 400
dx = L / Nref
xref = collect(0:(Nref - 1)) .* dx

Hmat = zeros(Float64, Nref, Nref)
kin_diag = 1.0 / dx^2
kin_off = -0.5 / dx^2
for i in 1:Nref
    Hmat[i, i] = kin_diag + V0 * cos(k * xref[i])
end
for i in 1:(Nref - 1)
    Hmat[i, i + 1] = kin_off
    Hmat[i + 1, i] = kin_off
end
Hmat[1, Nref] = kin_off
Hmat[Nref, 1] = kin_off

evals, evecs = eigen(Symmetric(Hmat))
E0_ref = evals[1]
psi0 = evecs[:, 1]
rho0_ref = abs2.(psi0)
rho0_ref ./= (sum(rho0_ref) * dx)

rho_trial = exp.(2 .* lambda_ .* cos.(k .* xref))
rho_trial ./= (sum(rho_trial) * dx)

targetN = 2000
dt = 3.0e-3
nsteps = 3000
nequil = 400
# Three snapshot times for smoothed-density comparisons: initial, mid, final.
snapshot_steps = [0, div(nsteps, 2), nsteps]
energy_plot_quantile_window = (0.01, 0.99)
et_smoothing_window = 31

params_plain = DMCParams(; dt=dt, nsteps=nsteps, nequil=nequil, targetN=targetN, ET0=-0.2, branch_cap=8, nblocks=50)
params_guided = DMCParams(; dt=dt, nsteps=nsteps, nequil=nequil, targetN=targetN, ET0=-0.2, branch_cap=8, nblocks=50)

rng_init = MersenneTwister(1234)
base_positions = [[rand(rng_init) * L] for _ in 1:targetN]

rng_plain = MersenneTwister(52)
rng_guided = MersenneTwister(77)

sim_plain = run_dmc(H, params_plain, base_positions; rng=rng_plain, snapshot_steps=snapshot_steps)
sim_guided = run_dmc(H, params_guided, base_positions; rng=rng_guided, guiding=guiding, snapshot_steps=snapshot_steps)

t = (0:nsteps) .* dt

Ehist_plain = sim_plain.energy_mean_history
Ehist_guided = sim_guided.energy_mean_history
Varhist_plain = sim_plain.energy_variance_history
Varhist_guided = sim_guided.energy_variance_history
EThist_plain = sim_plain.ET_history
EThist_guided = sim_guided.ET_history

start_idx = min(nequil + 1, length(Ehist_plain))
post_ET_plain = EThist_plain[start_idx:end]
post_ET_guided = EThist_guided[start_idx:end]
Ebar_plain = mean(post_ET_plain)
Ebar_guided = mean(post_ET_guided)
SEM_plain = (length(post_ET_plain) > 1) ? std(post_ET_plain) / sqrt(length(post_ET_plain)) : NaN
SEM_guided = (length(post_ET_guided) > 1) ? std(post_ET_guided) / sqrt(length(post_ET_guided)) : NaN
varbar_plain = mean(Varhist_plain[start_idx:end])
varbar_guided = mean(Varhist_guided[start_idx:end])

Eplain_str, dEplain_str = format_with_uncertainty(Ebar_plain, SEM_plain)
Eguid_str, dEguid_str = format_with_uncertainty(Ebar_guided, SEM_guided)
Eref_str = format_sigfig(E0_ref; sigfigs=7)

display_lo, display_hi = padded_quantile_limits(
    vcat(EThist_plain[start_idx:end], EThist_guided[start_idx:end]);
    qlo=energy_plot_quantile_window[1],
    qhi=energy_plot_quantile_window[2],
    pad_frac=0.06
)
EThist_plain_plot, kept_plain = mask_outside_range(EThist_plain, display_lo, display_hi)
EThist_guided_plot, kept_guided = mask_outside_range(EThist_guided, display_lo, display_hi)
EThist_plain_smooth = centered_moving_average_ignore_nan(EThist_plain_plot; window=et_smoothing_window)
EThist_guided_smooth = centered_moving_average_ignore_nan(EThist_guided_plot; window=et_smoothing_window)

println("Cosine PBC DMC (unguided, ET post-eq): E = $(Eplain_str) +/- $(dEplain_str)")
println("Cosine PBC DMC (guided, ET post-eq):   E = $(Eguid_str) +/- $(dEguid_str)")
println("Finite-difference periodic reference E0 ~ $(Eref_str)")
println(@sprintf("Variance reduction factor (unguided/guided) ~ %.3f", varbar_plain / varbar_guided))
println(@sprintf(
    "Energy-plot display filter: [%.6f, %.6f] | kept points (unguided/guided): %d/%d, %d/%d",
    display_lo, display_hi, kept_plain, length(EThist_plain), kept_guided, length(EThist_guided)
))
println("ET plot smoothing window (points): ", et_smoothing_window)

figdir = abspath(joinpath(@__DIR__, "..", "outputs", "figures"))
mkpath(figdir)

xplot = range(0.0, L; length=600)
Vplot = [V0 * cos(k * x) for x in xplot]
ELplot = [local_energy_theory(x) for x in xplot]

Logging.with_logger(Logging.NullLogger()) do
    redirect_stderr(devnull) do
        p1 = plot(
            xplot,
            Vplot;
            xlabel="x",
            ylabel="energy",
            title="Cosine lattice on ring: potential and trial local energy",
            linewidth=2.2,
            color=:navy,
            label="V(x)=V0 cos(kx)"
        )
        plot!(p1, xplot, ELplot; linewidth=2.2, color=:darkorange, linestyle=:dash, label="E_L(x) from trial")
        hline!(p1, [E0_ref]; color=:black, linestyle=:dot, linewidth=1.8, label="FD reference E0")
        savefig(p1, joinpath(figdir, "cosine_potential_and_trial_local_energy.png"))

        p2 = plot(
            t,
            EThist_plain_smooth;
            xlabel="imaginary time tau",
            ylabel="E_T(tau)",
            title="Cosine lattice: reference energy comparison (smoothed)",
            linewidth=2.1,
            color=:firebrick,
            label="unguided (smoothed)",
            ylims=(display_lo, display_hi)
        )
        plot!(p2, t, EThist_guided_smooth; linewidth=2.1, color=:teal, label="guided (smoothed)")
        hline!(p2, [E0_ref]; color=:black, linestyle=:dash, linewidth=1.8, label="FD reference E0")
        savefig(p2, joinpath(figdir, "cosine_ET_vs_tau_guided_vs_unguided.png"))

        p3 = plot(
            t,
            Ehist_plain;
            xlabel="imaginary time tau",
            ylabel="mean estimator",
            title="Cosine lattice: mean estimator history",
            linewidth=2.1,
            color=:firebrick,
            label="unguided mean(V)"
        )
        plot!(p3, t, Ehist_guided; linewidth=2.1, color=:teal, label="guided mean(E_L)")
        hline!(p3, [E0_ref]; color=:black, linestyle=:dash, linewidth=1.8, label="FD reference E0")
        savefig(p3, joinpath(figdir, "cosine_mean_energy_vs_tau_guided_vs_unguided.png"))

        p4 = plot(
            t,
            Varhist_plain;
            xlabel="imaginary time tau",
            ylabel="Var(E estimator)",
            title="Cosine lattice: variance reduction with guiding",
            linewidth=2.1,
            color=:firebrick,
            label="unguided"
        )
        plot!(p4, t, Varhist_guided; linewidth=2.1, color=:teal, label="guided")
        savefig(p4, joinpath(figdir, "cosine_variance_vs_tau_guided_vs_unguided.png"))

        snap_plain = sim_plain.walker_positions_history[end]
        snap_guided = sim_guided.walker_positions_history[end]
        p5a = plot_snapshot_1d_density(
            snap_plain;
            nbins=120,
            xmin=0.0,
            xmax=L,
            normalize=true,
            title="Final density (unguided)",
            curve_label="DMC unguided"
        )
        plot!(p5a, xref, rho0_ref; color=:black, linewidth=2.2, linestyle=:dash, label="FD |psi0|^2")
        plot!(p5a, xref, rho_trial; color=:gray30, linewidth=1.8, linestyle=:dot, label="|psi_T|^2 (guide)")

        p5b = plot_snapshot_1d_density(
            snap_guided;
            nbins=120,
            xmin=0.0,
            xmax=L,
            normalize=true,
            title="Final density (guided)",
            curve_label="DMC guided"
        )
        plot!(p5b, xref, rho0_ref; color=:black, linewidth=2.2, linestyle=:dash, label="FD |psi0|^2")
        plot!(p5b, xref, rho_trial; color=:gray30, linewidth=1.8, linestyle=:dot, label="|psi_T|^2 (guide)")

        p5 = plot(p5a, p5b; layout=(2, 1), size=(900, 900))
        savefig(p5, joinpath(figdir, "cosine_final_density_guided_vs_unguided.png"))

        p6a = plot(
            xlabel="x",
            ylabel="density",
            xlims=(0.0, L),
            title="Cosine lattice smoothed snapshots (unguided)",
            legend=:topright
        )
        p6b = plot(
            xlabel="x",
            ylabel="density",
            xlims=(0.0, L),
            title="Cosine lattice smoothed snapshots (guided)",
            legend=:topright
        )
        snapshot_colors = [:royalblue, :darkorange, :forestgreen]
        for (i, s) in enumerate(snapshot_steps)
            tau_label = "tau=$(format_sigfig(s * dt; sigfigs=3))"
            curve_color = snapshot_colors[mod1(i, length(snapshot_colors))]
            xc_plain, yc_plain = smoothed_density_curve(
                sim_plain.walker_positions_history[i];
                nbins=120,
                smoothing_window=11,
                xmin=0.0,
                xmax=L
            )
            xc_guided, yc_guided = smoothed_density_curve(
                sim_guided.walker_positions_history[i];
                nbins=120,
                smoothing_window=11,
                xmin=0.0,
                xmax=L
            )
            plot!(p6a, xc_plain, yc_plain; linewidth=2.2, color=curve_color, linestyle=:solid, label=tau_label)
            plot!(p6b, xc_guided, yc_guided; linewidth=2.2, color=curve_color, linestyle=:solid, label=tau_label)
        end
        plot!(p6a, xref, rho0_ref; color=:black, linewidth=2.5, linestyle=:dash, label="FD |psi0|^2")
        plot!(p6b, xref, rho0_ref; color=:black, linewidth=2.5, linestyle=:dash, label="FD |psi0|^2")
        savefig(p6a, joinpath(figdir, "cosine_snapshot_smoothed_unguided_vs_reference.png"))
        savefig(p6b, joinpath(figdir, "cosine_snapshot_smoothed_guided_vs_reference.png"))
    end
end

println("Saved figures to: ", figdir)
