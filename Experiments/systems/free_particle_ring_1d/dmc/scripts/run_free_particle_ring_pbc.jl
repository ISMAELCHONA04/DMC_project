using Random
using Statistics
using Printf
using Plots
using Logging

include("../notebooks/formatting.jl")

if !isdefined(Main, :System1D)
    include("../../../../../src/System1D.jl")
end
using .System1D: Hamiltonian, DMCParams, run_dmc, plot_snapshot_1d_density
using .System1D: PeriodicBoundary1D

# System: free particle on a periodic ring
L = 10.0
bc = PeriodicBoundary1D(0.0, L)
V_free(R) = 0.0
H = Hamiltonian(1, 0.5, V_free, bc)

targetN = 2000
dt = 1.0e-3
nsteps = 1200
nequil = 200
params = DMCParams(; dt=dt, nsteps=nsteps, nequil=nequil, targetN=targetN, ET0=0.0, branch_cap=8, nblocks=50)

rng_init = MersenneTwister(1234)
initial_positions = [[rand(rng_init) * L] for _ in 1:targetN]

snapshot_steps = [0, 40, 120, 400, 800, nsteps]
rng_sim = MersenneTwister(42)
sim = run_dmc(H, params, initial_positions; rng=rng_sim, snapshot_steps=snapshot_steps)

ET_history = sim.ET_history
population_history = sim.population_history
energy_mean_history = sim.energy_mean_history
energy_variance_history = sim.energy_variance_history

t = (0:params.nsteps) .* params.dt
start_idx = min(params.nequil + 1, length(energy_mean_history))
postE = energy_mean_history[start_idx:end]
Ebar = mean(postE)
SEM = (length(postE) > 1) ? std(postE) / sqrt(length(postE)) : NaN
Ebar_str, SEM_str = format_with_uncertainty(Ebar, SEM)
println("Free-ring DMC energy (post-eq, nequil=$(params.nequil)): E = $(Ebar_str) +/- $(SEM_str)")

figdir = abspath(joinpath(@__DIR__, "..", "outputs", "figures"))
mkpath(figdir)

xgrid = range(0.0, L; length=400)
uniform_density = fill(1.0 / L, length(xgrid))

Logging.with_logger(Logging.NullLogger()) do
    redirect_stderr(devnull) do
        p1 = plot(
            t,
            ET_history;
            xlabel="imaginary time tau",
            ylabel="E_T(tau)",
            title="Free particle ring: reference energy",
            linewidth=2.1,
            color=:navy,
            label="DMC E_T"
        )
        hline!(p1, [0.0]; linestyle=:dash, color=:black, linewidth=1.8, label="Theory E0 = 0")
        savefig(p1, joinpath(figdir, "free_ring_ET_vs_tau.png"))

        p2 = plot(
            t,
            population_history;
            xlabel="imaginary time tau",
            ylabel="N_w(tau)",
            title="Free particle ring: population control",
            linewidth=2.1,
            color=:darkgreen,
            label="walker population"
        )
        hline!(p2, [targetN]; linestyle=:dash, color=:black, linewidth=1.8, label="target N")
        savefig(p2, joinpath(figdir, "free_ring_population_vs_tau.png"))

        p3 = plot(
            t,
            energy_mean_history;
            xlabel="imaginary time tau",
            ylabel="<V>(tau)",
            title="Free particle ring: mean sampled energy",
            linewidth=2.1,
            color=:darkorange,
            label="DMC mean energy"
        )
        hline!(p3, [0.0]; linestyle=:dash, color=:black, linewidth=1.8, label="Theory E0 = 0")
        savefig(p3, joinpath(figdir, "free_ring_mean_energy_vs_tau.png"))

        snap_final = sim.walker_positions_history[end]
        p4 = plot_snapshot_1d_density(
            snap_final;
            nbins=120,
            xmin=0.0,
            xmax=L,
            normalize=true,
            title="Free particle ring: final density vs uniform theory",
            curve_label="DMC density"
        )
        plot!(
            p4,
            xgrid,
            uniform_density;
            linestyle=:dash,
            linewidth=2.3,
            color=:black,
            label="Theory rho(x)=1/L"
        )
        savefig(p4, joinpath(figdir, "free_ring_density_final_vs_uniform.png"))

        p5 = plot(
            xlabel="x",
            ylabel="density",
            xlims=(0.0, L),
            title="Free particle ring: density snapshots",
            legend=:topright
        )
        for (i, s) in enumerate(snapshot_steps)
            snap = sim.walker_positions_history[i]
            xs = [r[1] for r in snap]
            histogram!(
                p5,
                xs;
                bins=80,
                normalize=:pdf,
                alpha=0.2,
                linealpha=0.8,
                linewidth=1.2,
                label="tau=$(format_sigfig(s * dt; sigfigs=3))"
            )
        end
        plot!(
            p5,
            xgrid,
            uniform_density;
            linestyle=:dash,
            linewidth=2.5,
            color=:black,
            label="Theory rho(x)=1/L"
        )
        savefig(p5, joinpath(figdir, "free_ring_density_snapshots_vs_uniform.png"))
    end
end

println("Saved figures to: ", figdir)
