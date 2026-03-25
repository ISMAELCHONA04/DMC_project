using Random
using Statistics
using Printf
using LinearAlgebra
using Plots

const GFMC_BENCHMARK_DIR = @__DIR__
const EXPERIMENTS_DIR = normpath(joinpath(GFMC_BENCHMARK_DIR, "..", ".."))
const PROJECT_ROOT = normpath(joinpath(EXPERIMENTS_DIR, ".."))

include(joinpath(EXPERIMENTS_DIR, "common", "notebook_helpers.jl"))
if !isdefined(Main, :System1D)
    include(joinpath(PROJECT_ROOT, "src", "System1D.jl"))
end
using .System1D
include(joinpath(EXPERIMENTS_DIR, "systems", "periodic_ion_ring_1d", "gfmc", "notebooks", "periodic_ion_ring_helpers.jl"))

default(; dpi=170)

const BENCHMARK_CASE_IDS = [
    "free_particle_ring",
    "harmonic_oscillator_unguided",
    "harmonic_oscillator_guided",
    "cosine_lattice_ring",
    "hydrogen_unguided",
    "hydrogen_fixed_node",
    "two_particle_ho_guided",
    "two_particle_ho_fixed_node",
    "fermion_ring_fixed_node",
    "hardcore_boson_ring",
    "periodic_ion_single_particle",
    "periodic_ion_spinless_fermions",
    "periodic_ion_bosons_jastrow",
    "periodic_ion_bosons_tg",
]

const ACCURACY_CASE_IDS = Set([
    "free_particle_ring",
    "harmonic_oscillator_unguided",
    "harmonic_oscillator_guided",
    "cosine_lattice_ring",
    "two_particle_ho_guided",
    "two_particle_ho_fixed_node",
    "periodic_ion_single_particle",
    "periodic_ion_spinless_fermions",
])

const STRESS_CASE_IDS = Set(filter(id -> !(id in ACCURACY_CASE_IDS), BENCHMARK_CASE_IDS))

const BENCHMARK_COLORS = (
    navy=:navy,
    teal=:teal,
    orange=:darkorange,
    red=:firebrick,
    green=:forestgreen,
    purple=:purple,
    gold=:goldenrod,
)

function benchmark_case_selection(selection::AbstractString)
    token = lowercase(strip(selection))
    if token == "all"
        return copy(BENCHMARK_CASE_IDS)
    elseif token == "accuracy"
        return [id for id in BENCHMARK_CASE_IDS if id in ACCURACY_CASE_IDS]
    elseif token == "stress"
        return [id for id in BENCHMARK_CASE_IDS if id in STRESS_CASE_IDS]
    end

    ids = filter(!isempty, strip.(split(selection, ',')))
    isempty(ids) && error("No benchmark cases selected from: $selection")
    unknown = filter(id -> !(id in BENCHMARK_CASE_IDS), ids)
    isempty(unknown) || error("Unknown benchmark case id(s): $(join(unknown, ", "))")
    return ids
end

function benchmark_paths(case_id::AbstractString, tier::Symbol)
    root = joinpath(GFMC_BENCHMARK_DIR, "outputs", case_id, String(tier))
    figures_dir = joinpath(root, "figures")
    tables_dir = joinpath(root, "tables")
    logs_dir = joinpath(root, "logs")
    mkpath(figures_dir)
    mkpath(tables_dir)
    mkpath(logs_dir)
    return (root=root, figures_dir=figures_dir, tables_dir=tables_dir, logs_dir=logs_dir)
end

function suite_paths(suite_name::AbstractString, tier::Symbol)
    root = joinpath(GFMC_BENCHMARK_DIR, "outputs", "suite", suite_name, String(tier))
    figures_dir = joinpath(root, "figures")
    tables_dir = joinpath(root, "tables")
    logs_dir = joinpath(root, "logs")
    mkpath(figures_dir)
    mkpath(tables_dir)
    mkpath(logs_dir)
    return (root=root, figures_dir=figures_dir, tables_dir=tables_dir, logs_dir=logs_dir)
end

function _rounded_int(x::Real)
    return max(1, Int(round(Float64(x))))
end

function make_gfmc_params(;
    dt::Real,
    tau_total::Real,
    targetN::Integer,
    ET0::Real,
    feedback::Real,
    reconfiguration_interval::Integer,
    branch_cap::Real=5.0,
    energy_window::Integer=20,
    equil_fraction::Real=0.2,
)
    dt_f = Float64(dt)
    tau_f = Float64(tau_total)
    nsteps = max(20, Int(ceil(tau_f / dt_f)))
    nequil = clamp(_rounded_int(equil_fraction * nsteps), 1, max(1, nsteps - 1))
    return GFMCParams(
        dt_f,
        nsteps,
        nequil,
        Int(targetN),
        Float64(ET0),
        Float64(feedback),
        Int(reconfiguration_interval),
        Float64(branch_cap),
        Int(energy_window),
    )
end

function make_vmc_params(;
    dt::Real,
    tau_total::Real,
    targetN::Integer,
    ET0::Real,
)
    nsteps = max(20, Int(ceil(Float64(tau_total) / Float64(dt))))
    return VMCParams(; dt=Float64(dt), nsteps=nsteps, targetN=Int(targetN), ET0=Float64(ET0))
end

function pooled_coordinates(snapshot)
    xs = Float64[]
    for R in snapshot
        append!(xs, Float64.(R))
    end
    return xs
end

function pair_separations(snapshot, distance_fn)
    rs = Float64[]
    for R in snapshot
        for i in 1:(length(R) - 1)
            for j in (i + 1):length(R)
                push!(rs, distance_fn(R[i], R[j]))
            end
        end
    end
    return rs
end

function benchmark_summary(sim)
    start_idx = min(sim.params.nequil + 1, length(sim.energy_mean_history))
    mean_energy, sem_energy = nb_mean_sem(sim.energy_mean_history[start_idx:end])
    return (
        mean_energy=mean_energy,
        sem_energy=sem_energy,
        final_ET=sim.ET_history[end],
        final_population=sim.population_history[end],
        final_mean_weight=sim.mean_weight_history[end],
        final_effective_population=sim.effective_population_history[end],
        final_acceptance=sim.acceptance_history[end],
        final_variance=sim.energy_variance_history[end],
    )
end

function density_from_snapshot(
    snapshot;
    periodic::Bool,
    xmin::Real,
    xmax::Real,
    bandwidth::Real=0.1,
    grid_points::Integer=400,
    nbins::Integer=140,
)
    xs = pooled_coordinates(snapshot)
    if periodic
        return nb_periodic_kde_curve(
            xs;
            xmin=Float64(xmin),
            xmax=Float64(xmax),
            grid_points=Int(grid_points),
            bandwidth=Float64(bandwidth),
        )
    end
    return nb_density_curve(
        xs;
        nbins=Int(nbins),
        xmin=Float64(xmin),
        xmax=Float64(xmax),
        smoothing_window=9,
    )
end

function solve_periodic_onebody_fd(L::Real, D::Real, V1d; ngrid::Integer=600)
    n = Int(ngrid)
    n >= 64 || throw(ArgumentError("ngrid must be >= 64"))
    Lf = Float64(L)
    Df = Float64(D)
    dx = Lf / n
    xgrid = Float64[i * dx for i in 0:(n - 1)]

    H = zeros(Float64, n, n)
    kin_diag = 2.0 * Df / dx^2
    kin_off = -Df / dx^2
    @inbounds for i in 1:n
        H[i, i] = kin_diag + Float64(V1d(xgrid[i]))
    end
    @inbounds for i in 1:(n - 1)
        H[i, i + 1] = kin_off
        H[i + 1, i] = kin_off
    end
    H[1, n] = kin_off
    H[n, 1] = kin_off

    evals, evecs = eigen(Symmetric(H))
    return (
        x=xgrid,
        dx=dx,
        energies=Float64.(evals),
        vectors=Matrix{Float64}(evecs),
    )
end

function normalized_density(vec::AbstractVector{<:Real}, dx::Real)
    rho = abs2.(Float64.(vec))
    rho ./= sum(rho) * Float64(dx)
    return rho
end

function occupied_density(fd_ref, occupied::AbstractVector{<:Integer}; per_particle::Bool=true)
    density = zeros(Float64, length(fd_ref.x))
    for orb in occupied
        density .+= normalized_density(view(fd_ref.vectors, :, Int(orb)), fd_ref.dx)
    end
    if per_particle
        density ./= length(occupied)
    end
    return density
end

function save_history_figure(paths, case_title::AbstractString, variants, sims)
    fig = nb_plot_gfmc_history(
        sims;
        labels=[variant.label for variant in variants],
        colors=[variant.color for variant in variants],
        title_prefix=case_title,
    )
    path = joinpath(paths.figures_dir, "history.png")
    savefig(fig, path)
    return path
end

function save_density_figure(
    paths,
    case_title::AbstractString,
    variant,
    sim;
    periodic::Bool,
    xmin::Real,
    xmax::Real,
    bandwidth::Real=0.1,
    reference=nothing,
)
    snapshot = nb_last_snapshot(sim)
    centers, density = density_from_snapshot(
        snapshot;
        periodic=periodic,
        xmin=xmin,
        xmax=xmax,
        bandwidth=bandwidth,
    )
    fig = plot(
        centers,
        density;
        xlabel="x",
        ylabel="pooled one-body density",
        title="$case_title: final density",
        label=variant.label,
        color=variant.color,
        linewidth=2.4,
    )
    if reference !== nothing
        plot!(fig, reference.x, reference.y; label=reference.label, color=:black, linestyle=:dash, linewidth=2.2)
    end
    path = joinpath(paths.figures_dir, "density_$(variant.id).png")
    savefig(fig, path)
    return path
end

function save_pair_figure(paths, case_title::AbstractString, variant, sim, distance_fn; xmax::Real)
    rs = pair_separations(nb_last_snapshot(sim), distance_fn)
    isempty(rs) && return nothing
    centers, density = nb_density_curve(rs; nbins=100, xmin=0.0, xmax=Float64(xmax), smoothing_window=9)
    fig = plot(
        centers,
        density;
        xlabel="r",
        ylabel="pair-separation density",
        title="$case_title: final pair separations",
        label=variant.label,
        color=variant.color,
        linewidth=2.4,
    )
    path = joinpath(paths.figures_dir, "pair_$(variant.id).png")
    savefig(fig, path)
    return path
end

function save_sweep_figure(paths, case_title::AbstractString, variants, rows)
    sweep_rows = [row for row in rows if row.reference_kind != "none"]
    isempty(sweep_rows) && return nothing

    ids = unique(row.variant_family for row in sweep_rows)
    fig = plot(
        xscale=:log10,
        xlabel="GFMC time step",
        ylabel="post-equilibration mean energy",
        title="$case_title: sweep",
        legend=:best,
    )

    for family in ids
        family_rows = [row for row in sweep_rows if row.variant_family == family]
        sort!(family_rows, by=row -> row.dt)
        dts = Float64[row.dt for row in family_rows]
        means = Float64[row.mean_energy for row in family_rows]
        sems = Float64[row.sem_energy for row in family_rows]
        plot!(fig, dts, means; yerror=sems, marker=:circle, linewidth=2.2, label=family)
    end

    ref_vals = unique(filter(isfinite, Float64[row.reference_value for row in sweep_rows]))
    if length(ref_vals) == 1
        hline!(fig, [ref_vals[1]]; color=:black, linestyle=:dash, linewidth=1.8, label="reference")
    end

    path = joinpath(paths.figures_dir, "sweep.png")
    savefig(fig, path)
    return path
end

function run_direct_gfmc(H, initial_positions, params::GFMCParams;
    rng_seed::Integer=52,
    guiding::AbstractGuiding=NoGuiding(),
    nodepolicy::AbstractNodePolicy=NoNode(),
    reconfiguration::AbstractGFMCReconfiguration=SystematicReconfiguration(),
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout,
)
    sim = GFMCSim(
        H,
        params,
        initial_positions,
        MersenneTwister(Int(rng_seed));
        guiding=guiding,
        nodepolicy=nodepolicy,
        reconfiguration=reconfiguration,
    )
    return run_gfmc!(
        sim;
        snapshot_steps=snapshot_steps,
        show_progress=show_progress,
        progress_every=progress_every,
        progress_label=progress_label,
        debug=debug,
        debug_every=debug_every,
        debug_io=debug_io,
    )
end

function run_warmstart_gfmc(H, trial, initial_positions, gfmc_params::GFMCParams, vmc_params::VMCParams;
    vmc_seed::Integer=41,
    gfmc_seed::Integer=52,
    proposal::AbstractVMCProposal=DriftGaussianProposal(),
    guiding::AbstractGuiding=ImportanceGuiding(trial, H),
    nodepolicy::AbstractNodePolicy=NoNode(),
    reconfiguration::AbstractGFMCReconfiguration=SystematicReconfiguration(),
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    vmc_debug::Bool=false,
    vmc_debug_every::Integer=1,
    vmc_debug_io::IO=stdout,
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout,
)
    return run_gfmc_with_vmc_init(
        H,
        gfmc_params,
        initial_positions,
        trial,
        vmc_params;
        vmc_rng=MersenneTwister(Int(vmc_seed)),
        gfmc_rng=MersenneTwister(Int(gfmc_seed)),
        proposal=proposal,
        guiding=guiding,
        nodepolicy=nodepolicy,
        reconfiguration=reconfiguration,
        vmc_debug=vmc_debug,
        vmc_debug_every=vmc_debug_every,
        vmc_debug_io=vmc_debug_io,
        snapshot_steps=snapshot_steps,
        show_progress=show_progress,
        progress_every=progress_every,
        progress_label=progress_label,
        debug=debug,
        debug_every=debug_every,
        debug_io=debug_io,
    )
end

function run_spinless_fermion_warmstart(model, initial_positions, gfmc_params::GFMCParams, vmc_params::VMCParams;
    vmc_seed::Integer=41,
    gfmc_seed::Integer=52,
    proposal::AbstractVMCProposal=DriftGaussianProposal(),
    use_guiding::Bool=true,
    nodepolicy::AbstractNodePolicy=FixedNode(),
    reconfiguration::AbstractGFMCReconfiguration=SystematicReconfiguration(),
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    vmc_debug::Bool=false,
    vmc_debug_every::Integer=1,
    vmc_debug_io::IO=stdout,
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout,
)
    return run_gfmc_with_vmc_init(
        model,
        gfmc_params,
        initial_positions,
        vmc_params;
        vmc_rng=MersenneTwister(Int(vmc_seed)),
        gfmc_rng=MersenneTwister(Int(gfmc_seed)),
        proposal=proposal,
        use_guiding=use_guiding,
        nodepolicy=nodepolicy,
        reconfiguration=reconfiguration,
        vmc_debug=vmc_debug,
        vmc_debug_every=vmc_debug_every,
        vmc_debug_io=vmc_debug_io,
        snapshot_steps=snapshot_steps,
        show_progress=show_progress,
        progress_every=progress_every,
        progress_label=progress_label,
        debug=debug,
        debug_every=debug_every,
        debug_io=debug_io,
    )
end

function variant_record(;
    id::AbstractString,
    label::AbstractString,
    variant_family::AbstractString,
    color,
    theory_note::AbstractString,
    capability_note::AbstractString,
    gfmc_params::GFMCParams,
    vmc_params=nothing,
    reference_kind::AbstractString="none",
    reference_value::Real=NaN,
    reference_lower::Real=NaN,
    reference_upper::Real=NaN,
    density_reference=nothing,
    runner,
    metrics=(sim, extras) -> NamedTuple(),
    plotter=(paths, sim, extras) -> nothing,
)
    return (
        id=String(id),
        label=String(label),
        variant_family=String(variant_family),
        color=color,
        theory_note=String(theory_note),
        capability_note=String(capability_note),
        gfmc_params=gfmc_params,
        vmc_params=vmc_params,
        reference_kind=String(reference_kind),
        reference_value=Float64(reference_value),
        reference_lower=Float64(reference_lower),
        reference_upper=Float64(reference_upper),
        density_reference=density_reference,
        runner=runner,
        metrics=metrics,
        plotter=plotter,
    )
end

function tier_dt_values(tier::Symbol, values::Tuple)
    tier === :smoke && return [values[1]]
    tier === :final && return [values[end]]
    return collect(values)
end

function gaussian_ground_density(xgrid::AbstractVector{<:Real}, ω::Real)
    norm = sqrt(Float64(ω) / pi)
    return Float64[norm * exp(-Float64(ω) * Float64(x)^2) for x in xgrid]
end

function odd_ho_density(xgrid::AbstractVector{<:Real}, ω::Real)
    α = Float64(ω)
    prefac = 2.0 * α * sqrt(α / pi)
    return Float64[prefac * Float64(x)^2 * exp(-α * Float64(x)^2) for x in xgrid]
end

function build_free_particle_ring_case(tier::Symbol)
    L = 10.0
    bc = PeriodicBoundary1D(0.0, L)
    H = Hamiltonian(1, 0.5, _ -> 0.0, bc)
    targetN = tier === :smoke ? 128 : (tier === :sweep ? 512 : 1024)
    dt_values = tier_dt_values(tier, (2.0e-2, 1.0e-2, 5.0e-3))
    tau_total = tier === :smoke ? 0.8 : (tier === :sweep ? 1.2 : 2.0)

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=0.0,
            feedback=0.2,
            reconfiguration_interval=(tier === :final ? 10 : 5),
            energy_window=20,
            equil_fraction=0.2,
        )
        label = tier === :sweep ? @sprintf("dt=%.4f", dt) : "unguided"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="unguided",
            color=BENCHMARK_COLORS.red,
            theory_note="Exact ground state is uniform on the ring with E0 = 0.",
            capability_note="Exercises periodic wrapping with the generic unguided kernel and NoNode().",
            gfmc_params=params,
            reference_kind="exact_energy",
            reference_value=0.0,
            density_reference=(x=collect(range(0.0, L; length=400)), y=fill(1.0 / L, 400), label="uniform"),
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [[rand(rng_init) * L] for _ in 1:params.targetN]
                sim = run_direct_gfmc(
                    H,
                    initial_positions,
                    params;
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=(periodic=true, xmin=0.0, xmax=L, bandwidth=0.12))
            end,
            plotter=(paths, sim, extras) -> save_density_figure(
                paths,
                "Free particle ring GFMC",
                (id="density", label=label, color=BENCHMARK_COLORS.red),
                sim;
                periodic=extras.periodic,
                xmin=extras.xmin,
                xmax=extras.xmax,
                bandwidth=extras.bandwidth,
                reference=(x=collect(range(0.0, L; length=400)), y=fill(1.0 / L, 400), label="uniform"),
            ),
        )
    end

    return (
        id="free_particle_ring",
        title="Free Particle Ring 1D",
        benchmark_class=(tier === :smoke ? :accuracy_smoke : :accuracy),
        variants=variants,
    )
end

function build_harmonic_oscillator_unguided_case(tier::Symbol)
    ω = 1.0
    H = Hamiltonian(1, 0.5, R -> 0.5 * ω^2 * R[1]^2)
    targetN = tier === :smoke ? 192 : (tier === :sweep ? 768 : 2000)
    dt_values = tier_dt_values(tier, (1.0e-2, 5.0e-3, 2.5e-3))
    tau_total = tier === :smoke ? 0.8 : (tier === :sweep ? 1.5 : 2.5)

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=0.5 * ω,
            feedback=0.2,
            reconfiguration_interval=5,
            energy_window=30,
        )
        label = tier === :sweep ? @sprintf("dt=%.4f", dt) : "unguided"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="unguided",
            color=BENCHMARK_COLORS.red,
            theory_note="The one-dimensional harmonic oscillator has exact ground-state energy E0 = ω / 2.",
            capability_note="Open-boundary unguided reference for generic Hamiltonian propagation.",
            gfmc_params=params,
            reference_kind="exact_energy",
            reference_value=0.5 * ω,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [[2.0 * rand(rng_init) - 1.0] for _ in 1:params.targetN]
                sim = run_direct_gfmc(
                    H,
                    initial_positions,
                    params;
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=NamedTuple())
            end,
        )
    end

    return (
        id="harmonic_oscillator_unguided",
        title="Harmonic Oscillator 1D (Unguided)",
        benchmark_class=:accuracy,
        variants=variants,
    )
end

function build_harmonic_oscillator_guided_case(tier::Symbol)
    ω = 1.0
    H = Hamiltonian(1, 0.5, R -> 0.5 * ω^2 * R[1]^2)
    logpsi(R) = -0.5 * ω * R[1]^2
    gradlogpsi(R) = Float64[-ω * R[1]]
    lapllogpsi(R) = -ω
    trial = TrialWF(logpsi, gradlogpsi, lapllogpsi)
    guiding = ImportanceGuiding(trial, H)

    targetN = tier === :smoke ? 128 : (tier === :sweep ? 512 : 1500)
    dt_values = tier_dt_values(tier, (2.0e-2, 1.0e-2, 5.0e-3))
    tau_total = tier === :smoke ? 0.6 : (tier === :sweep ? 1.0 : 1.5)
    vmc_dt = tier === :smoke ? 2.0e-2 : 1.0e-2
    vmc_tau = tier === :final ? 0.8 : 0.4

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=0.5 * ω,
            feedback=0.05,
            reconfiguration_interval=(tier === :final ? 10 : 5),
            energy_window=20,
        )
        vmc_params = make_vmc_params(dt=vmc_dt, tau_total=vmc_tau, targetN=targetN, ET0=0.5 * ω)
        label = tier === :sweep ? @sprintf("guided dt=%.4f", dt) : "guided warm start"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="guided",
            color=BENCHMARK_COLORS.teal,
            theory_note="The exact Gaussian trial state makes this a near-zero-variance guided benchmark.",
            capability_note="Exercises importance sampling, VMC warm start, and low-variance population control.",
            gfmc_params=params,
            vmc_params=vmc_params,
            reference_kind="exact_energy",
            reference_value=0.5 * ω,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [[2.0 * rand(rng_init) - 1.0] for _ in 1:params.targetN]
                sim = run_warmstart_gfmc(
                    H,
                    trial,
                    initial_positions,
                    params,
                    vmc_params;
                    guiding=guiding,
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=NamedTuple())
            end,
        )
    end

    return (
        id="harmonic_oscillator_guided",
        title="Harmonic Oscillator 1D (Guided Exact Trial)",
        benchmark_class=:accuracy,
        variants=variants,
    )
end

function build_cosine_lattice_ring_case(tier::Symbol)
    M = 5
    a = 2.0
    L = M * a
    V0 = -5.0
    D = 0.5
    k_lat = 2pi / a
    λ_trial = -V0 / 2
    bc = PeriodicBoundary1D(0.0, L)
    H = Hamiltonian(1, D, R -> V0 * cos(k_lat * R[1]), bc)
    trial = TrialWF(
        R -> λ_trial * cos(k_lat * R[1]),
        R -> Float64[-λ_trial * k_lat * sin(k_lat * R[1])],
        R -> -λ_trial * k_lat^2 * cos(k_lat * R[1]),
    )
    guiding = ImportanceGuiding(trial, H)

    fd = solve_periodic_onebody_fd(L, D, x -> V0 * cos(k_lat * x); ngrid=700)
    exact_energy = fd.energies[1]
    exact_density = normalized_density(view(fd.vectors, :, 1), fd.dx)

    targetN = tier === :smoke ? 192 : (tier === :sweep ? 768 : 1500)
    dt_values = tier === :final ? [2.0e-3] : tier_dt_values(tier, (4.0e-3, 2.0e-3, 1.0e-3))
    tau_total = tier === :smoke ? 0.4 : (tier === :sweep ? 0.8 : 1.8)
    modes = tier === :final ? [("unguided", false, BENCHMARK_COLORS.red)] : [("unguided", false, BENCHMARK_COLORS.red), ("guided", true, BENCHMARK_COLORS.teal)]
    variants = NamedTuple[]

    for dt in dt_values, (family, use_guiding, color) in modes
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=exact_energy,
            feedback=0.15,
            reconfiguration_interval=5,
            energy_window=30,
        )
        label = tier === :sweep ? @sprintf("%s dt=%.4f", family, dt) : family
        push!(variants, variant_record(
            id=replace(label, " " => "_", "=" => "", "." => "p"),
            label=label,
            variant_family=family,
            color=color,
            theory_note="A single particle in a periodic cosine lattice has an exact one-body reference from the periodic finite-difference solver.",
            capability_note="Periodic external-potential benchmark comparing unguided and guided generic GFMC.",
            gfmc_params=params,
            reference_kind="exact_energy",
            reference_value=exact_energy,
            density_reference=(x=fd.x, y=exact_density, label="FD reference"),
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [[rand(rng_init) * L] for _ in 1:params.targetN]
                sim = run_direct_gfmc(
                    H,
                    initial_positions,
                    params;
                    guiding=(use_guiding ? guiding : NoGuiding()),
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=(periodic=true, xmin=0.0, xmax=L, bandwidth=0.12))
            end,
            plotter=(paths, sim, extras) -> save_density_figure(
                paths,
                "Cosine lattice ring GFMC",
                (id=replace(label, " " => "_", "=" => "", "." => "p"), label=label, color=color),
                sim;
                periodic=extras.periodic,
                xmin=extras.xmin,
                xmax=extras.xmax,
                bandwidth=extras.bandwidth,
                reference=(x=fd.x, y=exact_density, label="FD reference"),
            ),
        ))
    end

    return (
        id="cosine_lattice_ring",
        title="Cosine Lattice Ring 1D",
        benchmark_class=:accuracy,
        variants=variants,
    )
end

function _safe_nonzero_x(rng::AbstractRNG; x_min::Float64=1.0e-2)
    x = randn(rng)
    while abs(x) < x_min
        x = randn(rng)
    end
    return x
end

function build_hydrogen_unguided_case(tier::Symbol)
    H = Hamiltonian(1, 0.5, R -> -1 / abs(R[1]))
    targetN = tier === :smoke ? 128 : (tier === :sweep ? 384 : 1000)
    dt_values = tier_dt_values(tier, (2.0e-3, 1.0e-3, 5.0e-4))
    tau_total = tier === :smoke ? 0.08 : (tier === :sweep ? 0.2 : 0.35)
    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=-0.5,
            feedback=0.1,
            reconfiguration_interval=1,
            branch_cap=4.0,
            energy_window=20,
        )
        label = tier === :sweep ? @sprintf("dt=%.4f", dt) : "unguided singular"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="unguided",
            color=BENCHMARK_COLORS.red,
            theory_note="This singular 1D Coulomb-like problem is kept as a stability stress test rather than an accuracy regression.",
            capability_note="Exercises unguided branching near a singular potential.",
            gfmc_params=params,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [[_safe_nonzero_x(rng_init)] for _ in 1:params.targetN]
                sim = run_direct_gfmc(
                    H,
                    initial_positions,
                    params;
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=NamedTuple())
            end,
        )
    end
    return (id="hydrogen_unguided", title="Hydrogen 1D (Unguided Stress Test)", benchmark_class=:stress, variants=variants)
end

function build_hydrogen_fixed_node_case(tier::Symbol)
    α = 1.0
    H = Hamiltonian(1, 0.5, R -> -1 / abs(R[1]))
    logpsi(R) = begin
        x = R[1]
        ax = abs(x)
        ax == 0.0 ? -Inf : log(ax) - α * ax
    end
    gradlogpsi(R) = begin
        x = R[1]
        ax = abs(x)
        s = sign(x)
        ax == 0.0 ? Float64[0.0] : Float64[(s / ax) - α * s]
    end
    lapllogpsi(R) = begin
        x = R[1]
        ax = abs(x)
        ax == 0.0 ? -Inf : -1.0 / ax^2
    end
    signψ(R; tol::Float64=1.0e-10) = abs(R[1]) < tol ? 0.0 : sign(R[1])
    trial = TrialWF(logpsi, gradlogpsi, lapllogpsi, signψ)
    guiding = ImportanceGuiding(trial, H)

    targetN = tier === :smoke ? 128 : (tier === :sweep ? 384 : 800)
    dt_values = tier_dt_values(tier, (5.0e-3, 2.5e-3, 1.0e-3))
    tau_total = tier === :smoke ? 0.08 : (tier === :sweep ? 0.2 : 0.4)
    vmc_dt = 1.0e-3
    vmc_tau = 0.08
    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=-0.5,
            feedback=0.1,
            reconfiguration_interval=1,
            branch_cap=4.0,
            energy_window=20,
        )
        vmc_params = make_vmc_params(dt=vmc_dt, tau_total=vmc_tau, targetN=targetN, ET0=-0.5)
        label = tier === :sweep ? @sprintf("fixed-node dt=%.4f", dt) : "guided fixed node"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="fixed-node",
            color=BENCHMARK_COLORS.teal,
            theory_note="This odd-parity trial enforces a fixed node at x = 0 and is used as a singular fixed-node stress test.",
            capability_note="Exercises sign-aware importance sampling near a singular cusp.",
            gfmc_params=params,
            vmc_params=vmc_params,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [[_safe_nonzero_x(rng_init)] for _ in 1:params.targetN]
                sim = run_warmstart_gfmc(
                    H,
                    trial,
                    initial_positions,
                    params,
                    vmc_params;
                    guiding=guiding,
                    nodepolicy=FixedNode(),
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=NamedTuple())
            end,
        )
    end
    return (id="hydrogen_fixed_node", title="Hydrogen 1D (Fixed-Node Stress Test)", benchmark_class=:stress, variants=variants)
end

function build_two_particle_ho_guided_case(tier::Symbol)
    ω = 1.0
    κ = 0.7
    H = Hamiltonian(2, 0.5, R -> begin
        x1, x2 = R
        0.5 * ω^2 * (x1^2 + x2^2) + 0.5 * κ * (x1 - x2)^2
    end)
    ω_rel = sqrt(ω^2 + 2 * κ)
    trial = TrialWF(
        R -> begin
            x1, x2 = R
            S = x1 + x2
            Δ = x1 - x2
            -0.25 * ω * S^2 - 0.25 * ω_rel * Δ^2
        end,
        R -> begin
            x1, x2 = R
            S = x1 + x2
            Δ = x1 - x2
            Float64[
                -0.5 * (ω * S + ω_rel * Δ),
                -0.5 * (ω * S - ω_rel * Δ),
            ]
        end,
        _ -> -(ω + ω_rel),
    )
    guiding = ImportanceGuiding(trial, H)
    exact_energy = 0.5 * (ω + ω_rel)

    targetN = tier === :smoke ? 128 : (tier === :sweep ? 512 : 1200)
    dt_values = tier_dt_values(tier, (2.0e-2, 1.0e-2, 5.0e-3))
    tau_total = tier === :smoke ? 0.6 : (tier === :sweep ? 1.0 : 1.5)
    vmc_params = make_vmc_params(dt=1.0e-2, tau_total=0.4, targetN=targetN, ET0=exact_energy)

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=exact_energy,
            feedback=0.05,
            reconfiguration_interval=5,
            energy_window=20,
        )
        label = tier === :sweep ? @sprintf("guided dt=%.4f", dt) : "guided warm start"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="guided",
            color=BENCHMARK_COLORS.navy,
            theory_note="The coupled two-particle oscillator has an exact Gaussian ground state in center-of-mass and relative coordinates.",
            capability_note="Exercises multi-coordinate guided GFMC with an exact trial and warm start.",
            gfmc_params=params,
            vmc_params=vmc_params,
            reference_kind="exact_energy",
            reference_value=exact_energy,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [2 .* rand(rng_init, 2) .- 1 for _ in 1:params.targetN]
                sim = run_warmstart_gfmc(
                    H,
                    trial,
                    initial_positions,
                    params,
                    vmc_params;
                    guiding=guiding,
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=NamedTuple())
            end,
        )
    end
    return (id="two_particle_ho_guided", title="Coupled Two-Particle HO 1D", benchmark_class=:accuracy, variants=variants)
end

function build_two_particle_ho_fixed_node_case(tier::Symbol)
    ω = 1.0
    H = Hamiltonian(1, 0.5, R -> 0.5 * ω^2 * R[1]^2)
    trial = TrialWF(
        R -> begin
            x = R[1]
            ax = abs(x)
            ax == 0.0 ? -Inf : log(ax) - 0.5 * ω * x^2
        end,
        R -> begin
            x = R[1]
            ax = abs(x)
            s = sign(x)
            ax == 0.0 ? Float64[0.0] : Float64[(s / ax) - ω * x]
        end,
        R -> begin
            x = R[1]
            ax = abs(x)
            ax == 0.0 ? -Inf : -1.0 / ax^2 - ω
        end,
        R -> abs(R[1]) < 1.0e-10 ? 0.0 : sign(R[1]),
    )
    guiding = ImportanceGuiding(trial, H)
    exact_energy = 1.5 * ω

    targetN = tier === :smoke ? 128 : (tier === :sweep ? 512 : 1000)
    dt_values = tier_dt_values(tier, (1.0e-2, 5.0e-3, 2.5e-3))
    tau_total = tier === :smoke ? 0.8 : (tier === :sweep ? 1.0 : 1.5)
    vmc_params = make_vmc_params(dt=1.0e-2, tau_total=0.4, targetN=targetN, ET0=exact_energy)

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=exact_energy,
            feedback=0.1,
            reconfiguration_interval=5,
            energy_window=20,
        )
        label = tier === :sweep ? @sprintf("fixed-node dt=%.4f", dt) : "guided fixed node"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="fixed-node",
            color=BENCHMARK_COLORS.teal,
            theory_note="This is the exact odd-parity harmonic-oscillator sector, so fixed-node GFMC should recover E = 3ω / 2.",
            capability_note="Exercises exact fixed-node handling with a nodal surface at x = 0.",
            gfmc_params=params,
            vmc_params=vmc_params,
            reference_kind="exact_energy",
            reference_value=exact_energy,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = [[_safe_nonzero_x(rng_init)] for _ in 1:params.targetN]
                sim = run_warmstart_gfmc(
                    H,
                    trial,
                    initial_positions,
                    params,
                    vmc_params;
                    guiding=guiding,
                    nodepolicy=FixedNode(),
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=NamedTuple())
            end,
        )
    end
    return (id="two_particle_ho_fixed_node", title="Odd-Parity HO Fixed-Node Benchmark", benchmark_class=:accuracy, variants=variants)
end

function build_fermion_ring_fixed_node_case(tier::Symbol)
    N = 3
    M = 4
    a = 1.0
    L = M * a
    V0 = -1.0
    model = SpinlessFermionRing1D(N, a, L, V0; D=0.5, node_tol=1.0e-7, trig_eps=1.0e-10)
    fd = solve_periodic_onebody_fd(L, model.D, x -> V0 * cos(model.k_lat * x); ngrid=700)
    exact_energy = sum(@view fd.energies[1:N])
    exact_density = occupied_density(fd, collect(1:N); per_particle=true)

    targetN = tier === :smoke ? 160 : (tier === :sweep ? 512 : 1200)
    dt_values = tier_dt_values(tier, (4.0e-3, 2.0e-3, 1.0e-3))
    tau_total = tier === :smoke ? 0.3 : (tier === :sweep ? 0.8 : 1.5)
    vmc_dt = tier === :smoke ? 2.0e-2 : 1.0e-2
    vmc_tau = 0.6

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=exact_energy,
            feedback=0.08,
            reconfiguration_interval=2,
            energy_window=20,
        )
        vmc_params = make_vmc_params(dt=vmc_dt, tau_total=vmc_tau, targetN=targetN, ET0=exact_energy)
        label = tier === :sweep ? @sprintf("fixed-node dt=%.4f", dt) : "specialized fixed node"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="fixed-node",
            color=BENCHMARK_COLORS.navy,
            theory_note="This case compares the specialized fermion-ring kernel against the exact many-body energy of the underlying one-body cosine problem. It is currently best treated as a specialized fixed-node diagnostic rather than a locked accuracy benchmark.",
            capability_note="Exercises the specialized fermion-ring kernel, fixed-node handling, and VMC-backed initialization.",
            gfmc_params=params,
            vmc_params=vmc_params,
            reference_kind="exact_energy",
            reference_value=exact_energy,
            density_reference=(x=fd.x, y=exact_density, label="exact pooled density"),
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = sample_uniform_configurations(model, params.targetN, rng_init)
                sim = run_spinless_fermion_warmstart(
                    model,
                    initial_positions,
                    params,
                    vmc_params;
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=(periodic=true, xmin=0.0, xmax=L, bandwidth=0.10 * a, distance_fn=(x1, x2) -> distance_1d(model.bc, x1, x2)))
            end,
            metrics=(sim, extras) -> begin
                rs = pair_separations(nb_last_snapshot(sim), extras.distance_fn)
                return (
                    min_pair_separation=isempty(rs) ? NaN : minimum(rs),
                    pair_samples=length(rs),
                )
            end,
            plotter=(paths, sim, extras) -> begin
                save_density_figure(
                    paths,
                    "Spinless fermion ring GFMC",
                    (id=replace(label, " " => "_", "=" => "", "." => "p"), label=label, color=BENCHMARK_COLORS.navy),
                    sim;
                    periodic=extras.periodic,
                    xmin=extras.xmin,
                    xmax=extras.xmax,
                    bandwidth=extras.bandwidth,
                    reference=(x=fd.x, y=exact_density, label="exact pooled density"),
                )
                save_pair_figure(paths, "Spinless fermion ring GFMC", (id="pair", label=label, color=BENCHMARK_COLORS.navy), sim, extras.distance_fn; xmax=0.5 * L)
            end,
        )
    end
    return (id="fermion_ring_fixed_node", title="Spinless Fermion Ring 1D", benchmark_class=:stress, variants=variants)
end

function build_hardcore_boson_ring_case(tier::Symbol)
    N = 3
    M = 4
    a = 1.0
    L = M * a
    V0 = -1.5
    D = 0.5
    g_coul = 0.8
    coulomb_softening = 0.05 * a
    r_core = 0.18 * a
    hard_core_barrier = 1000.0
    λ_trial = -0.5 * V0
    pair_strength = 0.75
    pair_eps = 1.0e-8
    bc = PeriodicBoundary1D(0.0, L)
    k_lat = 2pi / a
    pair_alpha = pi / L

    pair_distance_fn(xi, xj) = abs(displacement(bc, xj, xi))
    pair_potential(r) = (r < r_core ? hard_core_barrier : 0.0) + g_coul / sqrt(r^2 + coulomb_softening^2)
    H = Hamiltonian(N, D, R -> begin
        total = 0.0
        @inbounds for i in 1:N
            total += V0 * cos(k_lat * R[i])
        end
        @inbounds for i in 1:(N - 1)
            for j in (i + 1):N
                total += pair_potential(pair_distance_fn(R[i], R[j]))
            end
        end
        total
    end, bc)

    function boson_trial_terms(R)
        logabs = 0.0
        grad = zeros(Float64, N)
        lapl = 0.0
        @inbounds for i in 1:N
            xi = wrap_coordinate(bc, R[i])
            coskx = cos(k_lat * xi)
            sinkx = sin(k_lat * xi)
            logabs += λ_trial * coskx
            grad[i] = -λ_trial * k_lat * sinkx
            lapl += -λ_trial * k_lat^2 * coskx
        end
        @inbounds for i in 1:(N - 1)
            for j in (i + 1):N
                dx_ij = displacement(bc, R[j], R[i])
                u = pair_alpha * dx_ij
                s = sin(u)
                s_abs = abs(s)
                s_eff = s_abs < pair_eps ? pair_eps : s_abs
                cot_u = cos(u) / (s >= 0 ? s_eff : -s_eff)
                logabs += pair_strength * log(s_eff)
                c = pair_strength * pair_alpha * cot_u
                grad[i] += c
                grad[j] -= c
                lapl += -2.0 * pair_strength * pair_alpha^2 / (s_eff * s_eff)
            end
        end
        return logabs, grad, lapl
    end

    trial = TrialWF(R -> boson_trial_terms(R)[1], R -> boson_trial_terms(R)[2], R -> boson_trial_terms(R)[3], _ -> 1.0)
    guiding = ImportanceGuiding(trial, H)

    function valid_hardcore_configuration(R)
        @inbounds for i in 1:(length(R) - 1)
            for j in (i + 1):length(R)
                pair_distance_fn(R[i], R[j]) > r_core || return false
            end
        end
        return true
    end

    function sample_hardcore_ring_configurations(nwalkers::Integer, rng::AbstractRNG; max_tries::Integer=10000)
        X = Matrix{Float64}(undef, N, Int(nwalkers))
        @inbounds for w in 1:Int(nwalkers)
            accepted = false
            for _ in 1:Int(max_tries)
                R = rand(rng, N) .* L
                if valid_hardcore_configuration(R)
                    X[:, w] = R
                    accepted = true
                    break
                end
            end
            accepted || error("Failed to sample a valid hard-core ring configuration after $(max_tries) tries")
        end
        return X
    end

    targetN = tier === :smoke ? 128 : (tier === :sweep ? 384 : 800)
    dt_values = tier_dt_values(tier, (1.0e-3, 5.0e-4, 2.5e-4))
    tau_total = tier === :smoke ? 0.08 : (tier === :sweep ? 0.2 : 0.4)
    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=0.0,
            feedback=0.1,
            reconfiguration_interval=2,
            branch_cap=4.0,
            energy_window=20,
        )
        label = tier === :sweep ? @sprintf("dt=%.4f", dt) : "guided hard-core bosons"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="jastrow",
            color=BENCHMARK_COLORS.teal,
            theory_note="This case is kept as a correlated periodic stress test: the key observable is suppressed short-range pair weight rather than an exact energy.",
            capability_note="Exercises periodic pair interactions, hard-core initialization, and bosonic guiding.",
            gfmc_params=params,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = sample_hardcore_ring_configurations(params.targetN, rng_init)
                sim = run_direct_gfmc(
                    H,
                    initial_positions,
                    params;
                    guiding=guiding,
                    nodepolicy=NoNode(),
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=(distance_fn=pair_distance_fn, pair_cutoff=r_core))
            end,
            metrics=(sim, extras) -> begin
                rs = pair_separations(nb_last_snapshot(sim), extras.distance_fn)
                return (
                    min_pair_separation=isempty(rs) ? NaN : minimum(rs),
                    count_below_core=count(r -> r < extras.pair_cutoff, rs),
                )
            end,
            plotter=(paths, sim, extras) -> save_pair_figure(
                paths,
                "Hard-core boson ring GFMC",
                (id=replace(label, " " => "_", "=" => "", "." => "p"), label=label, color=BENCHMARK_COLORS.teal),
                sim,
                extras.distance_fn;
                xmax=0.5 * L,
            ),
        )
    end
    return (id="hardcore_boson_ring", title="Hard-Core Boson Ring 1D", benchmark_class=:stress, variants=variants)
end

function build_periodic_ion_single_particle_case(tier::Symbol)
    M = 3
    a = 1.0
    ring = build_periodic_ion_ring_model(M, a; D=0.5, ion_strength=1.0, ion_softening=0.35 * a, kmax=6, quad_points=1536)
    H = periodic_ion_hamiltonian(ring, 1)
    trial = single_particle_trial_wavefunction(ring; orbital_index=1)
    guiding = ImportanceGuiding(trial, H)
    exact_energy = ring.energies[1]
    xgrid = Float64[i * (ring.L / 600) for i in 0:599]
    phi_curve, _, _ = orbital_curve(ring, 1, xgrid)
    density_ref = phi_curve .^ 2
    density_ref ./= sum(density_ref) * (ring.L / length(xgrid))

    targetN = tier === :smoke ? 128 : (tier === :sweep ? 512 : 1024)
    dt_values = tier_dt_values(tier, (6.0e-3, 3.0e-3, 1.5e-3))
    tau_total = tier === :smoke ? 0.3 : (tier === :sweep ? 0.6 : 1.2)
    vmc_dt = tier === :smoke ? 2.0e-2 : 1.0e-2
    vmc_tau = 0.6

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=exact_energy,
            feedback=0.05,
            reconfiguration_interval=5,
            energy_window=20,
        )
        vmc_params = make_vmc_params(dt=vmc_dt, tau_total=vmc_tau, targetN=targetN, ET0=exact_energy)
        label = tier === :sweep ? @sprintf("guided dt=%.4f", dt) : "guided warm-started single particle"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="guided",
            color=BENCHMARK_COLORS.navy,
            theory_note="The lowest orbital of the periodic ion-ring one-body Hamiltonian is exact, so guided GFMC should be near zero variance.",
            capability_note="Exercises generic Hamiltonian + exact guiding + VMC warm start on a periodic soft-Coulomb ion lattice.",
            gfmc_params=params,
            vmc_params=vmc_params,
            reference_kind="exact_energy",
            reference_value=exact_energy,
            density_reference=(x=xgrid, y=density_ref, label="exact orbital density"),
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = sample_uniform_ring_configurations(1, ring.L, params.targetN, rng_init)
                sim = run_warmstart_gfmc(
                    H,
                    trial,
                    initial_positions,
                    params,
                    vmc_params;
                    guiding=guiding,
                    nodepolicy=NoNode(),
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=(periodic=true, xmin=0.0, xmax=ring.L, bandwidth=0.08 * a))
            end,
            plotter=(paths, sim, extras) -> save_density_figure(
                paths,
                "Periodic ion ring GFMC (single particle)",
                (id=replace(label, " " => "_", "=" => "", "." => "p"), label=label, color=BENCHMARK_COLORS.navy),
                sim;
                periodic=extras.periodic,
                xmin=extras.xmin,
                xmax=extras.xmax,
                bandwidth=extras.bandwidth,
                reference=(x=xgrid, y=density_ref, label="exact orbital density"),
            ),
        )
    end
    return (id="periodic_ion_single_particle", title="Periodic Ion Ring 1D (Single Particle)", benchmark_class=:accuracy, variants=variants)
end

function build_periodic_ion_spinless_fermions_case(tier::Symbol)
    N = 3
    M = 3
    a = 1.0
    ring = build_periodic_ion_ring_model(M, a; D=0.5, ion_strength=1.0, ion_softening=0.35 * a, kmax=max(6, N + 2), quad_points=1536)
    H = periodic_ion_hamiltonian(ring, N)
    trial = fermion_determinant_trial_wavefunction(ring, N; node_tol=1.0e-10)
    guiding = ImportanceGuiding(trial, H)
    exact_energy = exact_manybody_energy(ring, N)
    xgrid = Float64[i * (ring.L / 600) for i in 0:599]
    density_ref = occupied_pooled_density(ring, collect(1:N), xgrid; per_particle=true)
    density_ref ./= sum(density_ref) * (ring.L / length(xgrid))

    targetN = tier === :smoke ? 128 : (tier === :sweep ? 512 : 1024)
    dt_values = tier_dt_values(tier, (4.0e-3, 2.0e-3, 1.0e-3))
    tau_total = tier === :smoke ? 0.3 : (tier === :sweep ? 0.7 : 1.2)
    vmc_dt = 3.0e-3
    vmc_tau = 0.24

    variants = map(enumerate(dt_values)) do (k, dt)
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=exact_energy,
            feedback=0.05,
            reconfiguration_interval=2,
            energy_window=20,
        )
        vmc_params = make_vmc_params(dt=vmc_dt, tau_total=vmc_tau, targetN=targetN, ET0=exact_energy)
        label = tier === :sweep ? @sprintf("fixed-node dt=%.4f", dt) : "exact-determinant fixed node"
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family="fixed-node",
            color=BENCHMARK_COLORS.navy,
            theory_note="The determinant of the lowest periodic ion-ring orbitals is the exact fermionic ground state, so fixed-node GFMC should reproduce the exact many-body energy.",
            capability_note="Exercises exact determinant guiding, fixed-node sign handling, and VMC warm start on the generic kernel path.",
            gfmc_params=params,
            vmc_params=vmc_params,
            reference_kind="exact_energy",
            reference_value=exact_energy,
            density_reference=(x=xgrid, y=density_ref, label="exact pooled density"),
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = sample_uniform_ring_configurations(N, ring.L, params.targetN, rng_init; min_separation=0.02 * a)
                sim = run_warmstart_gfmc(
                    H,
                    trial,
                    initial_positions,
                    params,
                    vmc_params;
                    guiding=guiding,
                    nodepolicy=FixedNode(),
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=(periodic=true, xmin=0.0, xmax=ring.L, bandwidth=0.08 * a, distance_fn=(x1, x2) -> distance_1d(ring.bc, x1, x2)))
            end,
            metrics=(sim, extras) -> begin
                rs = pair_separations(nb_last_snapshot(sim), extras.distance_fn)
                return (
                    min_pair_separation=isempty(rs) ? NaN : minimum(rs),
                    count_close_pairs=count(r -> r < 0.02 * a, rs),
                )
            end,
            plotter=(paths, sim, extras) -> begin
                save_density_figure(
                    paths,
                    "Periodic ion ring GFMC (spinless fermions)",
                    (id=replace(label, " " => "_", "=" => "", "." => "p"), label=label, color=BENCHMARK_COLORS.navy),
                    sim;
                    periodic=extras.periodic,
                    xmin=extras.xmin,
                    xmax=extras.xmax,
                    bandwidth=extras.bandwidth,
                    reference=(x=xgrid, y=density_ref, label="exact pooled density"),
                )
                save_pair_figure(paths, "Periodic ion ring GFMC (spinless fermions)", (id="pair", label=label, color=BENCHMARK_COLORS.navy), sim, extras.distance_fn; xmax=0.5 * ring.L)
            end,
        )
    end
    return (id="periodic_ion_spinless_fermions", title="Periodic Ion Ring 1D (Spinless Fermions)", benchmark_class=:accuracy, variants=variants)
end

function _periodic_ion_boson_case(tier::Symbol; tg_scaffold::Bool)
    N = 3
    M = 3
    a = 1.0
    problem = build_periodic_ion_ring_boson_problem(
        N,
        M,
        a;
        L=M * a,
        D=0.5,
        ion_strength=1.0,
        ion_softening=0.35 * a,
        bb_strength=0.45,
        bb_softening=0.25 * a,
        hard_core_radius=(tg_scaffold ? 0.0 : 0.05 * a),
        hard_core_barrier=(tg_scaffold ? 0.0 : 12.0),
        pair_metric=:chord,
        kmax=max(6, N + 2),
        quad_points=1536,
    )
    H = bosonic_hamiltonian(problem)
    nonint_ref = noninteracting_boson_energy(problem.ring, N)
    tg_ref = tg_scaffold_energy(problem.ring, N)
    xgrid = Float64[i * (problem.ring.L / 600) for i in 0:599]
    density_ref = if tg_scaffold
        rho = occupied_pooled_density(problem.ring, collect(1:N), xgrid; per_particle=true)
        rho ./= sum(rho) * (problem.ring.L / length(xgrid))
        rho
    else
        phi_curve, _, _ = orbital_curve(problem.ring, 1, xgrid)
        rho = phi_curve .^ 2
        rho ./= sum(rho) * (problem.ring.L / length(xgrid))
        rho
    end
    trial = if tg_scaffold
        bosonic_tg_scaffold_trial_wavefunction(
            problem;
            node_tol=1.0e-10,
            smooth_pair_strength=0.04,
            smooth_pair_softening=max(problem.bb_softening, 0.25 * a),
        )
    else
        bosonic_orbital_jastrow_trial_wavefunction(
            problem;
            orbital_index=1,
            pair_power=1.0,
            pair_eps=1.0e-8,
            smooth_pair_strength=0.08,
            smooth_pair_softening=max(problem.bb_softening, 0.25 * a),
        )
    end
    guiding = ImportanceGuiding(trial, H)
    targetN = tier === :smoke ? 128 : (tier === :sweep ? 384 : 768)
    dt_values = tier_dt_values(tier, (2.0e-3, 1.0e-3, 5.0e-4))
    tau_total = tier === :smoke ? 0.2 : (tier === :sweep ? 0.4 : 0.8)
    vmc_dt = 2.0e-3
    vmc_tau = 0.16
    family = tg_scaffold ? "tg-scaffold" : "orbital-jastrow"
    color = tg_scaffold ? BENCHMARK_COLORS.navy : BENCHMARK_COLORS.teal
    title = tg_scaffold ? "Periodic Ion Ring 1D (Bosons, TG Scaffold)" : "Periodic Ion Ring 1D (Bosons, Jastrow)"
    theory_note = tg_scaffold ?
        "For repulsive bosons the energy should lie between the noninteracting boson energy and the Tonks-Girardeau scaffold energy; the |det| scaffold targets the strong-correlation side of that bracket." :
        "For repulsive bosons the energy should lie between the noninteracting boson energy and the Tonks-Girardeau scaffold energy; the orbital-Jastrow trial targets the weak-to-moderate correlation side of that bracket."

    variants = map(enumerate(dt_values)) do (k, dt)
        ET0 = tg_scaffold ? tg_ref : nonint_ref
        params = make_gfmc_params(
            dt=dt,
            tau_total=tau_total,
            targetN=targetN,
            ET0=ET0,
            feedback=0.05,
            reconfiguration_interval=2,
            energy_window=20,
        )
        vmc_params = make_vmc_params(dt=vmc_dt, tau_total=vmc_tau, targetN=targetN, ET0=ET0)
        label = tier === :sweep ? @sprintf("%s dt=%.4f", family, dt) : family
        variant_record(
            id=tier === :sweep ? "dt$(k)" : "main",
            label=label,
            variant_family=family,
            color=color,
            theory_note=theory_note * " The noninteracting boson energy and TG scaffold energy are reported as context, not as strict bounds, because the smooth boson-boson repulsion can push the exact energy above the pure TG scaffold.",
            capability_note="Exercises generic periodic bosonic GFMC with pair interactions, guiding, and VMC warm start.",
            gfmc_params=params,
            vmc_params=vmc_params,
            reference_kind="theory_context",
            reference_lower=nonint_ref,
            reference_upper=tg_ref,
            runner=(; debug, debug_every, debug_io, show_progress, progress_every) -> begin
                rng_init = MersenneTwister(1234)
                initial_positions = sample_boson_ring_configurations(
                    problem,
                    params.targetN,
                    rng_init;
                    min_separation=(tg_scaffold ? 0.04 * a : problem.hard_core_radius),
                )
                sim = run_warmstart_gfmc(
                    H,
                    trial,
                    initial_positions,
                    params,
                    vmc_params;
                    guiding=guiding,
                    nodepolicy=NoNode(),
                    snapshot_steps=nb_default_snapshot_steps(params.nsteps),
                    show_progress=show_progress,
                    progress_every=progress_every,
                    progress_label=label,
                    vmc_debug=debug,
                    vmc_debug_every=debug_every,
                    vmc_debug_io=debug_io,
                    debug=debug,
                    debug_every=debug_every,
                    debug_io=debug_io,
                )
                return (sim=sim, extras=(periodic=true, xmin=0.0, xmax=problem.ring.L, bandwidth=0.08 * a, distance_fn=(x1, x2) -> pair_distance(problem, x1, x2)))
            end,
            metrics=(sim, extras) -> begin
                rs = pair_separations(nb_last_snapshot(sim), extras.distance_fn)
                return (
                    min_pair_separation=isempty(rs) ? NaN : minimum(rs),
                    noninteracting_reference=nonint_ref,
                    tg_reference=tg_ref,
                )
            end,
            plotter=(paths, sim, extras) -> begin
                save_density_figure(
                    paths,
                    title,
                    (id=replace(label, " " => "_", "=" => "", "." => "p"), label=label, color=color),
                    sim;
                    periodic=extras.periodic,
                    xmin=extras.xmin,
                    xmax=extras.xmax,
                    bandwidth=extras.bandwidth,
                    reference=(x=xgrid, y=density_ref, label=(tg_scaffold ? "TG pooled density" : "orbital density")),
                )
                save_pair_figure(paths, title, (id="pair", label=label, color=color), sim, extras.distance_fn; xmax=0.5 * problem.ring.L)
            end,
        )
    end

    return (
        id=(tg_scaffold ? "periodic_ion_bosons_tg" : "periodic_ion_bosons_jastrow"),
        title=title,
        benchmark_class=:stress,
        variants=variants,
    )
end

build_periodic_ion_bosons_jastrow_case(tier::Symbol) = _periodic_ion_boson_case(tier; tg_scaffold=false)
build_periodic_ion_bosons_tg_case(tier::Symbol) = _periodic_ion_boson_case(tier; tg_scaffold=true)

function benchmark_case(case_id::AbstractString, tier::Symbol)
    case_id == "free_particle_ring" && return build_free_particle_ring_case(tier)
    case_id == "harmonic_oscillator_unguided" && return build_harmonic_oscillator_unguided_case(tier)
    case_id == "harmonic_oscillator_guided" && return build_harmonic_oscillator_guided_case(tier)
    case_id == "cosine_lattice_ring" && return build_cosine_lattice_ring_case(tier)
    case_id == "hydrogen_unguided" && return build_hydrogen_unguided_case(tier)
    case_id == "hydrogen_fixed_node" && return build_hydrogen_fixed_node_case(tier)
    case_id == "two_particle_ho_guided" && return build_two_particle_ho_guided_case(tier)
    case_id == "two_particle_ho_fixed_node" && return build_two_particle_ho_fixed_node_case(tier)
    case_id == "fermion_ring_fixed_node" && return build_fermion_ring_fixed_node_case(tier)
    case_id == "hardcore_boson_ring" && return build_hardcore_boson_ring_case(tier)
    case_id == "periodic_ion_single_particle" && return build_periodic_ion_single_particle_case(tier)
    case_id == "periodic_ion_spinless_fermions" && return build_periodic_ion_spinless_fermions_case(tier)
    case_id == "periodic_ion_bosons_jastrow" && return build_periodic_ion_bosons_jastrow_case(tier)
    case_id == "periodic_ion_bosons_tg" && return build_periodic_ion_bosons_tg_case(tier)
    error("Unknown benchmark case id: $case_id")
end

function benchmark_summary_row(case_id::AbstractString, case_title::AbstractString, tier::Symbol, variant, sim, extra_metrics::NamedTuple)
    stats = benchmark_summary(sim)
    return merge((
        case_id=case_id,
        case_title=case_title,
        tier=String(tier),
        variant_id=variant.id,
        variant_label=variant.label,
        variant_family=variant.variant_family,
        dt=sim.params.dt,
        nsteps=sim.params.nsteps,
        nequil=sim.params.nequil,
        targetN=sim.params.targetN,
        reconfiguration_interval=sim.params.reconfiguration_interval,
        reference_kind=variant.reference_kind,
        reference_value=variant.reference_value,
        reference_lower=variant.reference_lower,
        reference_upper=variant.reference_upper,
        mean_energy=stats.mean_energy,
        sem_energy=stats.sem_energy,
        energy_error=(isfinite(variant.reference_value) ? (stats.mean_energy - variant.reference_value) : NaN),
        final_ET=stats.final_ET,
        final_population=stats.final_population,
        final_mean_weight=stats.final_mean_weight,
        final_effective_population=stats.final_effective_population,
        final_acceptance=stats.final_acceptance,
        final_variance=stats.final_variance,
        theory_note=variant.theory_note,
        capability_note=variant.capability_note,
    ), extra_metrics)
end

function run_benchmark_case(case_id::AbstractString, tier::Symbol)
    case = benchmark_case(case_id, tier)
    paths = benchmark_paths(case.id, tier)
    rows = NamedTuple[]
    sims = Any[]
    density_figure_paths = String[]

    debug_default = tier === :smoke
    debug_every = tier === :smoke ? 10 : 25
    show_progress = false
    progress_every = 0

    for variant in case.variants
        debug_log_path = joinpath(paths.logs_dir, "$(variant.id)_debug.log")
        open(debug_log_path, "w") do debug_io
            result = variant.runner(
                ;
                debug=debug_default,
                debug_every=debug_every,
                debug_io=debug_io,
                show_progress=show_progress,
                progress_every=progress_every,
            )
            sim = result.sim
            extras = result.extras
            push!(sims, sim)
            extra_metrics = variant.metrics(sim, extras)
            push!(rows, benchmark_summary_row(case.id, case.title, tier, variant, sim, extra_metrics))

            history_csv = joinpath(paths.tables_dir, "history_$(variant.id).csv")
            nb_write_csv(history_csv, nb_gfmc_rows(variant.label, sim))

            fig_path = variant.plotter(paths, sim, extras)
            fig_path === nothing || push!(density_figure_paths, fig_path)
        end
    end

    summary_csv = joinpath(paths.tables_dir, "summary.csv")
    nb_write_csv(summary_csv, rows)
    history_figure = save_history_figure(paths, case.title, case.variants, sims)
    sweep_figure = tier === :sweep ? save_sweep_figure(paths, case.title, case.variants, rows) : nothing

    return (
        case_id=case.id,
        case_title=case.title,
        summary_csv=summary_csv,
        history_figure=history_figure,
        sweep_figure=sweep_figure,
        rows=rows,
        density_figures=density_figure_paths,
    )
end

function run_benchmark_suite(case_ids::AbstractVector{<:AbstractString}, tier::Symbol; suite_name::AbstractString="default")
    results = [run_benchmark_case(case_id, tier) for case_id in case_ids]
    all_rows = reduce(vcat, [result.rows for result in results]; init=NamedTuple[])
    suite_output = suite_paths(suite_name, tier)
    summary_csv = joinpath(suite_output.tables_dir, "summary.csv")
    nb_write_csv(summary_csv, all_rows)

    exact_rows = [row for row in all_rows if row.reference_kind == "exact_energy" && isfinite(row.energy_error)]
    summary_figure = nothing
    if !isempty(exact_rows)
        ordered = sort(exact_rows, by=row -> abs(row.energy_error), rev=true)
        labels = [row.case_id * ":" * row.variant_family for row in ordered]
        errs = [row.energy_error for row in ordered]
        fig = bar(
            labels,
            errs;
            xlabel="benchmark case",
            ylabel="mean energy error",
            title="GFMC benchmark summary ($(tier))",
            legend=false,
            xrotation=45,
            size=(1500, 650),
        )
        hline!(fig, [0.0]; color=:black, linestyle=:dash, linewidth=1.5)
        summary_figure = joinpath(suite_output.figures_dir, "summary_energy_error.png")
        savefig(fig, summary_figure)
    end

    return (
        results=results,
        summary_csv=summary_csv,
        summary_figure=summary_figure,
    )
end
