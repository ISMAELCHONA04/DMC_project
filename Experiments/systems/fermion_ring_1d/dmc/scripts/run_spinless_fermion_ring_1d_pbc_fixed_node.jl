#=
Spinless-fermion DMC on a 1D periodic ring with a cosine lattice.

Physical model:
    H = -D * sum_i d^2/dx_i^2 + sum_i V0 * cos(2*pi*x_i/a),   D = 0.5
for N spinless fermions on [0, L) with periodic boundaries.

Trial wavefunction:
    Psi_T(R) = A(R) * G(R)
    A(R) = prod_{i<j} sin(pi * Delta x_ij / L)          (antisymmetric nodal factor)
    G(R) = exp(lambda * sum_i cos(2*pi*x_i/a))          (positive guiding factor)
with Delta x_ij from minimum-image signed displacement.

Fermionic nodes are enforced through FixedNode() + TrialWF.signpsi, where signpsi
returns 0.0 near nodes (|sin| < node_tol), otherwise +/-1.

This script produces:
1) console summary (parameters, ET estimate, SEM, variance, population stats),
2) ET history plot (smoothed + clipped),
3) mean-energy estimator history plot,
4) variance history plot,
5) final pooled one-body density vs trial one-body shape (+ N=1 exact reference),
6) smoothed snapshot density evolution,
7) external potential plot,
8) one-body drift-field illustration.
=#

using Random
using Statistics
using Printf
using LinearAlgebra

# Use a headless GR backend by default for cloud/batch runs.
if !haskey(ENV, "GKSwstype")
    ENV["GKSwstype"] = "100"
end

using Plots
using Logging

if !isdefined(Main, :System1D)
    include("../../../../../src/System1D.jl")
end

using .System1D: Hamiltonian, TrialWF, DMCParams, ImportanceGuiding, run_dmc
using .System1D: PeriodicBoundary1D, FixedNode, NoNode, NoGuiding
using .System1D: cell_length, wrap_coordinate, wrap_position, displacement
using .System1D: local_energy

# --------------------------------------------------------------------------
# User parameters 
# --------------------------------------------------------------------------
function parse_env_bool(key::AbstractString, default::Bool)
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return default
    low = lowercase(raw)
    if low in ("1", "true", "t", "yes", "y", "on")
        return true
    elseif low in ("0", "false", "f", "no", "n", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean in ENV[$key] = '$raw'. Use true/false, 1/0, yes/no."))
end

function parse_env_int(key::AbstractString, default::Int)
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return default
    try
        return parse(Int, raw)
    catch
        throw(ArgumentError("Invalid integer in ENV[$key] = '$raw'."))
    end
end

function parse_env_float(key::AbstractString, default::Real)
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return Float64(default)
    try
        return parse(Float64, raw)
    catch
        throw(ArgumentError("Invalid floating-point value in ENV[$key] = '$raw'."))
    end
end

function parse_env_optional_float(key::AbstractString, default::Union{Nothing,Float64}=nothing)
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return default
    low = lowercase(raw)
    if low in ("none", "nothing", "auto", "default")
        return nothing
    end
    try
        return parse(Float64, raw)
    catch
        throw(ArgumentError("Invalid optional floating-point value in ENV[$key] = '$raw'."))
    end
end

function parse_env_int_list(key::AbstractString, default::Vector{Int})
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return copy(default)
    vals = Int[]
    for piece in split(raw, ',')
        s = strip(piece)
        isempty(s) && continue
        try
            push!(vals, parse(Int, s))
        catch
            throw(ArgumentError("Invalid integer list in ENV[$key] = '$raw' (bad token '$s')."))
        end
    end
    isempty(vals) && return copy(default)
    return vals
end

N = parse_env_int("DMC_N", 3)                          # number of spinless fermions
a = parse_env_float("DMC_A", 1.0)                      # lattice period
L = parse_env_float("DMC_L", N * a)                    # ring length (common case: one particle per unit cell)
V0 = parse_env_float("DMC_V0", 1.0)                    # lattice amplitude
lambda_override = parse_env_optional_float("DMC_LAMBDA_OVERRIDE", nothing)   # override lambda = -V0/2

D = 0.5
node_tol = parse_env_float("DMC_NODE_TOL", 1e-7)      # near-node threshold for signpsi == 0
trig_eps = parse_env_float("DMC_TRIG_EPS", 1e-10)     # denominator guard for cot/csc^2 evaluations

targetN = parse_env_int("DMC_TARGET_N", 300)
dt = parse_env_float("DMC_DT", 0.003)
nsteps = parse_env_int("DMC_NSTEPS", 6000)
nequil = parse_env_int("DMC_NEQUIL", 400)
ET0 = parse_env_float("DMC_ET0", -0.2)
branch_cutoff = parse_env_int("DMC_BRANCH_CUTOFF", 2)
nblocks = haskey(ENV, "DMC_NBLOCKS") ? parse_env_int("DMC_NBLOCKS", 50) : parse_env_int("DMC_SAVE_EVERY", 50)

seed_init = parse_env_int("DMC_SEED_INIT", 1234)
seed_guided_fixed = parse_env_int("DMC_SEED_GUIDED_FIXED", 52)
seed_guided_no_node = parse_env_int("DMC_SEED_GUIDED_NO_NODE", 77)
seed_unguided = parse_env_int("DMC_SEED_UNGUIDED", 99)

snapshot_steps_user = parse_env_int_list("DMC_SNAPSHOT_STEPS", Int[0, nsteps ÷ 2, nsteps])

# Main mode is fixed-node guided. Debug comparators are optional.
run_fixed_node_guided = parse_env_bool("DMC_RUN_FIXED_GUIDED", true)
run_guided_no_node = parse_env_bool("DMC_RUN_GUIDED_NO_NODE", false)
run_unguided = parse_env_bool("DMC_RUN_UNGUIDED", false)

# Progress reporting
show_progress = parse_env_bool("DMC_SHOW_PROGRESS", true)
progress_every_steps = parse_env_int("DMC_PROGRESS_EVERY_STEPS", max(1, nsteps ÷ 20))

# Cloud/batch toggles
make_plots = parse_env_bool("DMC_MAKE_PLOTS", true)
output_root = abspath(get(ENV, "DMC_OUTPUT_DIR", joinpath(@__DIR__, "..", "outputs")))
summary_out = strip(get(ENV, "DMC_SUMMARY_OUT", ""))

# Plot/analysis settings
energy_clip_quantiles = (parse_env_float("DMC_ENERGY_QLO", 0.01), parse_env_float("DMC_ENERGY_QHI", 0.99))
energy_smoothing_window = parse_env_int("DMC_ENERGY_SMOOTH_WINDOW", 31)
density_nbins = parse_env_int("DMC_DENSITY_NBINS", 140)
density_smoothing_window = parse_env_int("DMC_DENSITY_SMOOTH_WINDOW", 11)
fd_grid_points = parse_env_int("DMC_FD_GRID_POINTS", 500)   # only used for N == 1


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
function sanitize_tag_component(x::Real)
    s = @sprintf("%.6f", Float64(x))
    s = replace(s, "-" => "m", "." => "p")
    s = replace(s, r"0+$" => "")
    s = replace(s, r"p$" => "")
    return s
end

function make_run_tag(N::Int, a::Real, L::Real, V0::Real)
    return "N$(N)_a$(sanitize_tag_component(a))_L$(sanitize_tag_component(L))_V0$(sanitize_tag_component(V0))"
end

function validate_inputs(N::Int, a::Real, L::Real)
    N >= 1 || throw(ArgumentError("N must be >= 1, got N=$N"))
    a > 0 || throw(ArgumentError("a must be > 0, got a=$a"))
    L > 0 || throw(ArgumentError("L must be > 0, got L=$L"))
end

function normalized_snapshot_steps(nsteps::Int, steps_user::AbstractVector{<:Integer})
    steps = sort(unique(vcat(Int[0, nsteps], Int.(steps_user))))
    for s in steps
        0 <= s <= nsteps || throw(ArgumentError("snapshot step $s is outside [0, $nsteps]"))
    end
    return steps
end

function assert_configuration_length(R::AbstractVector{<:Real}, N_expected::Int)
    length(R) == N_expected || throw(ArgumentError("Expected configuration length $N_expected, got $(length(R))"))
    return nothing
end

function wrapped_configuration(R::AbstractVector{<:Real}, bc::PeriodicBoundary1D, N_expected::Int)
    assert_configuration_length(R, N_expected)
    return wrap_position(bc, R)
end

function centered_moving_average_ignore_nan(values::AbstractVector{<:Real}, window::Int)
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

function quantile_clip_limits(series_list::AbstractVector{<:AbstractVector{<:Real}};
    qlo::Float64=0.01,
    qhi::Float64=0.99,
    pad_frac::Float64=0.08)
    0.0 <= qlo < qhi <= 1.0 || throw(ArgumentError("Need 0 <= qlo < qhi <= 1"))
    vals = Float64[]
    for y in series_list
        for v in y
            vf = Float64(v)
            isfinite(vf) && push!(vals, vf)
        end
    end
    isempty(vals) && return (-1.0, 1.0)
    lo, hi = quantile(vals, (qlo, qhi))
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
    for i in eachindex(values)
        v = Float64(values[i])
        out[i] = (isfinite(v) && lo <= v <= hi) ? v : NaN
    end
    return out
end

function format_energy_with_uncertainty_safe(E::Real, sem::Real)
    E_str = isfinite(E) ? @sprintf("%.10f", Float64(E)) : "NaN"
    sem_str = (isfinite(sem) && sem >= 0) ? @sprintf("%.3e", Float64(sem)) : "NaN"
    return E_str, sem_str
end

function post_equilibration_stats(sim::DMCSim, nequil::Int)
    start_idx = min(nequil + 1, length(sim.ET_history))
    post_ET = sim.ET_history[start_idx:end]
    post_var = sim.energy_variance_history[start_idx:end]

    Ebar = isempty(post_ET) ? NaN : mean(post_ET)
    sem = (length(post_ET) > 1) ? std(post_ET) / sqrt(length(post_ET)) : NaN
    varbar = isempty(post_var) ? NaN : mean(post_var)

    pop = sim.population_history
    pop_min = isempty(pop) ? 0 : minimum(pop)
    pop_max = isempty(pop) ? 0 : maximum(pop)
    pop_mean = isempty(pop) ? NaN : mean(pop)
    pop_final = isempty(pop) ? 0 : pop[end]

    return (
        Ebar=Ebar,
        sem=sem,
        varbar=varbar,
        start_idx=start_idx,
        pop_min=pop_min,
        pop_max=pop_max,
        pop_mean=pop_mean,
        pop_final=pop_final
    )
end

function pooled_onebody_coordinates(snap::AbstractVector{<:AbstractVector{<:Real}},
    N_expected::Int,
    bc::PeriodicBoundary1D)
    xs = Float64[]
    sizehint!(xs, length(snap) * N_expected)
    for R in snap
        assert_configuration_length(R, N_expected)
        for i in 1:N_expected
            push!(xs, wrap_coordinate(bc, R[i]))
        end
    end
    return xs
end

function density_histogram(xs::AbstractVector{<:Real}; nbins::Int=120, xmin::Real=0.0, xmax::Real=1.0)
    nbins >= 2 || throw(ArgumentError("nbins must be >= 2"))
    xmax > xmin || throw(ArgumentError("xmax must be greater than xmin"))
    isempty(xs) && throw(ArgumentError("Cannot build density histogram from an empty sample set"))

    width = (xmax - xmin) / nbins
    counts = zeros(Float64, nbins)
    for x in xs
        idx = floor(Int, (Float64(x) - xmin) / width) + 1
        idx = clamp(idx, 1, nbins)
        counts[idx] += 1.0
    end
    dens = counts ./ (length(xs) * width)
    centers = Float64[xmin + (i - 0.5) * width for i in 1:nbins]
    return centers, dens
end

function pooled_onebody_density_histogram(snap::AbstractVector{<:AbstractVector{<:Real}},
    N_expected::Int,
    bc::PeriodicBoundary1D;
    nbins::Int=120,
    xmin::Real=0.0,
    xmax::Real=cell_length(bc))
    xs = pooled_onebody_coordinates(snap, N_expected, bc)
    return density_histogram(xs; nbins=nbins, xmin=xmin, xmax=xmax)
end

function smoothed_density_curve(snap::AbstractVector{<:AbstractVector{<:Real}},
    N_expected::Int,
    bc::PeriodicBoundary1D;
    nbins::Int=120,
    smoothing_window::Int=9,
    xmin::Real=0.0,
    xmax::Real=cell_length(bc))
    centers, dens = pooled_onebody_density_histogram(
        snap,
        N_expected,
        bc;
        nbins=nbins,
        xmin=xmin,
        xmax=xmax
    )
    smooth = centered_moving_average_ignore_nan(dens, smoothing_window)
    return centers, smooth
end

function trapz(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(ArgumentError("trapz requires equal-length x and y"))
    n = length(x)
    n >= 2 || throw(ArgumentError("trapz requires at least two points"))
    s = 0.0
    for i in 1:(n - 1)
        dx = Float64(x[i + 1] - x[i])
        s += 0.5 * dx * (Float64(y[i]) + Float64(y[i + 1]))
    end
    return s
end

function normalized_trial_onebody_density(xgrid::AbstractVector{<:Real}, a::Real, lambda_::Real)
    k_lat = 2pi / a
    rho = Float64[exp(2 * lambda_ * cos(k_lat * x)) for x in xgrid]
    norm = trapz(xgrid, rho)
    norm > 0 || throw(ArgumentError("Trial density normalization failed (non-positive norm)."))
    rho ./= norm
    return rho
end

function finite_difference_reference_single_particle(L::Real, a::Real, V0::Real, ngrid::Int)
    ngrid >= 8 || throw(ArgumentError("ngrid must be >= 8"))
    dx = L / ngrid
    x = collect(0:(ngrid - 1)) .* dx
    k_lat = 2pi / a

    Hfd = zeros(Float64, ngrid, ngrid)
    kin_diag = 1.0 / dx^2
    kin_off = -0.5 / dx^2
    for i in 1:ngrid
        Hfd[i, i] = kin_diag + V0 * cos(k_lat * x[i])
    end
    for i in 1:(ngrid - 1)
        Hfd[i, i + 1] = kin_off
        Hfd[i + 1, i] = kin_off
    end
    Hfd[1, ngrid] = kin_off
    Hfd[ngrid, 1] = kin_off

    evals, evecs = eigen(Symmetric(Hfd))
    E0 = evals[1]
    psi0 = evecs[:, 1]
    rho0 = abs2.(psi0)
    rho0 ./= (sum(rho0) * dx)
    return x, E0, rho0
end

function build_base_positions(N::Int, targetN::Int, L::Real, seed::Int, bc::PeriodicBoundary1D)
    rng = MersenneTwister(seed)
    base = Vector{Vector{Float64}}(undef, targetN)
    for w in 1:targetN
        R = [rand(rng) * L for _ in 1:N]
        base[w] = wrap_position(bc, R)
    end
    return base
end

function build_node_policy(use_fixed_node::Bool)
    return use_fixed_node ? FixedNode() : NoNode()
end


# --------------------------------------------------------------------------
# Model setup and validation
# --------------------------------------------------------------------------
validate_inputs(N, a, L)

ncell_float = L / a
ncell_int = round(Int, ncell_float)
if abs(ncell_float - ncell_int) > 1e-9
    @warn "L/a is not close to an integer. External potential may be incompatible with strict ring periodicity." L a ratio=ncell_float
end

lambda_ = lambda_override === nothing ? (-V0 / 2) : Float64(lambda_override)
snapshot_steps = normalized_snapshot_steps(nsteps, snapshot_steps_user)

bc = PeriodicBoundary1D(0.0, L)
k_lat = 2pi / a
alpha_pair = pi / L

println("============================================================")
println("Spinless Fermion Ring DMC (1D, periodic)")
println("------------------------------------------------------------")
println("N               = ", N)
println("a               = ", @sprintf("%.6f", a))
println("L               = ", @sprintf("%.6f", L))
println("L/a             = ", @sprintf("%.6f", ncell_float), "  (~ ", ncell_int, " cells)")
println("V0              = ", @sprintf("%.6f", V0))
println("lambda          = ", @sprintf("%.6f", lambda_))
println("D               = ", @sprintf("%.6f", D))
println("node_tol        = ", @sprintf("%.2e", node_tol))
println("trig_eps        = ", @sprintf("%.2e", trig_eps))
println("targetN         = ", targetN)
println("dt              = ", @sprintf("%.6f", dt))
println("nsteps          = ", nsteps)
println("nequil          = ", nequil)
println("snapshot_steps  = ", snapshot_steps)
println("show_progress   = ", show_progress, " (every ", progress_every_steps, " steps)")
println("make_plots      = ", make_plots)
println("output_root     = ", output_root)
if !isempty(summary_out)
    println("summary_out     = ", summary_out)
end
println("run modes       = fixed-guided=", run_fixed_node_guided,
    ", guided-no-node=", run_guided_no_node,
    ", unguided=", run_unguided)
println("============================================================")


# --------------------------------------------------------------------------
# Hamiltonian and fermionic trial wavefunction
# --------------------------------------------------------------------------
function onebody_potential_sum(R::AbstractVector{<:Real})
    assert_configuration_length(R, N)
    Rw = wrapped_configuration(R, bc, N)
    V = 0.0
    @inbounds for i in 1:N
        V += V0 * cos(k_lat * Rw[i])
    end
    return V
end

function nodal_log_grad_lapl(R::AbstractVector{<:Real})
    Rw = wrapped_configuration(R, bc, N)
    grad = zeros(Float64, N)
    lapl = 0.0
    logabs = 0.0

    @inbounds for i in 1:(N - 1)
        for j in (i + 1):N
            dx_ij = displacement(bc, Rw[j], Rw[i])   # minimum-image signed (x_i - x_j)
            u = alpha_pair * dx_ij
            s = sin(u)
            abs_s = abs(s)

            # Guarded sine to avoid inf/NaN in log/cot/csc^2 near nodes.
            s_eff = abs_s < trig_eps ? (s >= 0 ? trig_eps : -trig_eps) : s
            logabs += log(abs(s_eff))

            cot_u = cos(u) / s_eff
            c = alpha_pair * cot_u
            grad[i] += c
            grad[j] -= c

            inv_s2 = inv(s_eff * s_eff)           # csc^2(u), guarded
            second = -(alpha_pair^2) * inv_s2     # d^2/dx_i^2 log|sin(u)|
            lapl += 2 * second                    # i and j each contribute equally
        end
    end
    return logabs, grad, lapl
end

function onebody_log_grad_lapl(R::AbstractVector{<:Real})
    Rw = wrapped_configuration(R, bc, N)
    grad = zeros(Float64, N)
    lapl = 0.0
    logg = 0.0
    @inbounds for i in 1:N
        xi = Rw[i]
        logg += lambda_ * cos(k_lat * xi)
        grad[i] = -lambda_ * k_lat * sin(k_lat * xi)
        lapl += -lambda_ * k_lat^2 * cos(k_lat * xi)
    end
    return logg, grad, lapl
end

function fermion_signpsi(R::AbstractVector{<:Real})
    Rw = wrapped_configuration(R, bc, N)
    sign_prod = 1.0
    @inbounds for i in 1:(N - 1)
        for j in (i + 1):N
            dx_ij = displacement(bc, Rw[j], Rw[i])  # x_i - x_j
            s = sin(alpha_pair * dx_ij)
            if abs(s) < node_tol
                return 0.0
            elseif s < 0
                sign_prod = -sign_prod
            end
        end
    end
    return sign_prod
end

trial_logpsi(R) = begin
    lnA, _, _ = nodal_log_grad_lapl(R)
    lnG, _, _ = onebody_log_grad_lapl(R)
    return lnA + lnG
end

trial_gradlogpsi(R) = begin
    _, gA, _ = nodal_log_grad_lapl(R)
    _, gG, _ = onebody_log_grad_lapl(R)
    return gA .+ gG
end

trial_lapllogpsi(R) = begin
    _, _, lA = nodal_log_grad_lapl(R)
    _, _, lG = onebody_log_grad_lapl(R)
    return lA + lG
end

H = Hamiltonian(N, D, onebody_potential_sum, bc)
trial = TrialWF(trial_logpsi, trial_gradlogpsi, trial_lapllogpsi, fermion_signpsi)
guiding = ImportanceGuiding(trial, H)


# --------------------------------------------------------------------------
# Optional one-particle finite-difference reference (N == 1 only)
# --------------------------------------------------------------------------
fd_ref = nothing
if N == 1
    x_ref, E0_ref, rho0_ref = finite_difference_reference_single_particle(L, a, V0, fd_grid_points)
    fd_ref = (x=x_ref, E0=E0_ref, rho=rho0_ref)
else
    @warn "Skipping exact reference diagonalization for N > 1 (no many-body exact solver in this script)." N
end


# --------------------------------------------------------------------------
# Run one or more DMC modes
# --------------------------------------------------------------------------
run_specs = NamedTuple[]
if run_fixed_node_guided
    push!(run_specs, (name="guided_fixed_node", use_guiding=true, use_fixed_node=true, seed=seed_guided_fixed, color=:teal))
end
if run_guided_no_node
    push!(run_specs, (name="guided_no_node", use_guiding=true, use_fixed_node=false, seed=seed_guided_no_node, color=:darkorange))
end
if run_unguided
    push!(run_specs, (name="unguided_no_node", use_guiding=false, use_fixed_node=false, seed=seed_unguided, color=:firebrick))
end
isempty(run_specs) && throw(ArgumentError("No run mode selected. Enable at least one of run_fixed_node_guided/run_guided_no_node/run_unguided."))

params = DMCParams(; dt=dt, nsteps=nsteps, nequil=nequil, targetN=targetN, ET0=ET0, branch_cap=branch_cutoff, nblocks=nblocks)
base_positions = build_base_positions(N, targetN, L, seed_init, bc)

results = NamedTuple[]
for spec in run_specs
    rng_sim = MersenneTwister(spec.seed)
    guiding_mode = spec.use_guiding ? guiding : NoGuiding()
    nodepolicy = build_node_policy(spec.use_fixed_node)

    sim = run_dmc(
        H,
        params,
        base_positions;
        rng=rng_sim,
        guiding=guiding_mode,
        nodepolicy=nodepolicy,
        snapshot_steps=snapshot_steps,
        show_progress=show_progress,
        progress_every=progress_every_steps,
        progress_label=spec.name
    )
    push!(results, (spec=spec, sim=sim))
end

primary_idx = findfirst(r -> r.spec.name == "guided_fixed_node", results)
primary = results[primary_idx === nothing ? 1 : primary_idx]


# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
println("\nRun summaries:")
summary_rows = NamedTuple[]
for r in results
    stats = post_equilibration_stats(r.sim, nequil)
    push!(summary_rows, (
        mode=String(r.spec.name),
        use_guiding=Bool(r.spec.use_guiding),
        use_fixed_node=Bool(r.spec.use_fixed_node),
        seed=Int(r.spec.seed),
        Ebar=Float64(stats.Ebar),
        sem=Float64(stats.sem),
        varbar=Float64(stats.varbar),
        pop_min=Int(stats.pop_min),
        pop_mean=Float64(stats.pop_mean),
        pop_max=Int(stats.pop_max),
        pop_final=Int(stats.pop_final)
    ))
    E_str, dE_str = format_energy_with_uncertainty_safe(stats.Ebar, stats.sem)
    var_str = isfinite(stats.varbar) ? @sprintf("%.8f", stats.varbar) : "NaN"

    println("------------------------------------------------------------")
    println("Mode                  : ", r.spec.name)
    println("Guiding               : ", r.spec.use_guiding)
    println("Fixed node            : ", r.spec.use_fixed_node)
    println("Seed                  : ", r.spec.seed)
    println("Post-eq ET            : ", E_str, " +/- ", dE_str)
    println("Post-eq variance      : ", var_str)
    println("Population (min/mean/max/final): ",
        stats.pop_min, " / ",
        @sprintf("%.2f", stats.pop_mean), " / ",
        stats.pop_max, " / ",
        stats.pop_final)
    if fd_ref !== nothing
        println("FD reference E0       : ", @sprintf("%.10f", fd_ref.E0))
    end
end
println("------------------------------------------------------------")

if !isempty(summary_out)
    mkpath(dirname(summary_out))
    open(summary_out, "w") do io
        println(io, "mode,use_guiding,use_fixed_node,seed,Ebar,sem,varbar,pop_min,pop_mean,pop_max,pop_final,nsteps,nequil,targetN,dt,N,a,L,V0,lambda")
        for row in summary_rows
            @printf(
                io,
                "%s,%s,%s,%d,%.16e,%.16e,%.16e,%d,%.8f,%d,%d,%d,%d,%d,%.16e,%d,%.16e,%.16e,%.16e,%.16e\n",
                row.mode,
                string(row.use_guiding),
                string(row.use_fixed_node),
                row.seed,
                row.Ebar,
                row.sem,
                row.varbar,
                row.pop_min,
                row.pop_mean,
                row.pop_max,
                row.pop_final,
                nsteps,
                nequil,
                targetN,
                dt,
                N,
                a,
                L,
                V0,
                lambda_
            )
        end
    end
    println("Wrote summary CSV: ", summary_out)
end


# --------------------------------------------------------------------------
# Plotting and outputs
# --------------------------------------------------------------------------
if make_plots
    figdir = abspath(joinpath(output_root, "figures", "fermion_ring"))
    mkpath(figdir)
    tag = make_run_tag(N, a, L, V0)
    t = (0:nsteps) .* dt

    Logging.with_logger(Logging.NullLogger()) do
        redirect_stderr(devnull) do
            # 1) Reference-energy history ET(tau), clipped + smoothed
            et_post_series = [r.sim.ET_history[min(nequil + 1, length(r.sim.ET_history)):end] for r in results]
            et_lo, et_hi = quantile_clip_limits(
                et_post_series;
                qlo=energy_clip_quantiles[1],
                qhi=energy_clip_quantiles[2],
                pad_frac=0.08
            )
            p_et = plot(
                xlabel="imaginary time tau",
                ylabel="E_T(tau)",
                title="Fermion ring: reference energy history",
                ylims=(et_lo, et_hi),
                legend=:topright
            )
            for r in results
                masked = mask_outside_range(r.sim.ET_history, et_lo, et_hi)
                smooth = centered_moving_average_ignore_nan(masked, energy_smoothing_window)
                plot!(p_et, t, smooth; linewidth=2.1, color=r.spec.color, label="$(r.spec.name) (smoothed)")
            end
            if fd_ref !== nothing
                hline!(p_et, [fd_ref.E0]; color=:black, linestyle=:dash, linewidth=1.8, label="FD E0")
            end
            savefig(p_et, joinpath(figdir, "$(tag)_ET_history.png"))

            # 2) Mean-estimator history
            mean_post_series = [r.sim.energy_mean_history[min(nequil + 1, length(r.sim.energy_mean_history)):end] for r in results]
            m_lo, m_hi = quantile_clip_limits(
                mean_post_series;
                qlo=energy_clip_quantiles[1],
                qhi=energy_clip_quantiles[2],
                pad_frac=0.08
            )
            p_mean = plot(
                xlabel="imaginary time tau",
                ylabel="mean estimator",
                title="Fermion ring: mean energy estimator history",
                ylims=(m_lo, m_hi),
                legend=:topright
            )
            for r in results
                masked = mask_outside_range(r.sim.energy_mean_history, m_lo, m_hi)
                smooth = centered_moving_average_ignore_nan(masked, energy_smoothing_window)
                plot!(p_mean, t, smooth; linewidth=2.1, color=r.spec.color, label="$(r.spec.name) (smoothed)")
            end
            if fd_ref !== nothing
                hline!(p_mean, [fd_ref.E0]; color=:black, linestyle=:dash, linewidth=1.8, label="FD E0")
            end
            savefig(p_mean, joinpath(figdir, "$(tag)_mean_estimator_history.png"))

            # 3) Variance history
            var_post_series = [r.sim.energy_variance_history[min(nequil + 1, length(r.sim.energy_variance_history)):end] for r in results]
            v_lo, v_hi = quantile_clip_limits(var_post_series; qlo=0.005, qhi=0.995, pad_frac=0.10)
            p_var = plot(
                xlabel="imaginary time tau",
                ylabel="Var(E estimator)",
                title="Fermion ring: variance history",
                ylims=(max(0.0, v_lo), v_hi),
                legend=:topright
            )
            for r in results
                masked = mask_outside_range(r.sim.energy_variance_history, v_lo, v_hi)
                smooth = centered_moving_average_ignore_nan(masked, energy_smoothing_window)
                plot!(p_var, t, smooth; linewidth=2.1, color=r.spec.color, label="$(r.spec.name) (smoothed)")
            end
            savefig(p_var, joinpath(figdir, "$(tag)_variance_history.png"))

            # 4) External potential over one ring period
            xplot = range(0.0, L; length=1200)
            Vext_plot = Float64[V0 * cos(k_lat * x) for x in xplot]
            p_vext = plot(
                xplot,
                Vext_plot;
                xlabel="x",
                ylabel="V_ext(x)",
                title="External periodic potential on ring",
                linewidth=2.2,
                color=:navy,
                label="V0 cos(2pi x / a)"
            )
            savefig(p_vext, joinpath(figdir, "$(tag)_external_potential.png"))

            # 5) One-body drift field illustration: F_1body(x) = 2 * d/dx log G
            F1_plot = Float64[2 * (-lambda_ * k_lat * sin(k_lat * x)) for x in xplot]
            p_drift = plot(
                xplot,
                F1_plot;
                xlabel="x",
                ylabel="F_1body(x)",
                title="One-body trial drift field",
                linewidth=2.2,
                color=:purple4,
                label="2 * (-lambda * 2pi/a * sin(2pi x/a))"
            )
            hline!(p_drift, [0.0]; color=:black, linestyle=:dot, linewidth=1.3, label=false)
            savefig(p_drift, joinpath(figdir, "$(tag)_onebody_drift_field.png"))

            # 6) Final pooled one-body density diagnostics (primary run)
            final_snap = primary.sim.walker_positions_history[end]
            xc_final, rho_final = pooled_onebody_density_histogram(
                final_snap,
                N,
                bc;
                nbins=density_nbins,
                xmin=0.0,
                xmax=L
            )
            rho_final_smooth = centered_moving_average_ignore_nan(rho_final, density_smoothing_window)

            xref = range(0.0, L; length=1200)
            rho_trial_1b = normalized_trial_onebody_density(xref, a, lambda_)

            p_final = plot(
                xc_final,
                rho_final_smooth;
                xlabel="x",
                ylabel="one-body density",
                title="Final pooled one-body density (primary run)",
                linewidth=2.3,
                color=:teal,
                label="DMC pooled density (smoothed)"
            )
            plot!(p_final, xref, rho_trial_1b; linewidth=2.0, color=:gray30, linestyle=:dash, label="trial one-body shape")
            if fd_ref !== nothing
                plot!(p_final, fd_ref.x, fd_ref.rho; linewidth=2.0, color=:black, linestyle=:dot, label="FD |psi0|^2")
            end
            savefig(p_final, joinpath(figdir, "$(tag)_final_onebody_density.png"))

            # 7) Smoothed snapshot density evolution (primary run)
            nsnaps = min(length(snapshot_steps), length(primary.sim.walker_positions_history))
            snap_colors = [:royalblue, :darkorange, :forestgreen, :firebrick, :purple, :goldenrod, :teal]
            p_snaps = plot(
                xlabel="x",
                ylabel="one-body density",
                xlims=(0.0, L),
                title="Pooled one-body density snapshots (primary run)",
                legend=:topright
            )
            for i in 1:nsnaps
                s = snapshot_steps[i]
                snap = primary.sim.walker_positions_history[i]
                xc, yc = smoothed_density_curve(
                    snap,
                    N,
                    bc;
                    nbins=density_nbins,
                    smoothing_window=density_smoothing_window,
                    xmin=0.0,
                    xmax=L
                )
                tau_label = @sprintf("tau=%.4f", s * dt)
                col = snap_colors[mod1(i, length(snap_colors))]
                plot!(p_snaps, xc, yc; linewidth=2.0, color=col, label=tau_label)
            end
            plot!(p_snaps, xref, rho_trial_1b; linewidth=2.0, color=:gray30, linestyle=:dash, label="trial one-body shape")
            if fd_ref !== nothing
                plot!(p_snaps, fd_ref.x, fd_ref.rho; linewidth=2.0, color=:black, linestyle=:dot, label="FD |psi0|^2")
            end
            savefig(p_snaps, joinpath(figdir, "$(tag)_snapshot_density_evolution.png"))
        end
    end

    println("Saved figures to: ", figdir)
else
    println("Plot generation disabled (DMC_MAKE_PLOTS=false).")
end
