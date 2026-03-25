# Helpers for periodic ion-ring GFMC benchmarks built on the generic Hamiltonian
# path. These utilities keep the notebooks compact while exposing the exact
# one-body orbital construction, determinant-based trial states, and general
# bosonic Jastrow/Tonks scaffolds.

using Random
using LinearAlgebra

struct PeriodicIonRingModel
    M::Int
    a::Float64
    L::Float64
    D::Float64
    ion_strength::Float64
    ion_softening::Float64
    bc::PeriodicBoundary1D{Float64}
    ion_positions::Vector{Float64}
    basis_kinds::Vector{Symbol}
    basis_indices::Vector{Int}
    coeffs::Matrix{Float64}
    energies::Vector{Float64}
    xquad::Vector{Float64}
    potential_quad::Vector{Float64}
end

struct PeriodicIonRingBosonProblem
    N::Int
    ring::PeriodicIonRingModel
    bb_strength::Float64
    bb_softening::Float64
    hard_core_radius::Float64
    hard_core_barrier::Float64
    pair_metric::Symbol
end

function _build_real_fourier_basis(kmax::Integer)
    kmax_int = Int(kmax)
    kmax_int >= 0 || throw(ArgumentError("kmax must be >= 0, got $kmax"))

    kinds = Symbol[:const]
    indices = Int[0]
    for k in 1:kmax_int
        push!(kinds, :cos)
        push!(indices, k)
        push!(kinds, :sin)
        push!(indices, k)
    end
    return kinds, indices
end

function _basis_value(kind::Symbol, k::Int, x::Real, L::Real)
    xw = Float64(mod(Float64(x), Float64(L)))
    if kind === :const
        return inv(sqrt(Float64(L)))
    end

    G = (2pi * k) / Float64(L)
    norm = sqrt(2.0 / Float64(L))
    if kind === :cos
        return norm * cos(G * xw)
    elseif kind === :sin
        return norm * sin(G * xw)
    end
    throw(ArgumentError("Unknown basis kind: $kind"))
end

function _basis_grad(kind::Symbol, k::Int, x::Real, L::Real)
    xw = Float64(mod(Float64(x), Float64(L)))
    if kind === :const
        return 0.0
    end

    G = (2pi * k) / Float64(L)
    norm = sqrt(2.0 / Float64(L))
    if kind === :cos
        return -norm * G * sin(G * xw)
    elseif kind === :sin
        return norm * G * cos(G * xw)
    end
    throw(ArgumentError("Unknown basis kind: $kind"))
end

function _basis_lapl(kind::Symbol, k::Int, x::Real, L::Real)
    if kind === :const
        return 0.0
    end
    G = (2pi * k) / Float64(L)
    return -(G * G) * _basis_value(kind, k, x, L)
end

function _basis_kinetic_eigenvalue(kind::Symbol, k::Int, D::Real, L::Real)
    if kind === :const
        return 0.0
    end
    G = (2pi * k) / Float64(L)
    return Float64(D) * G * G
end

function periodic_ion_lattice_potential(
    x::Real,
    ion_positions::AbstractVector{<:Real},
    bc::PeriodicBoundary1D;
    ion_strength::Real,
    ion_softening::Real,
)
    strength = Float64(ion_strength)
    soft = Float64(ion_softening)
    soft > 0 || throw(ArgumentError("ion_softening must be > 0, got $ion_softening"))

    xw = wrap_coordinate(bc, Float64(x))
    total = 0.0
    @inbounds for X in ion_positions
        r = distance_1d(bc, xw, X)
        total -= strength / sqrt(r * r + soft * soft)
    end
    return total
end

function _resolved_ion_positions(
    M::Integer,
    a::Real;
    ion_positions::Union{Nothing,AbstractVector{<:Real}}=nothing,
    L::Union{Nothing,Real}=nothing,
)
    M_int = Int(M)
    M_int >= 1 || throw(ArgumentError("M must be >= 1, got $M"))
    a_f = Float64(a)
    a_f > 0 || throw(ArgumentError("a must be > 0, got $a"))

    L_f = L === nothing ? (M_int * a_f) : Float64(L)
    L_f > 0 || throw(ArgumentError("L must be > 0, got $L_f"))
    bc = PeriodicBoundary1D(0.0, L_f)

    ions = if ion_positions === nothing
        Float64[m * (L_f / M_int) for m in 0:(M_int - 1)]
    else
        length(ion_positions) == M_int || throw(ArgumentError(
            "ion_positions length $(length(ion_positions)) must match M=$M_int",
        ))
        Float64[wrap_coordinate(bc, x) for x in ion_positions]
    end
    return L_f, bc, ions
end

function build_periodic_ion_ring_model(
    M::Integer,
    a::Real;
    D::Real=0.5,
    ion_strength::Real=1.0,
    ion_softening::Real=0.35 * Float64(a),
    ion_positions::Union{Nothing,AbstractVector{<:Real}}=nothing,
    L::Union{Nothing,Real}=nothing,
    kmax::Integer=6,
    quad_points::Integer=1536,
)
    M_int = Int(M)
    M_int >= 1 || throw(ArgumentError("M must be >= 1, got $M"))
    a_f = Float64(a)
    a_f > 0 || throw(ArgumentError("a must be > 0, got $a"))

    L_f, bc, ions = _resolved_ion_positions(
        M_int,
        a_f;
        ion_positions=ion_positions,
        L=L,
    )

    quad_int = Int(quad_points)
    quad_int >= 64 || throw(ArgumentError("quad_points must be >= 64, got $quad_points"))
    dx = L_f / quad_int
    xquad = Float64[i * dx for i in 0:(quad_int - 1)]
    potential_quad = Float64[
        periodic_ion_lattice_potential(
            x,
            ions,
            bc;
            ion_strength=ion_strength,
            ion_softening=ion_softening,
        ) for x in xquad
    ]

    basis_kinds, basis_indices = _build_real_fourier_basis(kmax)
    nbasis = length(basis_kinds)
    B = Matrix{Float64}(undef, quad_int, nbasis)
    Hkin = zeros(Float64, nbasis, nbasis)

    @inbounds for j in 1:nbasis
        kind = basis_kinds[j]
        k = basis_indices[j]
        Hkin[j, j] = _basis_kinetic_eigenvalue(kind, k, D, L_f)
        for i in 1:quad_int
            B[i, j] = _basis_value(kind, k, xquad[i], L_f)
        end
    end

    weighted_basis = similar(B)
    @inbounds for j in 1:nbasis
        for i in 1:quad_int
            weighted_basis[i, j] = potential_quad[i] * B[i, j]
        end
    end
    Hpot = dx .* transpose(B) * weighted_basis

    evals, evecs = eigen(Symmetric(Hkin + Hpot))
    coeffs = Matrix{Float64}(evecs)

    return PeriodicIonRingModel(
        M_int,
        a_f,
        L_f,
        Float64(D),
        Float64(ion_strength),
        Float64(ion_softening),
        bc,
        ions,
        basis_kinds,
        basis_indices,
        coeffs,
        Float64.(evals),
        xquad,
        potential_quad,
    )
end

function build_periodic_ion_ring_boson_problem(
    N::Integer,
    M::Integer,
    a::Real;
    D::Real=0.5,
    ion_strength::Real=1.0,
    ion_softening::Real=0.35 * Float64(a),
    ion_positions::Union{Nothing,AbstractVector{<:Real}}=nothing,
    L::Union{Nothing,Real}=nothing,
    bb_strength::Real=0.0,
    bb_softening::Real=0.25 * Float64(a),
    hard_core_radius::Real=0.0,
    hard_core_barrier::Real=0.0,
    pair_metric::Symbol=:chord,
    kmax::Integer=6,
    quad_points::Integer=1536,
)
    N_int = Int(N)
    N_int >= 1 || throw(ArgumentError("N must be >= 1, got $N"))
    pair_metric in (:minimum_image, :chord) || throw(ArgumentError(
        "pair_metric must be :minimum_image or :chord, got $pair_metric",
    ))

    ring = build_periodic_ion_ring_model(
        M,
        a;
        D=D,
        ion_strength=ion_strength,
        ion_softening=ion_softening,
        ion_positions=ion_positions,
        L=L,
        kmax=kmax,
        quad_points=quad_points,
    )

    return PeriodicIonRingBosonProblem(
        N_int,
        ring,
        Float64(bb_strength),
        Float64(bb_softening),
        Float64(hard_core_radius),
        Float64(hard_core_barrier),
        pair_metric,
    )
end

function onebody_potential(model::PeriodicIonRingModel, x::Real)
    return periodic_ion_lattice_potential(
        x,
        model.ion_positions,
        model.bc;
        ion_strength=model.ion_strength,
        ion_softening=model.ion_softening,
    )
end

function manybody_potential(model::PeriodicIonRingModel, R::AbstractVector{<:Real})
    total = 0.0
    @inbounds for x in R
        total += onebody_potential(model, x)
    end
    return total
end

function periodic_ion_hamiltonian(model::PeriodicIonRingModel, N::Integer)
    N_int = Int(N)
    N_int >= 1 || throw(ArgumentError("N must be >= 1, got $N"))
    return Hamiltonian(N_int, model.D, R -> manybody_potential(model, R), model.bc)
end

function pair_distance(
    bc::PeriodicBoundary1D,
    x1::Real,
    x2::Real;
    metric::Symbol=:minimum_image,
)
    if metric === :minimum_image
        return distance_1d(bc, x1, x2)
    elseif metric === :chord
        dx = displacement(bc, x2, x1)
        return (cell_length(bc) / pi) * abs(sin((pi / cell_length(bc)) * dx))
    end
    throw(ArgumentError("Unknown pair metric: $metric"))
end

function boson_pair_potential(problem::PeriodicIonRingBosonProblem, x1::Real, x2::Real)
    bc = problem.ring.bc
    r_metric = pair_distance(bc, x1, x2; metric=problem.pair_metric)
    repulsion = problem.bb_strength == 0 ? 0.0 :
        problem.bb_strength / sqrt(r_metric^2 + problem.bb_softening^2)

    if problem.hard_core_radius > 0
        r_min = distance_1d(bc, x1, x2)
        if r_min < problem.hard_core_radius
            repulsion += problem.hard_core_barrier
        end
    end
    return repulsion
end

pair_distance(problem::PeriodicIonRingBosonProblem, x1::Real, x2::Real) =
    pair_distance(problem.ring.bc, x1, x2; metric=problem.pair_metric)

function bosonic_manybody_potential(problem::PeriodicIonRingBosonProblem, R::AbstractVector{<:Real})
    length(R) == problem.N || throw(ArgumentError("Expected $(problem.N) coordinates, got $(length(R))"))
    total = manybody_potential(problem.ring, R)
    @inbounds for i in 1:(problem.N - 1)
        for j in (i + 1):problem.N
            total += boson_pair_potential(problem, R[i], R[j])
        end
    end
    return total
end

function bosonic_hamiltonian(problem::PeriodicIonRingBosonProblem)
    return Hamiltonian(problem.N, problem.ring.D, R -> bosonic_manybody_potential(problem, R), problem.ring.bc)
end

noninteracting_boson_energy(model::PeriodicIonRingModel, N::Integer) = Int(N) * model.energies[1]
tg_scaffold_energy(model::PeriodicIonRingModel, N::Integer) = sum(@view model.energies[1:Int(N)])

function _pair_kernel_terms(
    bc::PeriodicBoundary1D,
    xi::Real,
    xj::Real;
    pair_power::Real=0.0,
    pair_eps::Real=1.0e-8,
    smooth_pair_strength::Real=0.0,
    smooth_pair_softening::Real=1.0,
)
    α = pi / cell_length(bc)
    Δ = displacement(bc, xj, xi)
    u = α * Δ
    s = sin(u)
    c = cos(u)

    logterm = 0.0
    grad_i = 0.0
    lapl_pair = 0.0

    γ = Float64(pair_power)
    if γ != 0.0
        abs_s = abs(s)
        s_eff = abs_s < pair_eps ? (s >= 0 ? pair_eps : -pair_eps) : s
        logterm += γ * log(abs(s_eff))
        grad_i += γ * α * c / s_eff
        lapl_pair += -2.0 * γ * α^2 / (s_eff * s_eff)
    end

    η = Float64(smooth_pair_strength)
    if η != 0.0
        σ = Float64(smooth_pair_softening)
        σ > 0 || throw(ArgumentError("smooth_pair_softening must be > 0"))
        q = σ^2 + (s / α)^2
        ρ = sqrt(q)
        logterm += -η / ρ
        grad_i += η * s * c / (α * ρ^3)
        lapl_pair += 2.0 * (
            η * (c^2 - s^2) / (ρ^3) -
            3.0 * η * s^2 * c^2 / (α^2 * ρ^5)
        )
    end

    return logterm, grad_i, lapl_pair
end

function pair_metric_radius(L::Real, r::Real; metric::Symbol=:minimum_image)
    if metric === :minimum_image
        return Float64(r)
    elseif metric === :chord
        return (Float64(L) / pi) * sin((pi / Float64(L)) * Float64(r))
    end
    throw(ArgumentError("Unknown pair metric: $metric"))
end

function pair_potential_curve(problem::PeriodicIonRingBosonProblem, rgrid::AbstractVector{<:Real})
    curve = Vector{Float64}(undef, length(rgrid))
    @inbounds for i in eachindex(rgrid)
        r = Float64(rgrid[i])
        r_metric = pair_metric_radius(problem.ring.L, r; metric=problem.pair_metric)
        repulsion = problem.bb_strength == 0 ? 0.0 :
            problem.bb_strength / sqrt(r_metric^2 + problem.bb_softening^2)
        if problem.hard_core_radius > 0 && r < problem.hard_core_radius
            repulsion += problem.hard_core_barrier
        end
        curve[i] = repulsion
    end
    return curve
end

function jastrow_pair_factor_curve(
    L::Real,
    rgrid::AbstractVector{<:Real};
    pair_power::Real=1.0,
    smooth_pair_strength::Real=0.0,
    smooth_pair_softening::Real=1.0,
)
    α = pi / Float64(L)
    out = Vector{Float64}(undef, length(rgrid))
    @inbounds for i in eachindex(rgrid)
        r = Float64(rgrid[i])
        s = sin(α * r)
        pair_factor = abs(s)^Float64(pair_power)
        chord = abs(s) / α
        smooth_factor = exp(-Float64(smooth_pair_strength) / sqrt(chord^2 + Float64(smooth_pair_softening)^2))
        out[i] = pair_factor * smooth_factor
    end
    return out
end

function smooth_pair_factor_curve(
    L::Real,
    rgrid::AbstractVector{<:Real};
    smooth_pair_strength::Real=0.0,
    smooth_pair_softening::Real=1.0,
)
    return jastrow_pair_factor_curve(
        L,
        rgrid;
        pair_power=0.0,
        smooth_pair_strength=smooth_pair_strength,
        smooth_pair_softening=smooth_pair_softening,
    )
end

function _basis_triple(model::PeriodicIonRingModel, x::Real)
    nbasis = length(model.basis_kinds)
    vals = Vector{Float64}(undef, nbasis)
    grads = Vector{Float64}(undef, nbasis)
    lapls = Vector{Float64}(undef, nbasis)

    @inbounds for j in 1:nbasis
        kind = model.basis_kinds[j]
        k = model.basis_indices[j]
        vals[j] = _basis_value(kind, k, x, model.L)
        grads[j] = _basis_grad(kind, k, x, model.L)
        lapls[j] = _basis_lapl(kind, k, x, model.L)
    end
    return vals, grads, lapls
end

function orbital_values(model::PeriodicIonRingModel, x::Real; norb::Integer=size(model.coeffs, 2))
    norb_int = Int(norb)
    1 <= norb_int <= size(model.coeffs, 2) || throw(ArgumentError("norb must be in 1:$(size(model.coeffs, 2)), got $norb"))

    basis_vals, basis_grads, basis_lapls = _basis_triple(model, x)
    vals = Vector{Float64}(undef, norb_int)
    grads = Vector{Float64}(undef, norb_int)
    lapls = Vector{Float64}(undef, norb_int)

    @inbounds for a in 1:norb_int
        coeff_col = view(model.coeffs, :, a)
        vals[a] = dot(coeff_col, basis_vals)
        grads[a] = dot(coeff_col, basis_grads)
        lapls[a] = dot(coeff_col, basis_lapls)
    end
    return vals, grads, lapls
end

function orbital_curve(model::PeriodicIonRingModel, orbital_idx::Integer, xgrid::AbstractVector{<:Real})
    idx = Int(orbital_idx)
    vals = Vector{Float64}(undef, length(xgrid))
    grads = Vector{Float64}(undef, length(xgrid))
    lapls = Vector{Float64}(undef, length(xgrid))
    @inbounds for i in eachindex(xgrid)
        v, g, l = orbital_values(model, xgrid[i]; norb=idx)
        vals[i] = v[idx]
        grads[i] = g[idx]
        lapls[i] = l[idx]
    end
    return vals, grads, lapls
end

function occupied_pooled_density(
    model::PeriodicIonRingModel,
    occupied::AbstractVector{<:Integer},
    xgrid::AbstractVector{<:Real};
    per_particle::Bool=true,
)
    isempty(occupied) && throw(ArgumentError("occupied orbital list cannot be empty"))
    max_orb = maximum(occupied)
    density = zeros(Float64, length(xgrid))
    @inbounds for i in eachindex(xgrid)
        vals, _, _ = orbital_values(model, xgrid[i]; norb=max_orb)
        accum = 0.0
        for orb in occupied
            accum += vals[Int(orb)]^2
        end
        density[i] = accum
    end
    if per_particle
        density ./= length(occupied)
    end
    return density
end

function determinant_trial_terms(
    model::PeriodicIonRingModel,
    R::AbstractVector{<:Real},
    N::Integer;
    node_tol::Real=1.0e-10,
)
    N_int = Int(N)
    length(R) == N_int || throw(ArgumentError("Expected $(N_int) coordinates, got $(length(R))"))
    1 <= N_int <= size(model.coeffs, 2) || throw(ArgumentError("Need 1 <= N <= $(size(model.coeffs, 2)), got $N"))
    tol = Float64(node_tol)
    tol > 0 || throw(ArgumentError("node_tol must be > 0, got $node_tol"))

    Phi = Matrix{Float64}(undef, N_int, N_int)
    dPhi = Matrix{Float64}(undef, N_int, N_int)
    ddPhi = Matrix{Float64}(undef, N_int, N_int)
    @inbounds for i in 1:N_int
        vals, grads, lapls = orbital_values(model, R[i]; norb=N_int)
        for a in 1:N_int
            Phi[i, a] = vals[a]
            dPhi[i, a] = grads[a]
            ddPhi[i, a] = lapls[a]
        end
    end

    detPhi = det(Phi)
    absdet = abs(detPhi)
    sign = absdet < tol ? 0.0 : (detPhi >= 0 ? 1.0 : -1.0)
    logabs = log(max(absdet, tol))

    Phi_reg = absdet < tol ? (Phi + tol * I) : Phi
    A = inv(Phi_reg)

    grad = Vector{Float64}(undef, N_int)
    lapllog = 0.0
    @inbounds for i in 1:N_int
        gi = dot(view(A, :, i), view(dPhi, i, :))
        grad[i] = gi
        lapllog += dot(view(A, :, i), view(ddPhi, i, :)) - gi * gi
    end
    return logabs, grad, lapllog, sign
end

function single_particle_trial_wavefunction(
    model::PeriodicIonRingModel;
    orbital_index::Integer=1,
    amp_tol::Real=1.0e-10,
)
    idx = Int(orbital_index)
    tol = Float64(amp_tol)
    tol > 0 || throw(ArgumentError("amp_tol must be > 0, got $amp_tol"))

    function _single_terms(R)
        length(R) == 1 || throw(ArgumentError("Expected a single-particle configuration, got $(length(R)) coordinates"))
        vals, grads, lapls = orbital_values(model, R[1]; norb=idx)
        psi = vals[idx]
        grad_psi = grads[idx]
        lapl_psi = lapls[idx]
        sign = abs(psi) < tol ? 0.0 : (psi >= 0 ? 1.0 : -1.0)
        psi_eff = abs(psi) < tol ? (psi >= 0 ? tol : -tol) : psi
        gradlog = grad_psi / psi_eff
        lapllog = lapl_psi / psi_eff - gradlog^2
        return log(abs(psi_eff)), Float64[gradlog], lapllog, sign
    end

    return TrialWF(
        R -> _single_terms(R)[1],
        R -> _single_terms(R)[2],
        R -> _single_terms(R)[3],
        R -> _single_terms(R)[4],
    )
end

function bosonic_orbital_jastrow_trial_wavefunction(
    problem::PeriodicIonRingBosonProblem;
    orbital_index::Integer=1,
    orbital_tol::Real=1.0e-10,
    pair_power::Real=1.0,
    pair_eps::Real=1.0e-8,
    smooth_pair_strength::Real=0.0,
    smooth_pair_softening::Real=max(problem.bb_softening, 0.25 * problem.ring.a),
)
    idx = Int(orbital_index)
    1 <= idx <= size(problem.ring.coeffs, 2) || throw(ArgumentError(
        "orbital_index must be in 1:$(size(problem.ring.coeffs, 2)), got $orbital_index",
    ))
    tol = Float64(orbital_tol)
    tol > 0 || throw(ArgumentError("orbital_tol must be > 0, got $orbital_tol"))

    function _trial_terms(R)
        length(R) == problem.N || throw(ArgumentError("Expected $(problem.N) coordinates, got $(length(R))"))

        logabs = 0.0
        grad = zeros(Float64, problem.N)
        lapllog = 0.0

        @inbounds for i in 1:problem.N
            vals, grads, lapls = orbital_values(problem.ring, R[i]; norb=idx)
            psi = vals[idx]
            grad_psi = grads[idx]
            lapl_psi = lapls[idx]
            psi_eff = abs(psi) < tol ? (psi >= 0 ? tol : -tol) : psi
            grad_i = grad_psi / psi_eff
            logabs += log(abs(psi_eff))
            grad[i] += grad_i
            lapllog += lapl_psi / psi_eff - grad_i^2
        end

        @inbounds for i in 1:(problem.N - 1)
            for j in (i + 1):problem.N
                logterm, grad_i, lapl_pair = _pair_kernel_terms(
                    problem.ring.bc,
                    R[i],
                    R[j];
                    pair_power=pair_power,
                    pair_eps=pair_eps,
                    smooth_pair_strength=smooth_pair_strength,
                    smooth_pair_softening=smooth_pair_softening,
                )
                logabs += logterm
                grad[i] += grad_i
                grad[j] -= grad_i
                lapllog += lapl_pair
            end
        end

        return logabs, grad, lapllog, 1.0
    end

    return TrialWF(
        R -> _trial_terms(R)[1],
        R -> _trial_terms(R)[2],
        R -> _trial_terms(R)[3],
        R -> _trial_terms(R)[4],
    )
end

function bosonic_tg_scaffold_trial_wavefunction(
    problem::PeriodicIonRingBosonProblem;
    node_tol::Real=1.0e-10,
    smooth_pair_strength::Real=0.0,
    smooth_pair_softening::Real=max(problem.bb_softening, 0.25 * problem.ring.a),
)
    function _trial_terms(R)
        logabs, grad, lapllog, _ = determinant_trial_terms(
            problem.ring,
            R,
            problem.N;
            node_tol=node_tol,
        )
        grad_out = copy(grad)
        total_lapl = lapllog

        @inbounds for i in 1:(problem.N - 1)
            for j in (i + 1):problem.N
                logterm, grad_i, lapl_pair = _pair_kernel_terms(
                    problem.ring.bc,
                    R[i],
                    R[j];
                    pair_power=0.0,
                    smooth_pair_strength=smooth_pair_strength,
                    smooth_pair_softening=smooth_pair_softening,
                )
                logabs += logterm
                grad_out[i] += grad_i
                grad_out[j] -= grad_i
                total_lapl += lapl_pair
            end
        end

        return logabs, grad_out, total_lapl, 1.0
    end

    return TrialWF(
        R -> _trial_terms(R)[1],
        R -> _trial_terms(R)[2],
        R -> _trial_terms(R)[3],
        R -> _trial_terms(R)[4],
    )
end

function fermion_determinant_trial_wavefunction(
    model::PeriodicIonRingModel,
    N::Integer;
    node_tol::Real=1.0e-10,
)
    N_int = Int(N)
    1 <= N_int <= size(model.coeffs, 2) || throw(ArgumentError("Need 1 <= N <= $(size(model.coeffs, 2)), got $N"))
    tol = Float64(node_tol)
    tol > 0 || throw(ArgumentError("node_tol must be > 0, got $node_tol"))

    function _fermion_terms(R)
        return determinant_trial_terms(model, R, N_int; node_tol=tol)
    end

    return TrialWF(
        R -> _fermion_terms(R)[1],
        R -> _fermion_terms(R)[2],
        R -> _fermion_terms(R)[3],
        R -> _fermion_terms(R)[4],
    )
end

function sample_uniform_ring_configurations(
    N::Integer,
    L::Real,
    nwalkers::Integer,
    rng::AbstractRNG;
    min_separation::Real=0.0,
)
    N_int = Int(N)
    nwalkers_int = Int(nwalkers)
    L_f = Float64(L)
    min_sep_f = Float64(min_separation)

    N_int >= 1 || throw(ArgumentError("N must be >= 1, got $N"))
    nwalkers_int >= 1 || throw(ArgumentError("nwalkers must be >= 1, got $nwalkers"))
    L_f > 0 || throw(ArgumentError("L must be > 0, got $L"))
    min_sep_f >= 0 || throw(ArgumentError("min_separation must be >= 0, got $min_separation"))

    bc = PeriodicBoundary1D(0.0, L_f)
    X = Matrix{Float64}(undef, N_int, nwalkers_int)

    @inbounds for w in 1:nwalkers_int
        accepted = false
        for _ in 1:10_000
            for i in 1:N_int
                X[i, w] = rand(rng) * L_f
            end
            if min_sep_f == 0.0
                accepted = true
                break
            end
            ok = true
            for i in 1:(N_int - 1)
                for j in (i + 1):N_int
                    if distance_1d(bc, X[i, w], X[j, w]) <= min_sep_f
                        ok = false
                        break
                    end
                end
                ok || break
            end
            if ok
                accepted = true
                break
            end
        end
        accepted || error("Failed to sample a valid periodic-ring configuration after 10000 attempts")
    end

    return X
end

function sample_boson_ring_configurations(
    problem::PeriodicIonRingBosonProblem,
    nwalkers::Integer,
    rng::AbstractRNG;
    min_separation::Union{Nothing,Real}=nothing,
)
    sep = min_separation === nothing ? problem.hard_core_radius : Float64(min_separation)
    return sample_uniform_ring_configurations(
        problem.N,
        problem.ring.L,
        nwalkers,
        rng;
        min_separation=sep,
    )
end

exact_manybody_energy(model::PeriodicIonRingModel, N::Integer) = sum(@view model.energies[1:Int(N)])
