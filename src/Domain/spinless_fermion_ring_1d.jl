# Domain: Spinless fermions on a 1D periodic ring with a cosine lattice

"""
    SpinlessFermionRing1D(N, a, L, V0; kwargs...)

Built-in periodic-ring model for spinless fermions in a cosine lattice. The
model packages the Hamiltonian, analytic trial state, drift, and fixed-node
sign used by the fermion-ring DMC and GFMC experiments.
"""
struct SpinlessFermionRing1D{BC<:PeriodicBoundary1D}
    N::Int
    a::Float64
    L::Float64
    V0::Float64
    lambda::Float64
    D::Float64
    node_tol::Float64
    trig_eps::Float64
    bc::BC
    k_lat::Float64
    alpha_pair::Float64
end

function SpinlessFermionRing1D(
    N::Integer,
    a::Real,
    L::Real,
    V0::Real;
    lambda::Union{Nothing,Real}=nothing,
    D::Real=0.5,
    node_tol::Real=1e-7,
    trig_eps::Real=1e-10,
    twist::Union{Real,AbstractVector{<:Real}}=0.0,
)
    N_int = Int(N)
    N_int >= 1 || throw(ArgumentError("N must be >= 1, got $N"))
    a_f = Float64(a)
    L_f = Float64(L)
    a_f > 0 || throw(ArgumentError("a must be > 0, got $a"))
    L_f > 0 || throw(ArgumentError("L must be > 0, got $L"))

    lambda_f = lambda === nothing ? (-0.5 * Float64(V0)) : Float64(lambda)
    bc = PeriodicBoundary1D(0.0, L_f; twist=twist)
    return SpinlessFermionRing1D(
        N_int,
        a_f,
        L_f,
        Float64(V0),
        lambda_f,
        Float64(D),
        Float64(node_tol),
        Float64(trig_eps),
        bc,
        2pi / a_f,
        pi / L_f,
    )
end

nparticles(model::SpinlessFermionRing1D) = model.N
diffusion_constant(model::SpinlessFermionRing1D) = model.D
boundary(model::SpinlessFermionRing1D) = model.bc

function _assert_ring_configuration(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    length(R) == model.N || throw(ArgumentError("Expected configuration length $(model.N), got $(length(R))"))
    return nothing
end

function _wrapped_ring_configuration(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    _assert_ring_configuration(model, R)
    return wrap_position(model.bc, R)
end

function _ring_trial_terms(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    Rw = _wrapped_ring_configuration(model, R)
    grad = zeros(Float64, model.N)
    lapl = 0.0
    logabs = 0.0
    sign = 1.0

    @inbounds for i in 1:model.N
        xi = Rw[i]
        coskx = cos(model.k_lat * xi)
        sinkx = sin(model.k_lat * xi)
        logabs += model.lambda * coskx
        grad[i] = -model.lambda * model.k_lat * sinkx
        lapl += -model.lambda * model.k_lat^2 * coskx
    end

    @inbounds for i in 1:(model.N - 1)
        for j in (i + 1):model.N
            dx_ij = displacement(model.bc, Rw[j], Rw[i])
            u = model.alpha_pair * dx_ij
            s = sin(u)
            abs_s = abs(s)

            if abs_s < model.node_tol
                sign = 0.0
            elseif sign != 0.0 && s < 0
                sign = -sign
            end

            s_eff = abs_s < model.trig_eps ? (s >= 0 ? model.trig_eps : -model.trig_eps) : s
            logabs += log(abs(s_eff))

            cot_u = cos(u) / s_eff
            c = model.alpha_pair * cot_u
            grad[i] += c
            grad[j] -= c

            inv_s2 = inv(s_eff * s_eff)
            lapl += -2 * (model.alpha_pair^2) * inv_s2
        end
    end

    return logabs, grad, lapl, sign
end

function potential(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    Rw = _wrapped_ring_configuration(model, R)
    V = 0.0
    @inbounds for i in 1:model.N
        V += model.V0 * cos(model.k_lat * Rw[i])
    end
    return V
end

function trial_logpsi(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    logabs, _, _, _ = _ring_trial_terms(model, R)
    return logabs
end

function trial_gradlogpsi(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    _, grad, _, _ = _ring_trial_terms(model, R)
    return grad
end

function trial_lapllogpsi(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    _, _, lapl, _ = _ring_trial_terms(model, R)
    return lapl
end

function signpsi(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    _, _, _, sign = _ring_trial_terms(model, R)
    return sign
end

function drift(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    _, grad, _, _ = _ring_trial_terms(model, R)
    return 2.0 .* grad
end

function local_energy(model::SpinlessFermionRing1D, R::AbstractVector{<:Real})
    _, grad, lapl, _ = _ring_trial_terms(model, R)
    return -model.D * (lapl + sum(abs2, grad)) + potential(model, R)
end

"""Return a generic `Hamiltonian` view of the fermion-ring model."""
function hamiltonian(model::SpinlessFermionRing1D)
    return Hamiltonian(model.N, model.D, R -> potential(model, R), model.bc)
end

"""Return the built-in trial wavefunction used by the fermion-ring workflows."""
function trial_wavefunction(model::SpinlessFermionRing1D)
    return TrialWF(
        R -> trial_logpsi(model, R),
        R -> trial_gradlogpsi(model, R),
        R -> trial_lapllogpsi(model, R),
        R -> signpsi(model, R),
    )
end

"""Return the model's default `ImportanceGuiding` wrapper."""
function importance_guiding(model::SpinlessFermionRing1D)
    return ImportanceGuiding(trial_wavefunction(model), hamiltonian(model))
end

function sample_uniform_configurations(
    model::SpinlessFermionRing1D,
    nwalkers::Integer,
    rng::AbstractRNG,
)
    nwalkers_int = Int(nwalkers)
    nwalkers_int >= 1 || throw(ArgumentError("nwalkers must be >= 1, got $nwalkers"))

    X = Matrix{Float64}(undef, model.N, nwalkers_int)
    @inbounds for j in 1:nwalkers_int
        for i in 1:model.N
            X[i, j] = rand(rng) * model.L
        end
    end
    return X
end
