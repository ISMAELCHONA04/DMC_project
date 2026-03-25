# UseCases/gfmc: Fixed-population GFMC simulation container and state

mutable struct GFMCSim{
    K<:AbstractGFMCKernel,
    RNG,
    B<:AbstractGFMCBackend,
    R<:AbstractGFMCReconfiguration,
    M<:AbstractMatrix{Float64},
    V<:AbstractVector{Float64},
}
    kernel::K
    params::GFMCParams
    positions::M
    proposal_positions::M
    state_data::GFMCBatchData{M,V}
    proposal_data::GFMCBatchData{M,V}
    weights::V
    rng::RNG
    backend::B
    reconfiguration::R
    step::Int
    ET::Float64
    ET_history::Vector{Float64}
    population_history::Vector{Int}
    energy_mean_history::Vector{Float64}
    energy_variance_history::Vector{Float64}
    effective_population_history::Vector{Float64}
    mean_weight_history::Vector{Float64}
    acceptance_history::Vector{Float64}
    walker_positions_history::Vector{Vector{Vector{Float64}}}
end

function _position_matrix(positions::AbstractMatrix{<:Real})
    X = Matrix{Float64}(undef, size(positions, 1), size(positions, 2))
    @inbounds for j in axes(X, 2)
        for i in axes(X, 1)
            X[i, j] = Float64(positions[i, j])
        end
    end
    return X
end

function _position_matrix(positions::AbstractVector{<:AbstractVector{<:Real}})
    isempty(positions) && throw(ArgumentError("Need at least one walker configuration"))
    ncoords = length(first(positions))
    nwalkers = length(positions)
    X = Matrix{Float64}(undef, ncoords, nwalkers)
    @inbounds for j in 1:nwalkers
        pos = positions[j]
        length(pos) == ncoords || throw(ArgumentError("All walker configurations must have the same length"))
        for i in 1:ncoords
            X[i, j] = Float64(pos[i])
        end
    end
    return X
end

function _position_matrix(walkers::AbstractVector{<:Walker})
    return _position_matrix([position(w) for w in walkers])
end

function _snapshot_from_matrix(X::AbstractMatrix{<:Real})
    snap = Vector{Vector{Float64}}(undef, size(X, 2))
    @inbounds for j in axes(X, 2)
        snap[j] = collect(Float64, @view X[:, j])
    end
    return snap
end

function _effective_population(weights::AbstractVector{<:Real})
    s1 = 0.0
    s2 = 0.0
    @inbounds for w in weights
        wf = Float64(w)
        s1 += wf
        s2 += wf * wf
    end
    return s2 > 0 ? (s1 * s1 / s2) : 0.0
end

function estimate_energy(sim::GFMCSim)::Float64
    numerator = 0.0
    denom = 0.0
    @inbounds for j in eachindex(sim.weights)
        w = sim.weights[j]
        numerator += w * sim.state_data.branch_energy[j]
        denom += w
    end
    return denom > 0 ? numerator / denom : 0.0
end

function estimate_energy_variance(sim::GFMCSim)::Float64
    mean_E = estimate_energy(sim)
    numerator = 0.0
    denom = 0.0
    @inbounds for j in eachindex(sim.weights)
        w = sim.weights[j]
        diff = sim.state_data.branch_energy[j] - mean_E
        numerator += w * diff * diff
        denom += w
    end
    return denom > 0 ? numerator / denom : 0.0
end

function record_positions!(sim::GFMCSim)
    push!(sim.walker_positions_history, _snapshot_from_matrix(sim.positions))
    return sim
end

function record_state!(sim::GFMCSim, record_positions::Bool=true)
    push!(sim.population_history, size(sim.positions, 2))
    push!(sim.energy_mean_history, estimate_energy(sim))
    push!(sim.energy_variance_history, estimate_energy_variance(sim))
    push!(sim.ET_history, sim.ET)
    if record_positions
        record_positions!(sim)
    end
    return sim
end

function _reset_histories!(sim::GFMCSim)
    sim.step = 0
    sim.ET = sim.params.ET0
    sim.ET_history = Float64[]
    sim.population_history = Int[]
    sim.energy_mean_history = Float64[]
    sim.energy_variance_history = Float64[]
    sim.effective_population_history = Float64[]
    sim.mean_weight_history = Float64[]
    sim.acceptance_history = Float64[]
    sim.walker_positions_history = Vector{Vector{Vector{Float64}}}()
    return sim
end

function initialize!(sim::GFMCSim, positions)
    X = _position_matrix(positions)
    size(X, 2) == sim.params.targetN || throw(ArgumentError("Need exactly $(sim.params.targetN) walkers"))
    size(X, 1) == configuration_dimension(sim.kernel) || throw(ArgumentError("Configuration dimension mismatch"))
    sim.positions .= X
    wrap_configurations!(sim.kernel, sim.positions)
    evaluate_configuration_data!(sim.kernel, sim.state_data, sim.positions)
    fill_ones!(sim.backend, sim.weights)
    _reset_histories!(sim)
    return sim
end

function GFMCSim(
    kernel::K,
    params::GFMCParams,
    positions,
    rng::RNG;
    backend::B=CPUBackend(),
    reconfiguration::R=SystematicReconfiguration(),
    step::Integer=0,
) where {K<:AbstractGFMCKernel,RNG,B<:AbstractGFMCBackend,R<:AbstractGFMCReconfiguration}
    X = _position_matrix(positions)
    size(X, 2) == params.targetN || throw(ArgumentError("Need exactly $(params.targetN) walkers"))
    size(X, 1) == configuration_dimension(kernel) || throw(ArgumentError("Configuration dimension mismatch"))

    wrap_configurations!(kernel, X)
    proposal_positions = allocate_matrix(backend, size(X, 1), size(X, 2))
    proposal_positions .= X

    state_data = GFMCBatchData(backend, size(X, 1), size(X, 2))
    proposal_data = GFMCBatchData(backend, size(X, 1), size(X, 2))
    evaluate_configuration_data!(kernel, state_data, X)

    weights = allocate_vector(backend, size(X, 2))
    fill_ones!(backend, weights)

    sim = GFMCSim(
        kernel,
        params,
        X,
        proposal_positions,
        state_data,
        proposal_data,
        weights,
        rng,
        backend,
        reconfiguration,
        Int(step),
        params.ET0,
        Float64[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Vector{Vector{Vector{Float64}}}(),
    )
    _reset_histories!(sim)
    sim.step = Int(step)
    return sim
end

function GFMCSim(
    ham::Ham,
    params::GFMCParams,
    configs,
    rng::RNG;
    guiding::G=NoGuiding(),
    nodepolicy::NP=NoNode(),
    backend::B=CPUBackend(),
    reconfiguration::R=SystematicReconfiguration(),
    step::Integer=0,
) where {Ham,RNG,G<:AbstractGuiding,NP<:AbstractNodePolicy,B<:AbstractGFMCBackend,R<:AbstractGFMCReconfiguration}
    kernel = GenericGFMCKernel(ham; guiding=guiding, nodepolicy=nodepolicy)
    return GFMCSim(
        kernel,
        params,
        configs,
        rng;
        backend=backend,
        reconfiguration=reconfiguration,
        step=step,
    )
end

function GFMCSim(
    model::SpinlessFermionRing1D,
    params::GFMCParams,
    positions,
    rng::RNG;
    use_guiding::Bool=true,
    nodepolicy::NP=FixedNode(),
    backend::B=CPUBackend(),
    reconfiguration::R=SystematicReconfiguration(),
    step::Integer=0,
) where {NP<:AbstractNodePolicy,B<:AbstractGFMCBackend,R<:AbstractGFMCReconfiguration,RNG}
    kernel = SpinlessFermionRing1DKernel(model; use_guiding=use_guiding, nodepolicy=nodepolicy)
    return GFMCSim(
        kernel,
        params,
        positions,
        rng;
        backend=backend,
        reconfiguration=reconfiguration,
        step=step,
    )
end
