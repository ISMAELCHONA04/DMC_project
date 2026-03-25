# UseCases/gfmc: Fixed-population reconfiguration schemes

function _normalized_weights(weights::AbstractVector{<:Real})
    total = sum(weights)
    isfinite(total) && total > 0 || error("All GFMC walker weights vanished.")
    probs = Vector{Float64}(undef, length(weights))
    @inbounds for i in eachindex(weights)
        wi = Float64(weights[i])
        wi >= 0 || throw(ArgumentError("GFMC weights must be non-negative"))
        probs[i] = wi / total
    end
    return probs, total / length(weights)
end

function resample_indices(
    ::MultinomialReconfiguration,
    rng::AbstractRNG,
    weights::AbstractVector{<:Real},
)
    probs, mean_weight = _normalized_weights(weights)
    cumulative = cumsum(probs)
    out = Vector{Int}(undef, length(weights))
    @inbounds for i in eachindex(out)
        u = rand(rng)
        out[i] = searchsortedfirst(cumulative, u)
    end
    return out, mean_weight
end

function resample_indices(
    ::SystematicReconfiguration,
    rng::AbstractRNG,
    weights::AbstractVector{<:Real},
)
    probs, mean_weight = _normalized_weights(weights)
    cumulative = cumsum(probs)
    n = length(weights)
    out = Vector{Int}(undef, n)
    u0 = rand(rng) / n
    j = 1
    @inbounds for i in 1:n
        u = u0 + (i - 1) / n
        while j < n && u > cumulative[j]
            j += 1
        end
        out[i] = j
    end
    return out, mean_weight
end
