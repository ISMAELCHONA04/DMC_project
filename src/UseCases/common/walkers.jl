# UseCases/common: Walker coercion and boundary wrapping helpers

function _walker_from_position(R::AbstractVector{<:Real})
    return Walker(Float64[x for x in R])
end

function _coerce_walkers(walkers::AbstractVector{<:Walker})
    isempty(walkers) && throw(ArgumentError("Need at least one walker"))
    return [clone(w) for w in walkers]
end

function _coerce_walkers(positions::AbstractVector{<:AbstractVector{<:Real}})
    isempty(positions) && throw(ArgumentError("Need at least one walker configuration"))
    return [_walker_from_position(R) for R in positions]
end

function _coerce_walkers(position::AbstractVector{<:Real})
    return Walker[_walker_from_position(position)]
end

function _coerce_walkers(positions::AbstractMatrix{<:Real})
    size(positions, 2) >= 1 || throw(ArgumentError("Need at least one walker configuration"))
    walkers = Vector{Walker}(undef, size(positions, 2))
    @inbounds for j in axes(positions, 2)
        walkers[j] = _walker_from_position(@view positions[:, j])
    end
    return walkers
end

function _wrapped_walker(h, w::Walker)
    Rw, phase_factor = wrap_position_with_phase(boundary(h), position(w))
    return Walker(Rw, phase(w) * phase_factor)
end

function _wrapped_walkers(h, configs)
    return [_wrapped_walker(h, w) for w in _coerce_walkers(configs)]
end
