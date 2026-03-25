# UseCases/gfmc: Abstract kernel, backend, and resampling interfaces

abstract type AbstractGFMCKernel end
abstract type AbstractGFMCBackend end
abstract type AbstractGFMCReconfiguration end

struct CPUBackend <: AbstractGFMCBackend end

struct MultinomialReconfiguration <: AbstractGFMCReconfiguration end
struct SystematicReconfiguration <: AbstractGFMCReconfiguration end

configuration_dimension(::AbstractGFMCKernel) = error("configuration_dimension not implemented")
diffusion_constant(::AbstractGFMCKernel) = error("diffusion_constant not implemented")
boundary(::AbstractGFMCKernel) = error("boundary not implemented")
uses_importance_sampling(::AbstractGFMCKernel) = false
uses_fixed_node(::AbstractGFMCKernel) = false

allocate_matrix(::CPUBackend, nrows::Integer, ncols::Integer) = Matrix{Float64}(undef, Int(nrows), Int(ncols))
allocate_vector(::CPUBackend, n::Integer) = Vector{Float64}(undef, Int(n))

function randn!(::CPUBackend, rng::AbstractRNG, X::AbstractMatrix{<:Real})
    for idx in eachindex(X)
        X[idx] = randn(rng)
    end
    return X
end

function fill_ones!(::CPUBackend, v::AbstractVector{<:Real})
    fill!(v, 1.0)
    return v
end

function copy_columns!(
    ::CPUBackend,
    dest::AbstractMatrix{<:Real},
    src::AbstractMatrix{<:Real},
    indices::AbstractVector{<:Integer},
)
    size(dest, 2) == length(indices) || throw(ArgumentError("Destination column count and index count must match"))
    size(dest, 1) == size(src, 1) || throw(ArgumentError("Source/destination row counts must match"))
    @inbounds for j in eachindex(indices)
        src_col = Int(indices[j])
        for i in axes(dest, 1)
            dest[i, j] = src[i, src_col]
        end
    end
    return dest
end

function copy_values!(
    ::CPUBackend,
    dest::AbstractVector{<:Real},
    src::AbstractVector{<:Real},
    indices::AbstractVector{<:Integer},
)
    length(dest) == length(indices) || throw(ArgumentError("Destination length and index count must match"))
    @inbounds for j in eachindex(indices)
        dest[j] = src[Int(indices[j])]
    end
    return dest
end
