# UseCases/gfmc: Batched kernel evaluation

struct GenericGFMCKernel{H,G<:AbstractGuiding,NP<:AbstractNodePolicy} <: AbstractGFMCKernel
    H::H
    guiding::G
    nodepolicy::NP
end

function GenericGFMCKernel(
    ham::Ham;
    guiding::G=NoGuiding(),
    nodepolicy::NP=NoNode(),
) where {Ham,G<:AbstractGuiding,NP<:AbstractNodePolicy}
    return GenericGFMCKernel{Ham,G,NP}(ham, guiding, nodepolicy)
end

struct SpinlessFermionRing1DKernel{NP<:AbstractNodePolicy} <: AbstractGFMCKernel
    model::SpinlessFermionRing1D
    use_guiding::Bool
    nodepolicy::NP
end

function SpinlessFermionRing1DKernel(
    model::SpinlessFermionRing1D;
    use_guiding::Bool=true,
    nodepolicy::NP=FixedNode(),
) where {NP<:AbstractNodePolicy}
    return SpinlessFermionRing1DKernel{NP}(model, use_guiding, nodepolicy)
end

mutable struct GFMCBatchData{M<:AbstractMatrix{Float64},V<:AbstractVector{Float64}}
    drift::M
    logpsi::V
    branch_energy::V
    sign::V
end

function GFMCBatchData(backend::AbstractGFMCBackend, ncoords::Integer, nwalkers::Integer)
    return GFMCBatchData(
        allocate_matrix(backend, ncoords, nwalkers),
        allocate_vector(backend, nwalkers),
        allocate_vector(backend, nwalkers),
        allocate_vector(backend, nwalkers),
    )
end

configuration_dimension(kernel::GenericGFMCKernel) = nparticles(kernel.H)
diffusion_constant(kernel::GenericGFMCKernel) = diffusion_constant(kernel.H)
boundary(kernel::GenericGFMCKernel) = boundary(kernel.H)
uses_importance_sampling(kernel::GenericGFMCKernel) = kernel.guiding isa ImportanceGuiding
uses_fixed_node(kernel::GenericGFMCKernel) = kernel.nodepolicy isa FixedNode

configuration_dimension(kernel::SpinlessFermionRing1DKernel) = nparticles(kernel.model)
diffusion_constant(kernel::SpinlessFermionRing1DKernel) = diffusion_constant(kernel.model)
boundary(kernel::SpinlessFermionRing1DKernel) = boundary(kernel.model)
uses_importance_sampling(kernel::SpinlessFermionRing1DKernel) = kernel.use_guiding
uses_fixed_node(kernel::SpinlessFermionRing1DKernel) = kernel.nodepolicy isa FixedNode

function wrap_configurations!(kernel::AbstractGFMCKernel, X::AbstractMatrix{<:Real})
    bc = boundary(kernel)
    @inbounds for j in axes(X, 2)
        for i in axes(X, 1)
            X[i, j] = wrap_coordinate(bc, X[i, j])
        end
    end
    return X
end

function evaluate_configuration_data!(kernel::GenericGFMCKernel, data::GFMCBatchData, X::AbstractMatrix{<:Real})
    uses_guiding = uses_importance_sampling(kernel)
    @inbounds for j in axes(X, 2)
        R = @view X[:, j]
        if uses_guiding
            drift_j = drift(kernel.guiding, R)
            for i in axes(X, 1)
                data.drift[i, j] = Float64(drift_j[i])
            end
            data.logpsi[j] = Float64(kernel.guiding.trial.logpsi(R))
            data.branch_energy[j] = Float64(local_energy(kernel.guiding, R))
            data.sign[j] = Float64(signpsi(kernel.guiding, R))
        else
            for i in axes(X, 1)
                data.drift[i, j] = 0.0
            end
            data.logpsi[j] = 0.0
            data.branch_energy[j] = Float64(potential(kernel.H, R))
            data.sign[j] = 1.0
        end
    end
    return data
end

function evaluate_configuration_data!(kernel::SpinlessFermionRing1DKernel, data::GFMCBatchData, X::AbstractMatrix{<:Real})
    model = kernel.model
    ncoords = size(X, 1)
    ncoords == model.N || throw(ArgumentError("Expected $(model.N) coordinates per configuration, got $ncoords"))

    @inbounds for j in axes(X, 2)
        logpsi_j = 0.0
        lapl_j = 0.0
        potential_j = 0.0
        sign_j = 1.0

        if kernel.use_guiding
            for i in 1:model.N
                xi = X[i, j]
                coskx = cos(model.k_lat * xi)
                sinkx = sin(model.k_lat * xi)
                logpsi_j += model.lambda * coskx
                data.drift[i, j] = -2.0 * model.lambda * model.k_lat * sinkx
                lapl_j += -model.lambda * model.k_lat^2 * coskx
                potential_j += model.V0 * coskx
            end

            for i in 1:(model.N - 1)
                for k in (i + 1):model.N
                    dx_ik = displacement(model.bc, X[k, j], X[i, j])
                    u = model.alpha_pair * dx_ik
                    s = sin(u)
                    abs_s = abs(s)

                    if abs_s < model.node_tol
                        sign_j = 0.0
                    elseif sign_j != 0.0 && s < 0
                        sign_j = -sign_j
                    end

                    s_eff = abs_s < model.trig_eps ? (s >= 0 ? model.trig_eps : -model.trig_eps) : s
                    logpsi_j += log(abs(s_eff))

                    c = model.alpha_pair * cos(u) / s_eff
                    data.drift[i, j] += 2.0 * c
                    data.drift[k, j] -= 2.0 * c

                    lapl_j += -2.0 * (model.alpha_pair^2) / (s_eff * s_eff)
                end
            end

            grad_norm2 = 0.0
            for i in 1:model.N
                grad_i = 0.5 * data.drift[i, j]
                grad_norm2 += grad_i * grad_i
            end

            data.logpsi[j] = logpsi_j
            data.branch_energy[j] = -model.D * (lapl_j + grad_norm2) + potential_j
            data.sign[j] = sign_j
        else
            for i in 1:model.N
                xi = X[i, j]
                potential_j += model.V0 * cos(model.k_lat * xi)
                data.drift[i, j] = 0.0
            end
            data.logpsi[j] = 0.0
            data.branch_energy[j] = potential_j
            data.sign[j] = 1.0
        end
    end
    return data
end

function copy_batch_data!(
    backend::AbstractGFMCBackend,
    dest::GFMCBatchData,
    src::GFMCBatchData,
    indices::AbstractVector{<:Integer},
)
    copy_columns!(backend, dest.drift, src.drift, indices)
    copy_values!(backend, dest.logpsi, src.logpsi, indices)
    copy_values!(backend, dest.branch_energy, src.branch_energy, indices)
    copy_values!(backend, dest.sign, src.sign, indices)
    return dest
end
