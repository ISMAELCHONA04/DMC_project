# UseCases/common: VMC-backed initialization helpers for other methods

function _positions_from_walkers(walkers::AbstractVector{<:Walker})
    return [copy(position(w)) for w in walkers]
end

function _validate_vmc_targetN(params::VMCParams, count::Integer)
    count_int = Int(count)
    count_int == params.targetN || throw(ArgumentError(
        "VMC walker count mismatch: expected params.targetN=$(params.targetN), got $count_int",
    ))
    return count_int
end

function _validate_vmc_gfmc_targetN(vmc_params::VMCParams, gfmc_params::GFMCParams)
    vmc_params.targetN == gfmc_params.targetN || throw(ArgumentError(
        "VMC/GFMC targetN mismatch: vmc targetN=$(vmc_params.targetN), gfmc targetN=$(gfmc_params.targetN)",
    ))
    return gfmc_params.targetN
end

"""
    sample_vmc_positions(H, params, configs, trial_wf; kwargs...) -> Vector{Vector{Float64}}

Run a short VMC equilibration and return the final walker positions as raw
configuration vectors. This is intended as a warm-start helper for projector
methods such as GFMC.

Unlike `run_vmc`, this helper does not record the full VMC position history by
default. It only keeps the final evolved walker ensemble and returns deep
copies of those positions.
"""
function sample_vmc_positions(
    h,
    params::VMCParams,
    configs,
    trial_wf;
    rng::AbstractRNG=Random.default_rng(),
    proposal::Union{AbstractVMCProposal,Symbol}=DriftGaussianProposal(),
    step::Integer=0,
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout,
)
    sim = VMCSim(h, params, configs, trial_wf, rng; proposal=proposal, step=step)
    _validate_vmc_targetN(params, length(sim.walkers))

    run_vmc!(
        sim;
        snapshot_steps=[params.nsteps],
        show_progress=show_progress,
        progress_every=progress_every,
        progress_label=progress_label,
        debug=debug,
        debug_every=debug_every,
        debug_io=debug_io,
    )
    return _positions_from_walkers(sim.walkers)
end

"""
    run_gfmc_with_vmc_init(H, gfmc_params, init_configs, trial_wf, vmc_params; kwargs...) -> GFMCSim

Warm-start a generic Hamiltonian-based GFMC run by first sampling an initial
walker ensemble with VMC from the same trial wavefunction.

The VMC stage samples `|psi_T|^2`; the returned positions are then passed into
`GFMCSim(...)` and evolved by `run_gfmc!`.
"""
function run_gfmc_with_vmc_init(
    h,
    gfmc_params::GFMCParams,
    init_configs,
    trial_wf,
    vmc_params::VMCParams;
    vmc_rng::AbstractRNG=Random.default_rng(),
    gfmc_rng::AbstractRNG=Random.default_rng(),
    proposal::Union{AbstractVMCProposal,Symbol}=DriftGaussianProposal(),
    guiding::AbstractGuiding=ImportanceGuiding(trial_wf, h),
    nodepolicy::AbstractNodePolicy=NoNode(),
    backend::AbstractGFMCBackend=CPUBackend(),
    reconfiguration::AbstractGFMCReconfiguration=SystematicReconfiguration(),
    step::Integer=0,
    vmc_step::Integer=0,
    vmc_show_progress::Bool=false,
    vmc_progress_every::Integer=0,
    vmc_progress_label::AbstractString="",
    vmc_debug::Bool=false,
    vmc_debug_every::Integer=1,
    vmc_debug_io::IO=stdout,
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout,
)
    _validate_vmc_gfmc_targetN(vmc_params, gfmc_params)

    positions = sample_vmc_positions(
        h,
        vmc_params,
        init_configs,
        trial_wf;
        rng=vmc_rng,
        proposal=proposal,
        step=vmc_step,
        show_progress=vmc_show_progress,
        progress_every=vmc_progress_every,
        progress_label=vmc_progress_label,
        debug=vmc_debug,
        debug_every=vmc_debug_every,
        debug_io=vmc_debug_io,
    )

    sim = GFMCSim(
        h,
        gfmc_params,
        positions,
        gfmc_rng;
        guiding=guiding,
        nodepolicy=nodepolicy,
        backend=backend,
        reconfiguration=reconfiguration,
        step=step,
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

"""
    run_gfmc_with_vmc_init(model::SpinlessFermionRing1D, gfmc_params, init_configs, vmc_params; kwargs...) -> GFMCSim

Warm-start a `SpinlessFermionRing1D` GFMC run by first sampling initial walker
positions with VMC using the model's built-in trial wavefunction, then running
the specialized fermion-ring GFMC kernel.

If you need a custom VMC trial state for the warm start, use the generic
Hamiltonian overload with `hamiltonian(model)` and your custom `TrialWF`.
"""
function run_gfmc_with_vmc_init(
    model::SpinlessFermionRing1D,
    gfmc_params::GFMCParams,
    init_configs,
    vmc_params::VMCParams;
    vmc_rng::AbstractRNG=Random.default_rng(),
    gfmc_rng::AbstractRNG=Random.default_rng(),
    proposal::Union{AbstractVMCProposal,Symbol}=DriftGaussianProposal(),
    use_guiding::Bool=true,
    nodepolicy::AbstractNodePolicy=FixedNode(),
    backend::AbstractGFMCBackend=CPUBackend(),
    reconfiguration::AbstractGFMCReconfiguration=SystematicReconfiguration(),
    step::Integer=0,
    vmc_step::Integer=0,
    vmc_show_progress::Bool=false,
    vmc_progress_every::Integer=0,
    vmc_progress_label::AbstractString="",
    vmc_debug::Bool=false,
    vmc_debug_every::Integer=1,
    vmc_debug_io::IO=stdout,
    snapshot_steps::AbstractVector{<:Integer}=Int[],
    show_progress::Bool=false,
    progress_every::Integer=0,
    progress_label::AbstractString="",
    debug::Bool=false,
    debug_every::Integer=1,
    debug_io::IO=stdout,
)
    _validate_vmc_gfmc_targetN(vmc_params, gfmc_params)

    h = hamiltonian(model)
    trial_wf = trial_wavefunction(model)
    positions = sample_vmc_positions(
        h,
        vmc_params,
        init_configs,
        trial_wf;
        rng=vmc_rng,
        proposal=proposal,
        step=vmc_step,
        show_progress=vmc_show_progress,
        progress_every=vmc_progress_every,
        progress_label=vmc_progress_label,
        debug=vmc_debug,
        debug_every=vmc_debug_every,
        debug_io=vmc_debug_io,
    )

    sim = GFMCSim(
        model,
        gfmc_params,
        positions,
        gfmc_rng;
        use_guiding=use_guiding,
        nodepolicy=nodepolicy,
        backend=backend,
        reconfiguration=reconfiguration,
        step=step,
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
