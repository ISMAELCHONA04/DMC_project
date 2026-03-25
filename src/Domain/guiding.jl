# Domain: Guiding wavefunctions and importance sampling

abstract type AbstractGuiding end

"""Guiding policy used for unguided projector dynamics."""
struct NoGuiding <: AbstractGuiding end

"""
    ImportanceGuiding(trial, H)

Importance-sampling wrapper that combines a `TrialWF` with a `Hamiltonian`.
The method layers use it to compute drift vectors, local energies, and trial
signs.
"""
struct ImportanceGuiding{TW,H} <: AbstractGuiding
    trial::TW
    H::H
end

"""Return the quantum-force drift vector `2∇log|ψ_T|` for configuration `R`."""
drift(g::ImportanceGuiding, R) = 2 .* g.trial.gradlogpsi(R)

# Sign of the trial wavefunction for node handling.
"""
    signpsi(g, R) -> Real

Return the sign of the trial wavefunction at configuration `R`. Guiding
policies without a sign default to `+1`. For fixed-node, return `0` near
nodes (e.g., within a small tolerance).
"""
signpsi(::NoGuiding, R) = 1.0
signpsi(t::TrialWF, R) = t.signpsi(R)
signpsi(g::ImportanceGuiding, R) = signpsi(g.trial, R)

"""
    local_energy(g, R) -> Real

Return the trial-state local energy `(Hψ_T)/ψ_T` for configuration `R`.
"""
function local_energy(g::ImportanceGuiding, R)
    D = diffusion_constant(g.H)
    gradlog = g.trial.gradlogpsi(R)
    lapllog = g.trial.lapllogpsi(R)
    V = potential(g.H, R)
    return -D * (lapllog + sum(abs2, gradlog)) + V
end
