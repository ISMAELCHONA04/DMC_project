# Domain: Guiding wavefunctions and importance sampling

abstract type AbstractGuiding end

struct NoGuiding <: AbstractGuiding end

# Importance-sampling wrapper: combines trial WF with Hamiltonian.
struct ImportanceGuiding{TW,H} <: AbstractGuiding
    trial::TW
    H::H
end

# Quantum force (drift) used in drift-diffusion proposals.
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

# Local energy EL = (HψT)/ψT, used in branching.
function local_energy(g::ImportanceGuiding, R)
    D = diffusion_constant(g.H)
    gradlog = g.trial.gradlogpsi(R)
    lapllog = g.trial.lapllogpsi(R)
    V = potential(g.H, R)
    return -D * (lapllog + sum(abs2, gradlog)) + V
end
