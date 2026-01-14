# Domain: Hamiltonian for 1D systems

struct Hamiltonian{V}
    N::Int       # number of particles
    D::Float64   # diffusion constant
    V::V         # potential function: R -> V(R)
end



# Accessor methods
potential(H::Hamiltonian, R) = H.V(R)
nparticles(H::Hamiltonian) = H.N
diffusion_constant(H::Hamiltonian) = H.D

# Guiding potential
abstract type AbstractGuiding end
struct NoGuiding <: AbstractGuiding end

# Node policy
abstract type AbstractNodePolicy end
struct NoNode <: AbstractNodePolicy end
struct FixedNode <: AbstractNodePolicy end

# Guiding / trial wavefunction 
struct TrialWF{LP,GL,LL,SP}
    logpsi::LP # R -> log|psi_TT(R)|
    gradlogpsi::GL # R -> ∇ log|psi_T(R)| (Vector)
    lapllogpsi::LL # R -> ∇² log|psi_T(R)| (Scalar)
    signpsi::SP # R -> sign(psi_T(R)) ∈ {-1, 0, 1}
end

TrialWF(logpsi, gradlogpsi, lapllogpsi) = TrialWF(logpsi, gradlogpsi, lapllogpsi, R -> 1.0)

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
