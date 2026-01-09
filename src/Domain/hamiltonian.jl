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

# Guiding / trial wavefunction 
struct TrialWF{LP,GL,LL}
    logpsi::LP # R -> log|psi_TT(R)|
    gradlogpsi::GL # R -> ∇ log|psi_T(R)| (Vector)
    lapllogpsi::LL # R -> ∇² log|psi_T(R)| (Scalar)
end

# Importance-sampling wrapper: combines trial WF with Hamiltonian.
struct ImportanceGuiding{TW,H} <: AbstractGuiding
    trial::TW
    H::H
end

# Quantum force (drift) used in drift-diffusion proposals.
drift(g::ImportanceGuiding, R) = 2 .* g.trial.gradlogpsi(R)

# Local energy EL = (HψT)/ψT, used in branching.
function local_energy(g::ImportanceGuiding, R)
    D = diffusion_constant(g.H)
    gradlog = g.trial.gradlogpsi(R)
    lapllog = g.trial.lapllogpsi(R)
    V = potential(g.H, R)
    return -D * (lapllog + sum(abs2, gradlog)) + V
end
