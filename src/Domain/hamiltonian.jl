# Domain: Hamiltonian for 1D systems

"""
    Hamiltonian(N, D, V, bc=OpenBoundary1D())

Container for a one-dimensional many-body Hamiltonian represented by:

- `N`: configuration dimension / particle count
- `D`: diffusion constant used by projector methods
- `V`: potential callback with signature `R -> scalar`
- `bc`: boundary-condition policy
"""
struct Hamiltonian{V,BC<:AbstractBoundary1D}
    N::Int       # number of particles
    D::Float64   # diffusion constant
    V::V         # potential function: R -> V(R)
    bc::BC       # boundary-condition policy
end

# Backward-compatible constructor: defaults to open boundaries.
function Hamiltonian(N::Int, D::Real, V, bc::AbstractBoundary1D=OpenBoundary1D())
    return Hamiltonian{typeof(V),typeof(bc)}(N, Float64(D), V, bc)
end

"""Return the potential energy for configuration `R`."""
potential(H::Hamiltonian, R) = H.V(R)
"""Return the configuration dimension / particle count."""
nparticles(H::Hamiltonian) = H.N
"""Return the diffusion constant used by the stochastic propagators."""
diffusion_constant(H::Hamiltonian) = H.D
"""Return the boundary-condition policy attached to the Hamiltonian."""
boundary(H::Hamiltonian) = H.bc
