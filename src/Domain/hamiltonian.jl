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
