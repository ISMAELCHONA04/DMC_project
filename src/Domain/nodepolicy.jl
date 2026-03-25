# Domain: Node policies for fixed-node diffusion Monte Carlo

abstract type AbstractNodePolicy end

"""Node policy used when no fixed-node rejection is required."""
struct NoNode <: AbstractNodePolicy end

"""Node policy for fixed-node rejection against the trial-wavefunction sign."""
struct FixedNode <: AbstractNodePolicy end
