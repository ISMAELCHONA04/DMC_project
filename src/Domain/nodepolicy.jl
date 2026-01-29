# Domain: Node policies for fixed-node diffusion Monte Carlo

abstract type AbstractNodePolicy end

struct NoNode <: AbstractNodePolicy end

struct FixedNode <: AbstractNodePolicy end
