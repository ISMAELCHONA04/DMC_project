# Domain: Walker (particle state)

"""
    Walker(R)
    Walker(R, phase)

Mutable walker state containing a configuration vector `R` and a complex phase
accumulated from twisted periodic-boundary crossings.
"""
mutable struct Walker{X}
    R::X
    phase::ComplexF64
end

# Constructor and accessor methods
Walker(R::X) where {X} = Walker{X}(R, 1.0 + 0.0im)
Walker(R::X, phase::Number) where {X} = Walker{X}(R, ComplexF64(phase))
"""Return the configuration vector stored in walker `w`."""
position(w::Walker) = w.R
"""Return the accumulated complex phase stored in walker `w`."""
phase(w::Walker) = w.phase
"""Return a deep copy of walker `w`, including its phase."""
clone(w::Walker) = Walker(copy(w.R), w.phase)
