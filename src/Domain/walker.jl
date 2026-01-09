# Domain: Walker (particle state)

mutable struct Walker{X}
    R::X
end

# Constructor and accessor methods
Walker(R::X) where {X} = Walker{X}(R)
position(w::Walker) = w.R
clone(w::Walker) = Walker(copy(w.R)) # copy of the walker with same position
