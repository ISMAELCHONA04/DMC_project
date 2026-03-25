# UseCases/vmc: Runtime parameters for VMC

"""
    VMCParams(...)

Runtime parameters for variational Monte Carlo. See `docs/VMC_USER_GUIDE.md`
for the public constructor signatures.
"""
struct VMCParams
    dt::Float64
    nsteps::Int
    targetN::Int
    ET0::Float64
end

function VMCParams(
    dt::Real,
    nsteps::Integer,
    targetN::Integer,
    ET0::Real,
)
    dt_f = Float64(dt)
    nsteps_i = Int(nsteps)
    targetN_i = Int(targetN)

    dt_f > 0 || throw(ArgumentError("dt must be > 0, got $dt"))
    nsteps_i >= 1 || throw(ArgumentError("nsteps must be >= 1, got $nsteps"))
    targetN_i >= 1 || throw(ArgumentError("targetN must be >= 1, got $targetN"))

    return VMCParams(dt_f, nsteps_i, targetN_i, Float64(ET0))
end

function VMCParams(;
    dt::Real,
    nsteps::Integer,
    targetN::Integer,
    ET0::Real=0.0,
)
    return VMCParams(dt, nsteps, targetN, ET0)
end
