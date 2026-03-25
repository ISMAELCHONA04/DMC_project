# UseCases/dmc: Runtime parameters for DMC

"""
    DMCParams(...)

Runtime parameters for diffusion Monte Carlo. See `docs/DMC_USER_GUIDE.md` for
the constructor signatures and parameter semantics.
"""
struct DMCParams
    dt::Float64
    nsteps::Int
    nequil::Int
    targetN::Int     # target walker population
    ET0::Float64     # initial reference-energy guess
    population_control_gain::Float64  # strength of population-control feedback in the ET update
    branch_cap::Int  # max branching factor
    nblocks::Int     # averaging window used by the ET update
end

function DMCParams(
    dt::Real,
    nsteps::Integer,
    nequil::Integer,
    targetN::Integer,
    ET0::Real,
    population_control_gain::Real,
    branch_cap::Integer,
    nblocks::Integer,
)
    dt_f = Float64(dt)
    nsteps_i = Int(nsteps)
    nequil_i = Int(nequil)
    targetN_i = Int(targetN)
    population_control_gain_f = Float64(population_control_gain)
    branch_cap_i = Int(branch_cap)
    nblocks_i = Int(nblocks)

    dt_f > 0 || throw(ArgumentError("dt must be > 0, got $dt"))
    nsteps_i >= 1 || throw(ArgumentError("nsteps must be >= 1, got $nsteps"))
    0 <= nequil_i <= nsteps_i || throw(ArgumentError("nequil must satisfy 0 <= nequil <= nsteps; got nequil=$nequil, nsteps=$nsteps"))
    targetN_i >= 1 || throw(ArgumentError("targetN must be >= 1, got $targetN"))
    isfinite(population_control_gain_f) || throw(ArgumentError("population_control_gain must be finite, got $population_control_gain"))
    population_control_gain_f >= 0 || throw(ArgumentError("population_control_gain must be >= 0, got $population_control_gain"))
    branch_cap_i >= 1 || throw(ArgumentError("branch_cap must be >= 1, got $branch_cap"))
    nblocks_i >= 1 || throw(ArgumentError("nblocks must be >= 1, got $nblocks"))

    return DMCParams(dt_f, nsteps_i, nequil_i, targetN_i, Float64(ET0), population_control_gain_f, branch_cap_i, nblocks_i)
end

function DMCParams(
    dt::Real,
    nsteps::Integer,
    nequil::Integer,
    targetN::Integer,
    ET0::Real,
    branch_cap::Integer,
    nblocks::Integer,
)
    return DMCParams(dt, nsteps, nequil, targetN, ET0, 1.0, branch_cap, nblocks)
end

function DMCParams(;
    dt::Real,
    nsteps::Integer,
    targetN::Integer,
    nequil::Integer=0,
    ET0::Real=0.0,
    branch_cap::Integer=8,
    nblocks::Integer=50,
    population_control_gain::Real=1.0,
    pop_gain::Union{Nothing,Real}=nothing,
)
    population_control_gain_f = Float64(population_control_gain)
    if pop_gain !== nothing
        pop_gain_f = Float64(pop_gain)
        if population_control_gain_f != 1.0 && population_control_gain_f != pop_gain_f
            throw(ArgumentError("population_control_gain and pop_gain disagree; pass only one name or give them the same value"))
        end
        population_control_gain_f = pop_gain_f
    end

    return DMCParams(dt, nsteps, nequil, targetN, ET0, population_control_gain_f, branch_cap, nblocks)
end
