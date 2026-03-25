# UseCases/gfmc: Runtime parameters for fixed-population Green Function Monte Carlo

"""
    GFMCParams(...)

Runtime parameters for fixed-population GFMC. See `docs/GFMC_USER_GUIDE.md`
for field definitions and the expected usage pattern.
"""
struct GFMCParams
    dt::Float64
    nsteps::Int
    nequil::Int
    targetN::Int
    ET0::Float64
    feedback::Float64
    reconfiguration_interval::Int
    branch_cap::Float64
    energy_window::Int
end

function GFMCParams(
    dt::Real,
    nsteps::Integer,
    nequil::Integer,
    targetN::Integer,
    ET0::Real,
    feedback::Real,
    reconfiguration_interval::Integer,
    branch_cap::Real,
    energy_window::Integer,
)
    dt_f = Float64(dt)
    dt_f > 0 || throw(ArgumentError("dt must be > 0, got $dt"))
    nsteps_int = Int(nsteps)
    nequil_int = Int(nequil)
    targetN_int = Int(targetN)
    reconf_int = Int(reconfiguration_interval)
    window_int = Int(energy_window)
    nsteps_int >= 1 || throw(ArgumentError("nsteps must be >= 1, got $nsteps"))
    0 <= nequil_int <= nsteps_int || throw(ArgumentError("nequil must satisfy 0 <= nequil <= nsteps"))
    targetN_int >= 1 || throw(ArgumentError("targetN must be >= 1, got $targetN"))
    reconf_int >= 1 || throw(ArgumentError("reconfiguration_interval must be >= 1, got $reconfiguration_interval"))
    window_int >= 1 || throw(ArgumentError("energy_window must be >= 1, got $energy_window"))
    branch_cap_f = Float64(branch_cap)
    branch_cap_f > 0 || throw(ArgumentError("branch_cap must be > 0, got $branch_cap"))
    return GFMCParams(
        dt_f,
        nsteps_int,
        nequil_int,
        targetN_int,
        Float64(ET0),
        Float64(feedback),
        reconf_int,
        branch_cap_f,
        window_int,
    )
end
