# UseCases/vmc: Main VMC simulation loop (skeleton)

"""
    vmc_step!(sim::VMCSim)

Perform a single VMC step (skeleton implementation).
Currently not implemented.
"""
function vmc_step!(sim::VMCSim)
    error("VMC not implemented yet")
end

"""
    run_vmc!(sim::VMCSim; snapshot_steps=Int[])

Run the VMC simulation (skeleton implementation).
Currently not implemented.

# Arguments
- `sim::VMCSim`: The VMC simulation object
- `snapshot_steps::AbstractVector{<:Integer}`: Steps at which to record walker positions

# Raises
- `error("VMC not implemented yet")`
"""
function run_vmc!(sim::VMCSim; snapshot_steps::AbstractVector{<:Integer}=Int[])
    error("VMC not implemented yet")
end
