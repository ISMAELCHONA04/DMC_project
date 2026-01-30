# UseCases/vmc: VMC energy estimators (skeleton)

"""
    estimate_energy_vmc(sim::VMCSim) -> Float64

Estimate energy from current walker ensemble in VMC (skeleton implementation).
Currently not implemented.
"""
function estimate_energy_vmc(sim::VMCSim)::Float64
    if isempty(sim.walkers)
        return 0.0
    end
    vals = Float64[local_energy(sim.guiding, position(w)) for w in sim.walkers]
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end

"""
    compute_variance(sim::VMCSim) -> Float64

Compute energy variance for the walker ensemble (skeleton implementation).
Currently not implemented.
"""
function compute_variance(sim::VMCSim)::Float64
    if isempty(sim.walkers)
        return 0.0
    end
    energies = Float64[local_energy(sim.guiding, position(w)) for w in sim.walkers]
    mean_E = sum(energies) / length(energies)
    variance = sum((E - mean_E)^2 for E in energies) / length(energies)
    return variance
end
