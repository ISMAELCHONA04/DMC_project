# UseCases/vmc: VMC energy estimators

"""
    estimate_energy_vmc(sim::VMCSim) -> Float64

Estimate the current mean local energy of the walker ensemble.
"""
function estimate_energy_vmc(sim::VMCSim)::Float64
    return estimate_energy(sim)
end

"""
    compute_variance(sim::VMCSim) -> Float64

Compute the variance of the current local-energy samples.
"""
function compute_variance(sim::VMCSim)::Float64
    return estimate_energy_variance(sim)
end
