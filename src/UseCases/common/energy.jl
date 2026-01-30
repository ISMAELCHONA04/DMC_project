# UseCases/common: Energy estimation utilities shared by DMC and VMC

"""
    estimate_energy(sim) -> Float64

Estimate mean energy from current walker set using the appropriate guiding method.
"""
function estimate_energy(sim)::Float64
    return estimate_energy(sim, sim.guiding)
end

function estimate_energy(sim, ::NoGuiding)::Float64
    vals = Float64[potential(sim.H, position(w)) for w in sim.walkers]
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end

function estimate_energy(sim, g::ImportanceGuiding)::Float64
    vals = Float64[local_energy(g, position(w)) for w in sim.walkers]
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end

"""
    estimate_energy_variance(sim) -> Float64

Compute variance of local energy from current walker set using the appropriate guiding method.
"""
function estimate_energy_variance(sim)::Float64
    return estimate_energy_variance(sim, sim.guiding)
end

function estimate_energy_variance(sim, ::NoGuiding)::Float64
    vals = Float64[potential(sim.H, position(w)) for w in sim.walkers]
    if isempty(vals)
        return 0.0
    end
    mean_E = sum(vals) / length(vals)
    return sum((E - mean_E)^2 for E in vals) / length(vals)
end

function estimate_energy_variance(sim, g::ImportanceGuiding)::Float64
    vals = Float64[local_energy(g, position(w)) for w in sim.walkers]
    if isempty(vals)
        return 0.0
    end
    mean_E = sum(vals) / length(vals)
    return sum((E - mean_E)^2 for E in vals) / length(vals)
end
