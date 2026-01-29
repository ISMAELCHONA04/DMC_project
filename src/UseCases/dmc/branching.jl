# UseCases/dmc: Branching factor calculation for walker replication

# Branching factor for no importance sampling (uses potential energy).
function branching_factor(sim::DMCSim, Rold, Rnew, ::NoGuiding)
    dt = sim.params.dt
    ET = sim.ET
    Vold = potential(sim.H, Rold)
    Vnew = potential(sim.H, Rnew)
    P = exp(-0.5*dt * ((Vold + Vnew) - 2ET))
    return min(P, sim.params.branch_cap)
end

# Branching factor for importance sampling (uses local energy).
function branching_factor(sim::DMCSim, Rold, Rnew, g::ImportanceGuiding)
    dt = sim.params.dt
    ET = sim.ET
    ELold = local_energy(g, Rold)
    ELnew = local_energy(g, Rnew)
    P = exp(-0.5*dt * ((ELold + ELnew) - 2ET))
    return min(P, sim.params.branch_cap)
end
