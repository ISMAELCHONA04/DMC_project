# UseCases/dmc: Move proposal and node crossing detection

"""
    crosses_node(nodepolicy, guiding, Rold, Rnew) -> Bool

Return `true` if a proposed move crosses the nodal surface defined by the trial
wavefunction sign.
"""
crosses_node(::NoNode, guiding, Rold, Rnew) = false
function crosses_node(::FixedNode, guiding, Rold, Rnew)
    sold = signpsi(guiding, Rold)
    snew = signpsi(guiding, Rnew)
    return (sold == 0) || (snew == 0) || (sold != snew)
end

# Propose a diffusion move without importance sampling.
function propose_move(sim::DMCSim, w::Walker, ::NoGuiding)
    Dcoef = diffusion_constant(sim.H)   # from your Domain accessor
    dt = sim.params.dt
    Rold = position(w)                 # from your Domain accessor
    # Sample new position from diffusion PDF
    Rnew = Rold .+ sqrt(2 * Dcoef * dt) .* randn(sim.rng, length(Rold))
    return Rnew
end

# Drift-diffusion proposal with Metropolis accept/reject.
function propose_move(sim::DMCSim, w::Walker, g::ImportanceGuiding)
    D  = diffusion_constant(sim.H)
    dt = sim.params.dt
    Rold  = position(w)
    Fold = drift(g, Rold)
    # Forward proposal: drift + diffusion.
    Rnew = Rold .+ D*dt .* Fold .+ sqrt(2 * D * dt) .* randn(sim.rng, length(Rold))

    Fnew = drift(g, Rnew)
    denom = 4 * D * dt
    # Green's function ratio for detailed balance.
    log_gf = -sum(abs2, Rnew .- Rold .- D*dt .* Fold) / denom
    log_gb = -sum(abs2, Rold .- Rnew .- D*dt .* Fnew) / denom
    logpsi_old = g.trial.logpsi(Rold)
    logpsi_new = g.trial.logpsi(Rnew)
    log_ratio = log_gb - log_gf + 2 * (logpsi_new - logpsi_old)

    # Accept/reject step.
    if log(rand(sim.rng)) < min(0.0, log_ratio)
        return Rnew
    end
    return Rold
end
