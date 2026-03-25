# UseCases/vmc: Variational Monte Carlo move proposals with importance sampling

"""
    proposal_step!(sim::VMCSim, w::Walker) -> (Rnew, log_ratio)

Generate a VMC proposal move based on `sim.proposal`.
Returns the new position and the log acceptance ratio.
"""
function proposal_step!(sim::VMCSim, w::Walker)
    return proposal_step!(sim, w, sim.proposal)
end

"""
    proposal_step!(sim::VMCSim, w::Walker, ::DriftGaussianProposal) -> (Rnew, log_ratio)

Drift-diffusion proposal (default): Gaussian centered at `Rold + D*dt*F`.
"""
function proposal_step!(sim::VMCSim, w::Walker, ::DriftGaussianProposal)
    D  = diffusion_constant(sim.H)
    bc = boundary(sim.H)
    dt = sim.params.dt
    Rold  = position(w)
    
    # Compute drift from trial wavefunction
    Fold = drift(sim.guiding, Rold)
    
    # Forward proposal: drift + diffusion
    sigma = sqrt(2 * D * dt)
    mu_f = Rold .+ D*dt .* Fold
    Rnew = wrap_position(bc, mu_f .+ sigma .* randn(sim.rng, length(Rold)))
    
    # Backward drift for detailed balance
    Fnew = drift(sim.guiding, Rnew)
    denom = 4 * D * dt
    mu_b = Rnew .+ D*dt .* Fnew
    
    # Green's function ratio for detailed balance
    # Use minimum-image displacements so periodic moves are measured correctly.
    df = displacement(bc, mu_f, Rnew)
    db = displacement(bc, mu_b, Rold)
    log_gf = -sum(abs2, df) / denom
    log_gb = -sum(abs2, db) / denom
    
    logpsi_old = sim.trial_wf.logpsi(Rold)
    logpsi_new = sim.trial_wf.logpsi(Rnew)
    
    # Log acceptance ratio
    log_ratio = log_gb - log_gf + 2 * (logpsi_new - logpsi_old)
    
    return Rnew, log_ratio
end

"""
    proposal_step!(sim::VMCSim, w::Walker, ::GaussianProposal) -> (Rnew, log_ratio)

Pure diffusion proposal: Gaussian centered at current position.
"""
function proposal_step!(sim::VMCSim, w::Walker, ::GaussianProposal)
    D  = diffusion_constant(sim.H)
    bc = boundary(sim.H)
    dt = sim.params.dt
    Rold = position(w)

    sigma = sqrt(2 * D * dt)
    Rnew = wrap_position(bc, Rold .+ sigma .* randn(sim.rng, length(Rold)))

    # Symmetric proposal kernel q(Rold->Rnew) == q(Rnew->Rold)
    logpsi_old = sim.trial_wf.logpsi(Rold)
    logpsi_new = sim.trial_wf.logpsi(Rnew)
    log_ratio = 2 * (logpsi_new - logpsi_old)
    return Rnew, log_ratio
end

"""
    metropolis_step!(sim::VMCSim, w::Walker) -> Walker

Perform a Metropolis accept/reject step for a single walker.
Returns accepted or rejected walker position.
"""
function metropolis_step!(sim::VMCSim, w::Walker)
    Rnew, log_ratio = proposal_step!(sim, w)
    
    # Metropolis accept/reject
    if log(rand(sim.rng)) < min(0.0, log_ratio)
        sim.acceptance_count += 1
        return Walker(copy(Rnew), phase(w))
    else
        return w
    end
end
