# UseCases/vmc: Variational Monte Carlo move proposals with importance sampling

"""
    proposal_step!(sim::VMCSim, w::Walker) -> (Rnew, log_ratio)

Generate a drift-diffusion proposal move for VMC with importance sampling.
Returns the new position and the log acceptance ratio.
"""
function proposal_step!(sim::VMCSim, w::Walker)
    D  = diffusion_constant(sim.H)
    dt = sim.params.dt
    Rold  = position(w)
    
    # Compute drift from trial wavefunction
    Fold = drift(sim.guiding, Rold)
    
    # Forward proposal: drift + diffusion
    Rnew = Rold .+ D*dt .* Fold .+ sqrt(2 * D * dt) .* randn(sim.rng, length(Rold))
    
    # Backward drift for detailed balance
    Fnew = drift(sim.guiding, Rnew)
    denom = 4 * D * dt
    
    # Green's function ratio for detailed balance
    log_gf = -sum(abs2, Rnew .- Rold .- D*dt .* Fold) / denom
    log_gb = -sum(abs2, Rold .- Rnew .- D*dt .* Fnew) / denom
    
    logpsi_old = sim.trial_wf.logpsi(Rold)
    logpsi_new = sim.trial_wf.logpsi(Rnew)
    
    # Log acceptance ratio
    log_ratio = log_gb - log_gf + 2 * (logpsi_new - logpsi_old)
    
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
        return Walker(copy(Rnew))
    else
        return w
    end
end
