# UseCases/dmc: Move proposal and node crossing detection

struct DMCMove{R}
    position::R
    phase_factor::ComplexF64
end

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
function _propose_move(sim::DMCSim, w::Walker, ::NoGuiding)
    Dcoef = diffusion_constant(sim.H)   # from your Domain accessor
    dt = sim.params.dt
    Rold = position(w)                 # from your Domain accessor
    # Sample new position from diffusion PDF
    Rtrial = Rold .+ sqrt(2 * Dcoef * dt) .* randn(sim.rng, length(Rold))
    Rnew, phase_factor = wrap_position_with_phase(boundary(sim.H), Rtrial)
    return DMCMove(Rnew, phase_factor)
end

# Drift-diffusion proposal with Metropolis accept/reject.
function _propose_move(sim::DMCSim, w::Walker, g::ImportanceGuiding)
    D  = diffusion_constant(sim.H)
    dt = sim.params.dt
    Rold = position(w)
    Fold = drift(g, Rold)
    bc = boundary(sim.H)
    sigma = sqrt(2 * D * dt)

    # Forward proposal: drift + diffusion, then wrap into the simulation cell.
    mu_f = Rold .+ D * dt .* Fold
    Rnew, phase_factor = wrap_position_with_phase(bc, mu_f .+ sigma .* randn(sim.rng, length(Rold)))

    Fnew = drift(g, Rnew)
    denom = 4 * D * dt
    mu_b = Rnew .+ D * dt .* Fnew

    # Green's function ratio for detailed balance.
    # Use minimum-image displacements so periodic moves are measured correctly.
    df = displacement(bc, mu_f, Rnew)
    db = displacement(bc, mu_b, Rold)
    log_gf = -sum(abs2, df) / denom
    log_gb = -sum(abs2, db) / denom
    logpsi_old = g.trial.logpsi(Rold)
    logpsi_new = g.trial.logpsi(Rnew)
    log_ratio = log_gb - log_gf + 2 * (logpsi_new - logpsi_old)

    # Accept/reject step.
    if log(rand(sim.rng)) < min(0.0, log_ratio)
        return DMCMove(Rnew, phase_factor)
    end
    return DMCMove(Rold, 1.0 + 0.0im)
end

function propose_move(sim::DMCSim, w::Walker, guiding::AbstractGuiding)
    return _propose_move(sim, w, guiding).position
end
