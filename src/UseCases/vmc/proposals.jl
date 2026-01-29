# UseCases/vmc: Variational Monte Carlo move proposals (skeleton)

"""
    metropolis_step!(sim::VMCSim, w::Walker) -> Walker

Perform a Metropolis step for a single walker (skeleton implementation).
Currently not implemented.
"""
function metropolis_step!(sim::VMCSim, w::Walker)
    error("VMC not implemented yet")
end

"""
    proposal_step!(sim::VMCSim, w::Walker) -> Walker

Generate a proposal move for variational Monte Carlo (skeleton implementation).
Currently not implemented.
"""
function proposal_step!(sim::VMCSim, w::Walker)
    error("VMC not implemented yet")
end
