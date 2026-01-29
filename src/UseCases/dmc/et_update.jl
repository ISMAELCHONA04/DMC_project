# UseCases/dmc: Reference energy update based on population growth

"""
    update_ET!(sim::DMCSim)

Update reference energy based on current population (growth estimator).
Uses the average reference energy over the last `nblocks` steps.
"""
function update_ET!(sim::DMCSim)
    Nt = length(sim.walkers)
    N0 = sim.params.targetN
    dt = sim.params.dt
    window = sim.params.nblocks

    if isempty(sim.ET_history)
        Eblock = sim.ET
    else
        hi = lastindex(sim.ET_history)
        lo = max(1, hi - window + 1)
        block = @view sim.ET_history[lo:hi]
        Eblock = sum(block) / length(block)
    end

    sim.ET = Eblock - (1 / dt) * log(Nt / N0)
    return sim
end
