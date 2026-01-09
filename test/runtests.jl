using Test
using Random
using System1D


@testset "System1D basic" begin
    # Potential must return a scalar for a whole configuration (vector R)
    H = Hamiltonian(10, 0.5, R->sum(abs2, R)) # 
    W = Walker(rand(10))  # Distrute walkers along -1<x<1
    p = DMCParams(0.01, 5, 0, 10, 0.0, 0.1, 5, 5) # dt, nsteps, nequil, targetN, ET0, pop_gain, branch_cap, nblocks
    sim = DMCSim(H, p, [W], MersenneTwister(42), 0)

    @test isa(H, Hamiltonian)
    @test isa(W, Walker)
    @test isa(p, DMCParams)
    @test isa(sim, DMCSim)
    @test sim.guiding isa System1D.NoGuiding

    # run a short simulation and check diagnostics got recorded
    System1D.run_simulation!(sim)
    @test length(sim.population_history) == p.nsteps + 1
    @test length(sim.energy_mean_history) == p.nsteps + 1
    @test isa(System1D.estimate_energy(sim), Float64)

    trial = System1D.TrialWF(
        R -> -0.5 * sum(abs2, R),
        R -> -R,
        R -> -length(R)
    )
    guiding = System1D.ImportanceGuiding(trial, H)
    sim_guided = DMCSim(H, p, [Walker(rand(10))], guiding, MersenneTwister(24), 0)

    System1D.run_simulation!(sim_guided)
    @test length(sim_guided.population_history) == p.nsteps + 1
    @test length(sim_guided.energy_mean_history) == p.nsteps + 1
    @test isa(System1D.estimate_energy(sim_guided), Float64)

end
