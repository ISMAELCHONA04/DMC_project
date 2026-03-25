using Test
using Random
using System1D


@testset "System1D basic" begin
    # Potential must return a scalar for a whole configuration (vector R)
    H = Hamiltonian(10, 0.5, R->sum(abs2, R)) # 
    W = Walker(rand(10))  # Distrute walkers along -1<x<1
    p = DMCParams(0.01, 5, 0, 10, 0.0, 1.0, 5, 5) # dt, nsteps, nequil, targetN, ET0, population_control_gain, branch_cap, nblocks
    sim = DMCSim(H, p, [W], MersenneTwister(42), 0)

    @test isa(H, Hamiltonian)
    @test isa(W, Walker)
    @test isa(p, DMCParams)
    @test isa(sim, DMCSim)
    @test sim.guiding isa System1D.NoGuiding
    @test p.population_control_gain == 1.0

    p_legacy = DMCParams(0.01, 5, 0, 10, 0.0, 5, 5)
    @test p_legacy == p

    p_kw = DMCParams(; dt=0.01, nsteps=5, nequil=0, targetN=10, ET0=0.0, population_control_gain=1.0, branch_cap=5, nblocks=5)
    @test p_kw == p

    p_alias = DMCParams(; dt=0.01, nsteps=5, nequil=0, targetN=10, ET0=0.0, pop_gain=1.0, branch_cap=5, nblocks=5)
    @test p_alias == p

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

@testset "Periodic boundary behavior" begin
    bc = System1D.PeriodicBoundary1D(0.0, 1.0)

    @test System1D.isperiodic(bc)
    @test !System1D.isperiodic(System1D.OpenBoundary1D())
    @test System1D.cell_length(bc) ≈ 1.0
    @test System1D.cell_bounds(bc) == (0.0, 1.0)

    @test System1D.wrap_coordinate(bc, 1.2) ≈ 0.2
    @test System1D.wrap_coordinate(bc, -0.1) ≈ 0.9
    @test System1D.minimum_image(bc, 0.6) ≈ -0.4
    @test System1D.minimum_image(bc, -0.6) ≈ 0.4
    @test System1D.displacement(bc, 0.95, 0.05) ≈ 0.1
    @test System1D.distance_1d(bc, 0.95, 0.05) ≈ 0.1
    @test System1D.squared_distance_1d(bc, 0.95, 0.05) ≈ 0.01

    bc_twist = System1D.PeriodicBoundary1D(0.0, 1.0; twist=pi / 2)
    wrapped, phase_factor = System1D.wrap_position_with_phase(bc_twist, [1.25])
    @test wrapped[1] ≈ 0.25
    @test phase_factor ≈ cis(pi / 2)

    H = System1D.Hamiltonian(1, 0.5, R -> 0.0, bc)
    p = System1D.DMCParams(0.2, 2, 0, 1, 0.0, 1.0, 10, 1)
    w = System1D.Walker([0.99])
    sim = System1D.DMCSim(H, p, [w], MersenneTwister(11), 0)

    # No-guiding proposals should always wrap into the periodic cell.
    for _ in 1:200
        Rnew = System1D.propose_move(sim, w, System1D.NoGuiding())
        @test 0.0 <= Rnew[1] < 1.0
    end

    # Importance-sampled proposals should also stay in-cell.
    trial = System1D.TrialWF(
        R -> 0.0,
        R -> zeros(length(R)),
        R -> 0.0
    )
    guiding = System1D.ImportanceGuiding(trial, H)
    for _ in 1:200
        Rnew = System1D.propose_move(sim, w, guiding)
        @test 0.0 <= Rnew[1] < 1.0
    end

    H_twist = System1D.Hamiltonian(1, 0.5, R -> 0.0, bc_twist)
    sim_twist = System1D.DMCSim(H_twist, p, [[1.25]], MersenneTwister(12))
    @test System1D.position(sim_twist.walkers[1])[1] ≈ 0.25
    @test System1D.phase(sim_twist.walkers[1]) ≈ cis(pi / 2)
end

@testset "DMC plotting helpers" begin
    snap = [[-1.0], [0.0], [1.0], [2.0]]
    @test System1D.mean_walker_position_1d(snap) ≈ 0.5

    vals = [1.0, 3.0, 5.0, 7.0]
    avg2 = System1D.sliding_window_average(vals, 2)
    @test avg2 ≈ [1.0, 2.0, 4.0, 6.0]

    avg3 = System1D.sliding_window_average(vals, 3)
    @test avg3 ≈ [1.0, 2.0, 3.0, 5.0]

    @test_throws ArgumentError System1D.sliding_window_average(vals, 0)

    p = System1D.plot_snapshot_1d_smoothed_density(
        snap;
        nbins=12,
        smoothing_window=3,
        xmin=-2.0,
        xmax=3.0
    )
    @test p !== nothing
end

@testset "VMC proposal modes" begin
    bc = System1D.PeriodicBoundary1D(0.0, 1.0)
    H = System1D.Hamiltonian(1, 0.5, R -> 0.0, bc)
    trial = System1D.TrialWF(
        R -> -0.5 * sum(abs2, R),
        R -> -R,
        R -> -length(R)
    )
    params = System1D.VMCParams(; dt=0.02, nsteps=3, targetN=8, ET0=0.0)
    @test params == System1D.VMCParams(0.02, 3, 8, 0.0)
    walkers = [[0.97] for _ in 1:params.targetN]

    sim_default = System1D.VMCSim(H, params, walkers, trial, MersenneTwister(17))
    @test sim_default.proposal isa System1D.DriftGaussianProposal
    @test all(0.0 <= System1D.position(w)[1] < 1.0 for w in sim_default.walkers)

    for _ in 1:100
        Rnew, log_ratio = System1D.proposal_step!(sim_default, sim_default.walkers[1])
        @test 0.0 <= Rnew[1] < 1.0
        @test isfinite(log_ratio)
    end

    sim_gauss_sym = System1D.VMCSim(H, params, walkers, trial, MersenneTwister(18); proposal=:gaussian)
    @test sim_gauss_sym.proposal isa System1D.GaussianProposal
    Rnew, log_ratio = System1D.proposal_step!(sim_gauss_sym, sim_gauss_sym.walkers[1])
    @test 0.0 <= Rnew[1] < 1.0
    @test isfinite(log_ratio)

    sim_gauss_obj = System1D.VMCSim(H, params, walkers, trial, MersenneTwister(19); proposal=System1D.GaussianProposal())
    @test sim_gauss_obj.proposal isa System1D.GaussianProposal
    System1D.run_vmc!(sim_gauss_obj; snapshot_steps=[0, params.nsteps])
    @test sim_gauss_obj.step == params.nsteps
    @test length(sim_gauss_obj.energy_history) == params.nsteps + 1
    @test length(sim_gauss_obj.acceptance_history) == params.nsteps + 1
    @test length(sim_gauss_obj.walker_positions_history) == 2

    sim_run = System1D.run_vmc(
        H,
        params,
        walkers,
        trial;
        rng=MersenneTwister(21),
        proposal=:gaussian,
        snapshot_steps=[0, params.nsteps]
    )
    @test sim_run isa System1D.VMCSim
    @test sim_run.step == params.nsteps
    @test length(sim_run.acceptance_history) == params.nsteps + 1
    @test length(sim_run.walker_positions_history) == 2

    @test_throws ArgumentError System1D.VMCSim(H, params, walkers, trial, MersenneTwister(20); proposal=:bad_mode)
end

@testset "VMC-backed GFMC warm start" begin
    bc = System1D.PeriodicBoundary1D(0.0, 1.0)
    H = System1D.Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R), bc)
    trial = System1D.TrialWF(
        R -> -0.5 * sum(abs2, R),
        R -> -R,
        R -> -length(R)
    )

    vmc_params = System1D.VMCParams(; dt=0.02, nsteps=3, targetN=6, ET0=0.0)
    gfmc_params = System1D.GFMCParams(0.02, 3, 0, 6, 0.5, 0.2, 1, 5.0, 3)
    seed_positions = [[0.97] for _ in 1:vmc_params.targetN]

    warm_positions = System1D.sample_vmc_positions(
        H,
        vmc_params,
        seed_positions,
        trial;
        rng=MersenneTwister(41),
        proposal=:gaussian,
    )
    @test length(warm_positions) == vmc_params.targetN
    @test all(length(R) == 1 for R in warm_positions)
    @test all(0.0 <= R[1] < 1.0 for R in warm_positions)

    sim = System1D.run_gfmc_with_vmc_init(
        H,
        gfmc_params,
        seed_positions,
        trial,
        vmc_params;
        vmc_rng=MersenneTwister(42),
        gfmc_rng=MersenneTwister(43),
        proposal=:gaussian,
        snapshot_steps=[0, gfmc_params.nsteps],
    )
    @test sim isa System1D.GFMCSim
    @test sim.step == gfmc_params.nsteps
    @test length(sim.walker_positions_history) == 2
    @test all(==(gfmc_params.targetN), sim.population_history)
    @test sim.kernel.guiding isa System1D.ImportanceGuiding

    @test_throws ArgumentError System1D.sample_vmc_positions(
        H,
        vmc_params,
        seed_positions[1:end-1],
        trial;
        rng=MersenneTwister(44),
    )

    bad_vmc_params = System1D.VMCParams(; dt=0.02, nsteps=3, targetN=5, ET0=0.0)
    @test_throws ArgumentError System1D.run_gfmc_with_vmc_init(
        H,
        gfmc_params,
        seed_positions[1:end-1],
        trial,
        bad_vmc_params;
        vmc_rng=MersenneTwister(45),
        gfmc_rng=MersenneTwister(46),
    )

    model = System1D.SpinlessFermionRing1D(2, 1.0, 2.0, 1.0; D=0.5, node_tol=1.0e-7, trig_eps=1.0e-10)
    vmc_model_params = System1D.VMCParams(; dt=1.0e-3, nsteps=2, targetN=8, ET0=0.0)
    gfmc_model_params = System1D.GFMCParams(1.0e-3, 2, 0, 8, -0.1, 0.1, 1, 5.0, 2)
    model_seed_positions = System1D.sample_uniform_configurations(model, gfmc_model_params.targetN, MersenneTwister(47))

    model_sim = System1D.run_gfmc_with_vmc_init(
        model,
        gfmc_model_params,
        model_seed_positions,
        vmc_model_params;
        vmc_rng=MersenneTwister(48),
        gfmc_rng=MersenneTwister(49),
        proposal=:gaussian,
        snapshot_steps=[0, gfmc_model_params.nsteps],
    )
    @test model_sim isa System1D.GFMCSim
    @test model_sim.step == gfmc_model_params.nsteps
    @test length(model_sim.walker_positions_history) == 2
    @test all(==(gfmc_model_params.targetN), model_sim.population_history)
end

@testset "DMC convenience runner" begin
    H = System1D.Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R))
    params = System1D.DMCParams(; dt=0.01, nsteps=3, nequil=0, targetN=2, ET0=0.5, population_control_gain=1.0, branch_cap=4, nblocks=2)
    positions = [[0.1], [0.2]]

    sim = System1D.run_dmc(H, params, positions; rng=MersenneTwister(99), snapshot_steps=[0, params.nsteps])

    @test sim isa System1D.DMCSim
    @test sim.step == params.nsteps
    @test length(sim.walker_positions_history) == 2
end

@testset "Debug runner output" begin
    dmc_H = System1D.Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R))
    dmc_params = System1D.DMCParams(; dt=0.01, nsteps=2, nequil=0, targetN=2, ET0=0.5, population_control_gain=1.0, branch_cap=4, nblocks=2)
    dmc_buf = IOBuffer()
    System1D.run_dmc(
        dmc_H,
        dmc_params,
        [[0.1], [0.2]];
        rng=MersenneTwister(301),
        debug=true,
        debug_io=dmc_buf,
    )
    dmc_log = String(take!(dmc_buf))
    @test occursin("[debug]", dmc_log)
    @test occursin("population=", dmc_log)
    @test occursin("ET=", dmc_log)

    vmc_H = System1D.Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R))
    vmc_trial = System1D.TrialWF(
        R -> -0.5 * sum(abs2, R),
        R -> -R,
        R -> -length(R)
    )
    vmc_params = System1D.VMCParams(; dt=0.02, nsteps=2, targetN=4, ET0=0.5)
    vmc_buf = IOBuffer()
    System1D.run_vmc(
        vmc_H,
        vmc_params,
        [[0.1], [0.2], [0.3], [0.4]],
        vmc_trial;
        rng=MersenneTwister(302),
        proposal=:gaussian,
        debug=true,
        debug_io=vmc_buf,
    )
    vmc_log = String(take!(vmc_buf))
    @test occursin("[debug]", vmc_log)
    @test occursin("acceptance_rate=", vmc_log)
    @test occursin("energy_mean=", vmc_log)

    gfmc_bc = System1D.PeriodicBoundary1D(0.0, 1.0)
    gfmc_H = System1D.Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R), gfmc_bc)
    gfmc_trial = System1D.TrialWF(
        R -> -0.5 * sum(abs2, R),
        R -> -R,
        R -> -length(R)
    )
    gfmc_guiding = System1D.ImportanceGuiding(gfmc_trial, gfmc_H)
    gfmc_params = System1D.GFMCParams(0.02, 2, 0, 4, 0.5, 0.2, 1, 5.0, 2)
    gfmc_walkers = [System1D.Walker([0.1 * i]) for i in 1:gfmc_params.targetN]
    gfmc_sim = System1D.GFMCSim(gfmc_H, gfmc_params, gfmc_walkers, MersenneTwister(303); guiding=gfmc_guiding, nodepolicy=System1D.NoNode())
    gfmc_buf = IOBuffer()
    System1D.run_gfmc!(
        gfmc_sim;
        snapshot_steps=[0, gfmc_params.nsteps],
        debug=true,
        debug_io=gfmc_buf,
    )
    gfmc_log = String(take!(gfmc_buf))
    @test occursin("[debug]", gfmc_log)
    @test occursin("mean_weight=", gfmc_log)
    @test occursin("effective_population=", gfmc_log)
end

@testset "GFMC fixed walker population" begin
    bc = System1D.PeriodicBoundary1D(0.0, 1.0)
    H = System1D.Hamiltonian(1, 0.5, R -> 0.5 * sum(abs2, R), bc)
    trial = System1D.TrialWF(
        R -> -0.5 * sum(abs2, R),
        R -> -R,
        R -> -length(R)
    )
    guiding = System1D.ImportanceGuiding(trial, H)

    params = System1D.GFMCParams(0.02, 4, 0, 16, 0.5, 0.2, 1, 5.0, 4)
    walkers = [System1D.Walker([rand(MersenneTwister(44 + i))]) for i in 1:params.targetN]
    sim = System1D.GFMCSim(H, params, walkers, MersenneTwister(91); guiding=guiding, nodepolicy=System1D.NoNode())
    raw_positions = [[0.01 * i] for i in 1:params.targetN]
    sim_from_positions = System1D.GFMCSim(H, params, raw_positions, MersenneTwister(92); guiding=guiding, nodepolicy=System1D.NoNode())

    System1D.run_gfmc!(sim; snapshot_steps=[0, params.nsteps])
    System1D.run_gfmc!(sim_from_positions; snapshot_steps=[0, params.nsteps])

    @test sim.step == params.nsteps
    @test sim_from_positions.step == params.nsteps
    @test length(sim.population_history) == params.nsteps + 1
    @test all(==(params.targetN), sim.population_history)
    @test all(==(params.targetN), sim_from_positions.population_history)
    @test length(sim.energy_mean_history) == params.nsteps + 1
    @test length(sim.energy_variance_history) == params.nsteps + 1
    @test length(sim.mean_weight_history) == params.nsteps + 1
    @test length(sim.effective_population_history) == params.nsteps + 1
    @test length(sim.acceptance_history) == params.nsteps + 1
    @test length(sim.walker_positions_history) == 2
    @test all(isfinite, sim.ET_history)
end

@testset "GFMC spinless fermion ring kernel" begin
    model = System1D.SpinlessFermionRing1D(3, 1.0, 3.0, 1.0; lambda=-0.5, node_tol=1e-8, trig_eps=1e-10)
    params = System1D.GFMCParams(0.003, 3, 0, 12, -0.2, 0.1, 1, 5.0, 3)
    positions = System1D.sample_uniform_configurations(model, params.targetN, MersenneTwister(101))
    sim = System1D.GFMCSim(model, params, positions, MersenneTwister(202); use_guiding=true, nodepolicy=System1D.FixedNode())

    System1D.run_gfmc!(sim; snapshot_steps=[0, params.nsteps])

    @test sim.kernel isa System1D.SpinlessFermionRing1DKernel
    @test sim.step == params.nsteps
    @test all(==(params.targetN), sim.population_history)
    @test length(sim.walker_positions_history) == 2
    @test all(isfinite, sim.energy_mean_history)
    @test all(isfinite, sim.energy_variance_history)
    @test all(isfinite, sim.mean_weight_history)
    @test all(isfinite, sim.effective_population_history)
    @test all(0.0 <= x < model.L for snap in sim.walker_positions_history for R in snap for x in R)

    safe_R = [0.15, 1.05, 2.15]
    @test isfinite(System1D.trial_logpsi(model, safe_R))
    @test all(isfinite, System1D.trial_gradlogpsi(model, safe_R))
    @test isfinite(System1D.trial_lapllogpsi(model, safe_R))
    @test System1D.signpsi(model, safe_R) in (-1.0, 1.0)
    @test isfinite(System1D.local_energy(model, safe_R))
end
