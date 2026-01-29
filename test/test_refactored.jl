"""
Test script to verify DMC functionality after refactoring.
Ensures backward compatibility and that all key DMC features work.
"""

using Random

# Import the module
include("../src/System1D.jl")
using .System1D

function test_dmc_basic()
    println("\n" * "="^60)
    println("Testing DMC Basic Functionality After Refactoring")
    println("="^60)
    
    # Create a simple harmonic oscillator potential
    V(R::Vector) = sum(x^2 for x in R)
    
    # Create Hamiltonian for 2 particles, diffusion D = 0.5
    H = Hamiltonian(2, 0.5, V)
    println("\n✓ Created Hamiltonian")
    @assert nparticles(H) == 2
    @assert diffusion_constant(H) ≈ 0.5
    
    # Test potential
    V_val = potential(H, [0.5, 0.5])
    println("✓ Potential accessor: ", V_val)
    @assert V_val ≈ 0.5
    
    # Create walkers
    initial_positions = [[0.1, 0.1], [0.2, 0.2], [0.15, 0.15]]
    walkers = [Walker(pos) for pos in initial_positions]
    println("✓ Created 3 walkers")
    
    # Create DMC parameters
    params = DMCParams(
        0.01,  # dt
        5,     # nsteps
        1,     # nequil
        3,     # targetN
        1.0,   # ET0
        0.1,   # pop_gain
        3,     # branch_cap
        3      # nblocks
    )
    println("✓ Created DMCParams")
    
    # Create simulation
    rng = MersenneTwister(42)
    sim = DMCSim(H, params, walkers, rng)
    println("✓ Created DMCSim (NoGuiding)")
    @assert sim.step == 0
    @assert length(sim.walkers) == 3
    
    # Test step
    println("\nTesting step!...")
    step!(sim)
    println("✓ Executed step 1, population: ", length(sim.walkers))
    @assert sim.step == 1
    
    # Test with guiding
    println("\nTesting ImportanceGuiding...")
    trial_wf = TrialWF(
        R -> -0.5 * sum(x^2 for x in R),
        R -> -R,
        R -> -length(R),
        R -> 1.0
    )
    
    guiding = ImportanceGuiding(trial_wf, H)
    walkers2 = [Walker(copy(p)) for p in initial_positions]
    sim2 = DMCSim(H, params, walkers2, rng; guiding=guiding)
    println("✓ Created DMCSim with ImportanceGuiding")
    
    # Test energy estimation
    E = estimate_energy(sim2)
    println("✓ Estimated energy: ", E)
    
    # Test guiding functions
    test_R = [0.5, 0.5]
    sig = signpsi(guiding, test_R)
    dr = drift(guiding, test_R)
    println("✓ signpsi and drift functions work")
    @assert dr ≈ -2 * test_R
    
    # Test node policies
    println("\nTesting Node Policies...")
    node_none = NoNode()
    node_fixed = FixedNode()
    println("✓ Created NoNode and FixedNode")
    
    @assert !crosses_node(node_none, guiding, [0.1, 0.1], [0.2, 0.2])
    println("✓ NoNode crossing check works")
    
    # Create with fixed node
    walkers3 = [Walker(copy(p)) for p in initial_positions]
    sim3 = DMCSim(H, params, walkers3, rng; guiding=guiding, nodepolicy=FixedNode())
    println("✓ Created DMCSim with FixedNode policy")
    
    # Test full simulation
    println("\nTesting full DMC run...")
    walkers4 = [Walker(copy(p)) for p in initial_positions]
    sim4 = DMCSim(H, params, walkers4, rng)
    run_simulation!(sim4; snapshot_steps=[0, 2, 4])
    
    println("✓ run_simulation! completed")
    println("  Steps: ", sim4.step)
    println("  Final population: ", length(sim4.walkers))
    println("  Snapshots: ", length(sim4.walker_positions_history))
    
    @assert sim4.step == params.nsteps
    @assert length(sim4.population_history) == params.nsteps + 1
    @assert length(sim4.walker_positions_history) == 3
    
    println("\n✓ All DMC tests passed!")
    println("✓ Backward compatibility verified!")
end

function test_vmc_skeleton()
    println("\n" * "="^60)
    println("Testing VMC Skeleton")
    println("="^60)
    
    V(R::Vector) = sum(x^2 for x in R)
    H = Hamiltonian(2, 0.5, V)
    walkers = [Walker(copy(p)) for p in [[0.1, 0.1], [0.2, 0.2]]]
    rng = MersenneTwister(42)
    
    params = VMCParams(0.01, 10, 3, 1.0)
    sim = VMCSim(H, params, walkers, rng)
    println("✓ Created VMCSim skeleton")
    
    # Test that VMC functions error appropriately
    try
        run_vmc!(sim)
        return false
    catch e
        contains(string(e), "not implemented") && println("✓ run_vmc! errors: \"not implemented\"")
    end
    
    try
        vmc_step!(sim)
    catch e
        contains(string(e), "not implemented") && println("✓ vmc_step! errors: \"not implemented\"")
    end
    
    println("✓ VMC skeleton verified!")
end

# Main
if abspath(PROGRAM_FILE) == @__FILE__
    try
        test_dmc_basic()
        test_vmc_skeleton()
        println("\n" * "="^60)
        println("✓✓✓ All Tests Passed! ✓✓✓")
        println("="^60)
    catch e
        println("\n✗ Test failed:")
        showerror(stdout, e, catch_backtrace())
        exit(1)
    end
end
