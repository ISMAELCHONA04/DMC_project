using Random

include("../../../../../src/System1D.jl")
using .System1D

# Single-particle 1D harmonic oscillator
H = Hamiltonian(1, 0.5, R -> 0.5 * R[1]^2)

# Initial positions and params
targetN = 10000
rng = MersenneTwister(42)
initial_positions = [[rand(rng) - 0.5] for _ in 1:targetN]
params = DMCParams(; dt=0.005, nsteps=50, nequil=50, targetN=targetN, ET0=0.5, branch_cap=10, nblocks=50)

sim = run_dmc(H, params, initial_positions; rng=rng)  # records snapshots at every time step by default

output_dir = joinpath(dirname(@__DIR__), "outputs", "animations")
mkpath(output_dir)
output_file = joinpath(output_dir, "dmc_single_particle_walkers.gif")
animate_dmc_walker_distribution_1d(
    sim;
    output_path=output_file,
    fps=10,
    nbins=150,
    density_smoothing_window=15,
    analytic_density = x -> sqrt(1.0 / (2*pi)) * exp(-0.5 * 1.0 * x^2),
    analytic_energy=0.5,
    size=(1100, 950)
)

println("Saved DMC animation: ", output_file)
