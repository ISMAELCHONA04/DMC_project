using Random

include("../../../../../src/System1D.jl")
using .System1D

# Single-particle 1D harmonic oscillator with frequency omega.
omega = 1.0
H = Hamiltonian(1, 0.5, R -> 0.5 * omega^2 * R[1]^2)

# Exact guiding wavefunction for the 1D HO ground state in m = hbar = 1 units.
logpsi = R -> -0.5 * omega * R[1]^2
gradlogpsi = R -> [-omega * R[1]]
lapllogpsi = R -> -omega
trial = TrialWF(logpsi, gradlogpsi, lapllogpsi)
guiding = ImportanceGuiding(trial, H)

# Initial positions and params
targetN = 5000
rng = MersenneTwister(42)
initial_positions = [[2 * rand(rng) - 1] for _ in 1:targetN]
params = DMCParams(; dt=0.005, nsteps=100, nequil=50, targetN=targetN, ET0=0.5 * omega, branch_cap=10, nblocks=50)

sim = run_dmc(H, params, initial_positions; rng=rng, guiding=guiding)  # records snapshots at every time step by default

output_dir = joinpath(dirname(@__DIR__), "outputs", "animations")
mkpath(output_dir)
output_file = joinpath(output_dir, "dmc_single_particle_walkers_importance.gif")
animate_dmc_walker_distribution_1d(
    sim;
    output_path=output_file,
    fps=10,
    nbins=125,
    density_smoothing_window=15,
    analytic_density=x -> System1D.harmonic_oscillator_ground_density_1d(x; omega=omega),
    analytic_energy=0.5 * omega,
    size=(1100, 950),
    density_plot_style=:smoothed
)

println("Saved DMC importance-sampling animation: ", output_file)
