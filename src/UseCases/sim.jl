# Use-cases: DMC simulation container and state

mutable struct DMCSim{H,W,RNG,G<:AbstractGuiding,NP<:AbstractNodePolicy}
    H::H
    params::DMCParams
    walkers::Vector{W}
    guiding::G # guiding wavefunction / policy
    nodepolicy::NP # nodal surface policy
    rng::RNG
    step::Int
    ET::Float64 # Reference energy
    ET_history::Vector{Float64}
    population_history::Vector{Int}
    energy_mean_history::Vector{Float64}
    walker_positions_history::Vector{Vector{Vector{Float64}}}
end

# Outer constructor with sensible defaults for histories & ET.
function DMCSim(h::H, params::DMCParams, walkers::Vector{W}, rng::RNG;
    guiding::G=NoGuiding(),
    nodepolicy::NP=NoNode(),
    step::Int=0) where {H,W,RNG,G<:AbstractGuiding,NP<:AbstractNodePolicy}
    return DMCSim{H,W,RNG,G,NP}(h, params, walkers, guiding, nodepolicy, rng, step, params.ET0, Float64[], Int[], Float64[], Vector{Vector{Vector{Float64}}}())
end

function DMCSim(h::H, params::DMCParams, walkers::Vector{W}, guiding::G, rng::RNG, step::Int=0) where {H,W,G<:AbstractGuiding,RNG}
    return DMCSim(h, params, walkers, rng; guiding=guiding, nodepolicy=NoNode(), step=step)
end

function DMCSim(h::H, params::DMCParams, walkers::Vector{W}, rng::RNG, step::Int=0) where {H,W,RNG}
    return DMCSim(h, params, walkers, rng; guiding=NoGuiding(), nodepolicy=NoNode(), step=step)
end

# Initialize simulation with given walkers.
function initialize!(sim::DMCSim, positions::Vector{<:AbstractVector})
    sim.walkers = [Walker(pos) for pos in positions]
    sim.step = 0
    sim.ET = sim.params.ET0
    sim.ET_history = Float64[]
    sim.population_history = Int[]
    sim.energy_mean_history = Float64[]
    sim.walker_positions_history = Vector{Vector{Vector{Float64}}}()
    return sim
end

# Record current walker positions (deep copy to freeze snapshot).
function record_positions!(sim::DMCSim)
    snap = [copy(position(w)) for w in sim.walkers]   # deep copy each R
    push!(sim.walker_positions_history, snap)
    return sim
end

# Update reference energy based on current population (growth estimator).
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

# Estimate mean energy from current walker set.
function estimate_energy(sim::DMCSim)::Float64
    return estimate_energy(sim, sim.guiding)
end

function estimate_energy(sim::DMCSim, ::NoGuiding)::Float64
    vals = Float64[potential(sim.H, position(w)) for w in sim.walkers]
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end

function estimate_energy(sim::DMCSim, g::ImportanceGuiding)::Float64
    vals = Float64[local_energy(g, position(w)) for w in sim.walkers]
    return isempty(vals) ? 0.0 : sum(vals) / length(vals)
end

# Record current state of the simulation for history tracking.
function record_state!(sim::DMCSim, record_positions::Bool=true)
    push!(sim.population_history, length(sim.walkers))
    push!(sim.energy_mean_history, estimate_energy(sim))
    push!(sim.ET_history, sim.ET)
    if record_positions
        record_positions!(sim)
    end
    return sim
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

# Branching factor for no importance sampling (uses potential energy).
function branching_factor(sim::DMCSim, Rold, Rnew, ::NoGuiding)
    dt = sim.params.dt
    ET = sim.ET
    Vold = potential(sim.H, Rold)
    Vnew = potential(sim.H, Rnew)
    P = exp(-0.5*dt * ((Vold + Vnew) - 2ET))
    return min(P, sim.params.branch_cap)
end

# Branching factor for importance sampling (uses local energy).
function branching_factor(sim::DMCSim, Rold, Rnew, g::ImportanceGuiding)
    dt = sim.params.dt
    ET = sim.ET
    ELold = local_energy(g, Rold)
    ELnew = local_energy(g, Rnew)
    P = exp(-0.5*dt * ((ELold + ELnew) - 2ET))
    return min(P, sim.params.branch_cap)
end

# Propose a diffusion move without importance sampling.
function propose_move(sim::DMCSim, w::Walker, ::NoGuiding)
    Dcoef = diffusion_constant(sim.H)   # from your Domain accessor
    dt = sim.params.dt
    Rold = position(w)                 # from your Domain accessor
    # Sample new position from diffusion PDF
    Rnew = Rold .+ sqrt(2 * Dcoef * dt) .* randn(sim.rng, length(Rold))
    return Rnew
end

# Drift-diffusion proposal with Metropolis accept/reject.
function propose_move(sim::DMCSim, w::Walker, g::ImportanceGuiding)
    D  = diffusion_constant(sim.H)
    dt = sim.params.dt
    Rold  = position(w)
    Fold = drift(g, Rold)
    # Forward proposal: drift + diffusion.
    Rnew = Rold .+ D*dt .* Fold .+ sqrt(2D*dt) .* randn(sim.rng, length(Rold))

    Fnew = drift(g, Rnew)
    denom = 4 * D * dt
    # Green's function ratio for detailed balance.
    log_gf = -sum(abs2, Rnew .- Rold .- D*dt .* Fold) / denom
    log_gb = -sum(abs2, Rold .- Rnew .- D*dt .* Fnew) / denom
    logpsi_old = g.trial.logpsi(Rold)
    logpsi_new = g.trial.logpsi(Rnew)
    log_ratio = log_gb - log_gf + 2 * (logpsi_new - logpsi_old)

    # Accept/reject step.
    if log(rand(sim.rng)) < min(0.0, log_ratio)
        return Rnew
    end
    return Rold
end



# Perform a single DMC step: propose moves then branch.
function step!(sim::DMCSim)
    new_walkers = similar(sim.walkers, 0)  # same element type, empty

    for w in sim.walkers
        Rold = position(w)
        Rnew = propose_move(sim, w, sim.guiding)

        if crosses_node(sim.nodepolicy, sim.guiding, Rold, Rnew)
            continue
        end

        P = branching_factor(sim, Rold, Rnew, sim.guiding)
        m = floor(Int, P + rand(sim.rng))

        for _ in 1:m
            push!(new_walkers, Walker(copy(Rnew)))  # copy is crucial for Vector R
        end
    end

    isempty(new_walkers) && error("No walkers left (extinction).")

    sim.walkers = new_walkers
    sim.step += 1
    return sim
end

# Main simulation loop (includes equilibration).
"""
    run_simulation!(sim::DMCSim; snapshot_steps=Int[])

Run the DMC simulation. When `snapshot_steps` is non-empty, walker positions are
recorded only at those integer steps (including step 0).
"""
function run_simulation!(sim::DMCSim; snapshot_steps::AbstractVector{<:Integer}=Int[])
    snapshot_set = isempty(snapshot_steps) ? nothing : Set(snapshot_steps)
    record_state!(sim, snapshot_set === nothing ? true : (0 in snapshot_set))

    for s in 1:sim.params.nsteps
        step!(sim) # preforms dmc step

        if s > sim.params.nequil
            update_ET!(sim)   # update reference energy after equilibration
        end
        record_state!(sim, snapshot_set === nothing ? true : (s in snapshot_set))
    end
    return sim
end


# Some convenience accessors

# Access current walker positions as a vector of vectors
function plot_snapshot_1d_density(snap::AbstractVector{<:AbstractVector};
    nbins::Int=200,
    xmin=nothing,
    xmax=nothing,
    normalize::Bool=true,
    title::AbstractString="Walker density")
    xs = [R[1] for R in snap]

    if xmin === nothing || xmax === nothing
        lo, hi = minimum(xs), maximum(xs)
        pad = 0.05 * (hi - lo + eps())   # small padding; eps avoids zero width
        xmin = (xmin === nothing) ? lo - pad : xmin
        xmax = (xmax === nothing) ? hi + pad : xmax
    end

    # `normalize=:pdf` makes it a density
    norm = normalize ? :pdf : :none
    return histogram(xs;
        bins=nbins,
        normalize=norm,
        xlabel="x",
        ylabel=normalize ? "density" : "count",
        title=title,
        label=false,
        xlims=(xmin, xmax)
    )
end


function plot_snapshot_1d_points(snap::AbstractVector{<:AbstractVector};
    title::AbstractString="Walker positions")
    xs = [R[1] for R in snap]
    return scatter(xs, zeros(length(xs));
        xlabel="x",
        yticks=false,
        ylabel="",
        title=title,
        label=false
    )
end
