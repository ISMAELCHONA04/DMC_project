#!/usr/bin/env julia

using Printf
using Statistics

function parse_env_bool(key::AbstractString, default::Bool)
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return default
    low = lowercase(raw)
    if low in ("1", "true", "t", "yes", "y", "on")
        return true
    elseif low in ("0", "false", "f", "no", "n", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean in ENV[$key] = '$raw'."))
end

function parse_env_int(key::AbstractString, default::Int)
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return default
    try
        return parse(Int, raw)
    catch
        throw(ArgumentError("Invalid integer in ENV[$key] = '$raw'."))
    end
end

function parse_summary_energy(path::AbstractString, mode::AbstractString)
    isfile(path) || return NaN
    lines = readlines(path)
    length(lines) >= 2 || return NaN
    for line in lines[2:end]
        isempty(strip(line)) && continue
        cols = split(line, ',')
        length(cols) >= 5 || continue
        if cols[1] == mode
            try
                return parse(Float64, cols[5]) # Ebar column
            catch
                return NaN
            end
        end
    end
    return NaN
end

project_root = abspath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
base_script = abspath(joinpath(@__DIR__, "run_spinless_fermion_ring_1d_pbc_fixed_node.jl"))

nreplicas = parse_env_int("DMC_REPLICAS", max(1, min(Sys.CPU_THREADS, 8)))
nparallel = parse_env_int("DMC_PARALLEL_PROCS", max(1, min(nreplicas, Sys.CPU_THREADS)))
seed_base = parse_env_int("DMC_REPLICA_SEED_BASE", 10000)
show_progress = parse_env_bool("DMC_SHOW_PROGRESS", false)
make_plots = parse_env_bool("DMC_MAKE_PLOTS", false)
aggregate_mode = get(ENV, "DMC_AGGREGATE_MODE", "guided_fixed_node")

out_root = abspath(get(ENV, "DMC_REPLICA_OUT_ROOT", joinpath(@__DIR__, "..", "outputs", "replicas")))
mkpath(out_root)

println("============================================================")
println("Replica launcher: spinless fermion ring DMC")
println("project_root     = ", project_root)
println("base_script      = ", base_script)
println("nreplicas        = ", nreplicas)
println("parallel_procs   = ", nparallel)
println("seed_base        = ", seed_base)
println("make_plots       = ", make_plots)
println("show_progress    = ", show_progress)
println("out_root         = ", out_root)
println("aggregate_mode   = ", aggregate_mode)
println("============================================================")

nreplicas >= 1 || throw(ArgumentError("DMC_REPLICAS must be >= 1"))
nparallel >= 1 || throw(ArgumentError("DMC_PARALLEL_PROCS must be >= 1"))

function run_replica(replica_id::Int, project_root::String, base_script::String, out_root::String, seed_base::Int, make_plots::Bool, show_progress::Bool)
    replica_dir = joinpath(out_root, @sprintf("replica_%03d", replica_id))
    mkpath(replica_dir)
    log_path = joinpath(replica_dir, "run.log")
    summary_path = joinpath(replica_dir, "summary.csv")

    # Keep deterministic but separated RNG streams per replica.
    shift = 10_000 * (replica_id - 1)
    child_env = Dict(
        "GKSwstype" => get(ENV, "GKSwstype", "100"),
        "DMC_MAKE_PLOTS" => string(make_plots),
        "DMC_SHOW_PROGRESS" => string(show_progress),
        "DMC_REPLICA_ID" => string(replica_id),
        "DMC_OUTPUT_DIR" => replica_dir,
        "DMC_SUMMARY_OUT" => summary_path,
        "DMC_SEED_INIT" => string(seed_base + shift + 1),
        "DMC_SEED_GUIDED_FIXED" => string(seed_base + shift + 2),
        "DMC_SEED_GUIDED_NO_NODE" => string(seed_base + shift + 3),
        "DMC_SEED_UNGUIDED" => string(seed_base + shift + 4)
    )

    cmd = `julia --project=$project_root $base_script`
    merged_env = merge(Dict(string(k) => string(v) for (k, v) in ENV), child_env)

    t0 = time()
    process = nothing
    open(log_path, "w") do io
        process = run(pipeline(setenv(cmd, merged_env), stdout=io, stderr=io); wait=false)
        wait(process)
    end
    elapsed = time() - t0
    exit_code = process.exitcode
    ok = exit_code == 0
    return (
        replica=replica_id,
        ok=ok,
        exit_code=exit_code,
        elapsed_s=elapsed,
        log_path=log_path,
        summary_path=summary_path
    )
end

sem = Base.Semaphore(nparallel)
results = Vector{Any}(undef, nreplicas)
tasks = Task[]

for replica_id in 1:nreplicas
    task = @async begin
        Base.acquire(sem)
        try
            println(@sprintf("[launcher] starting replica %d/%d", replica_id, nreplicas))
            res = run_replica(replica_id, project_root, base_script, out_root, seed_base, make_plots, show_progress)
            results[replica_id] = res
            status = res.ok ? "ok" : "failed"
            println(@sprintf("[launcher] finished replica %d/%d (%s, exit=%d, %.1fs)", replica_id, nreplicas, status, res.exit_code, res.elapsed_s))
        finally
            Base.release(sem)
        end
    end
    push!(tasks, task)
end

foreach(wait, tasks)

ok_count = count(r -> r.ok, results)
println("============================================================")
println("Replica completion: ", ok_count, "/", nreplicas, " succeeded")

aggregate_rows = NamedTuple[]
energies = Float64[]
for r in results
    e = parse_summary_energy(r.summary_path, aggregate_mode)
    push!(aggregate_rows, (
        replica=r.replica,
        ok=r.ok,
        exit_code=r.exit_code,
        elapsed_s=r.elapsed_s,
        energy=e,
        summary_path=r.summary_path,
        log_path=r.log_path
    ))
    if r.ok && isfinite(e)
        push!(energies, e)
    end
end

aggregate_csv = joinpath(out_root, "aggregate.csv")
open(aggregate_csv, "w") do io
    println(io, "replica,ok,exit_code,elapsed_s,mode,energy,summary_path,log_path")
    for row in aggregate_rows
        @printf(
            io,
            "%d,%s,%d,%.3f,%s,%.16e,%s,%s\n",
            row.replica,
            string(row.ok),
            row.exit_code,
            row.elapsed_s,
            aggregate_mode,
            row.energy,
            row.summary_path,
            row.log_path
        )
    end
end

if !isempty(energies)
    emean = mean(energies)
    esem = length(energies) > 1 ? std(energies) / sqrt(length(energies)) : NaN
    println(@sprintf("Aggregated %s energy over %d replicas: %.10f +/- %.3e", aggregate_mode, length(energies), emean, esem))
else
    println("No finite energies found for mode='", aggregate_mode, "'. Check per-replica logs/summaries.")
end
println("Aggregate CSV: ", aggregate_csv)
println("============================================================")
