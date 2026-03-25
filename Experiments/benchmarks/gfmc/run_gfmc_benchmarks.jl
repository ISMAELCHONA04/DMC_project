include(joinpath(@__DIR__, "gfmc_benchmark_helpers.jl"))

function suite_name_from_selection(selection::AbstractString)
    raw = lowercase(strip(selection))
    raw = replace(raw, "," => "__")
    raw = replace(raw, " " => "_")
    isempty(raw) && return "default"
    return raw
end

function main(args::Vector{String})
    tier = length(args) >= 1 ? Symbol(lowercase(args[1])) : :smoke
    selection = length(args) >= 2 ? args[2] : "all"

    tier in (:smoke, :sweep, :final) || error("Tier must be smoke, sweep, or final")

    case_ids = benchmark_case_selection(selection)
    println("GFMC benchmark tier: ", tier)
    println("Selected cases: ", join(case_ids, ", "))

    results = run_benchmark_suite(case_ids, tier; suite_name=suite_name_from_selection(selection))

    println()
    println("Completed ", length(results.results), " benchmark case(s).")
    println("Summary CSV: ", abspath(results.summary_csv))
    if results.summary_figure !== nothing
        println("Summary figure: ", abspath(results.summary_figure))
    end
end

main(ARGS)
