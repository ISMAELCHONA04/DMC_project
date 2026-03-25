include(joinpath(@__DIR__, "gfmc_benchmark_helpers.jl"))
include(joinpath(@__DIR__, "gfmc_literature_report_helpers.jl"))

using Printf
using Statistics
using Plots

const THEORY_DOC_DIR = joinpath(PROJECT_ROOT, "docs")
const THEORY_ASSET_DIR = joinpath(THEORY_DOC_DIR, "assets", "gfmc_theory_comparisons")
const THEORY_DOC_PATH = joinpath(THEORY_DOC_DIR, "GFMC_BENCHMARK_THEORY_COMPARISONS.md")

mkpath(THEORY_ASSET_DIR)

default(; dpi=180)

function gallery_specs()
    return [
        (
            case_id="free_particle_ring",
            tier=:final,
            variant_id="main",
            variant_family="unguided",
            slug="free_ring",
            density_figure="density_density.png",
        ),
        (
            case_id="harmonic_oscillator_unguided",
            tier=:final,
            variant_id="main",
            variant_family="unguided",
            slug="ho_unguided",
            density_figure=nothing,
        ),
        (
            case_id="harmonic_oscillator_guided",
            tier=:final,
            variant_id="main",
            variant_family="guided",
            slug="ho_guided",
            density_figure=nothing,
        ),
        (
            case_id="two_particle_ho_guided",
            tier=:final,
            variant_id="main",
            variant_family="guided",
            slug="coupled_2ho_guided",
            density_figure=nothing,
        ),
        (
            case_id="two_particle_ho_fixed_node",
            tier=:final,
            variant_id="main",
            variant_family="fixed-node",
            slug="odd_ho_fixed_node",
            density_figure=nothing,
        ),
        (
            case_id="hydrogen_fixed_node",
            tier=:final,
            variant_id="main",
            variant_family="fixed-node",
            slug="hydrogen_fixed_node",
            density_figure=nothing,
        ),
        (
            case_id="hydrogen_unguided",
            tier=:final,
            variant_id="main",
            variant_family="unguided",
            slug="hydrogen_unguided",
            density_figure=nothing,
        ),
        (
            case_id="cosine_lattice_ring",
            tier=:final,
            variant_id="unguided",
            variant_family="unguided",
            slug="cosine_lattice_unguided",
            density_figure="density_unguided.png",
        ),
        (
            case_id="cosine_lattice_ring",
            tier=:sweep,
            variant_id="guided_dt0p0010",
            variant_family="guided",
            slug="cosine_lattice_guided_dt0p0010",
            density_figure="density_guided_dt0p0010.png",
        ),
    ]
end

function required_runs(specs)
    pairs = Tuple{String,Symbol}[]
    for spec in specs
        key = (spec.case_id, spec.tier)
        key in pairs || push!(pairs, key)
    end
    return pairs
end

function literature_meta(case_id::AbstractString, variant_family::AbstractString)
    for meta in literature_reference_catalog()
        if meta.case_id == case_id && meta.variant_family == variant_family
            return meta
        end
    end
    error("No literature metadata found for $(case_id):$(variant_family)")
end

function summary_row(case_id::AbstractString, tier::Symbol, variant_id::AbstractString)
    path = joinpath(benchmark_paths(case_id, tier).tables_dir, "summary.csv")
    rows = read_csv_rows(path)
    for row in rows
        row["variant_id"] == variant_id && return row
    end
    error("No summary row found for $(case_id):$(tier):$(variant_id)")
end

function history_rows(case_id::AbstractString, tier::Symbol, variant_id::AbstractString)
    path = joinpath(benchmark_paths(case_id, tier).tables_dir, "history_$(variant_id).csv")
    return read_csv_rows(path)
end

function density_figure_path(case_id::AbstractString, tier::Symbol, filename::AbstractString)
    return joinpath(benchmark_paths(case_id, tier).figures_dir, filename)
end

function save_history_comparison(spec, meta, row::Dict{String,String}, history::Vector{Dict{String,String}})
    taus = Float64[parse_float(entry, "tau") for entry in history]
    energies = Float64[parse_float(entry, "energy_mean") for entry in history]
    reference = reference_value(meta, row)
    errors = energies .- reference
    eq_tau = parse_int(row, "nequil") * parse_float(row, "dt")

    p_energy = plot(
        taus,
        energies;
        xlabel="imaginary time",
        ylabel="mean branch energy",
        title="$(meta.display_label): energy vs theory",
        label="simulation",
        color=:navy,
        linewidth=2.2,
    )
    hline!(p_energy, [reference]; color=:black, linestyle=:dash, linewidth=2.0, label="theory")
    vline!(p_energy, [eq_tau]; color=:gray50, linestyle=:dot, linewidth=1.6, label="equilibration cutoff")

    p_error = plot(
        taus,
        errors;
        xlabel="imaginary time",
        ylabel="energy minus theory",
        title="$(meta.display_label): error to theory",
        label="simulation - theory",
        color=:firebrick,
        linewidth=2.2,
    )
    hline!(p_error, [0.0]; color=:black, linestyle=:dash, linewidth=2.0, label="zero error")
    vline!(p_error, [eq_tau]; color=:gray50, linestyle=:dot, linewidth=1.6, label="equilibration cutoff")

    fig = plot(p_energy, p_error; layout=(1, 2), size=(1450, 480))
    path = joinpath(THEORY_ASSET_DIR, "$(spec.slug)_history.png")
    savefig(fig, path)
    return path
end

function copy_density_figure(spec)
    spec.density_figure === nothing && return nothing
    src = density_figure_path(spec.case_id, spec.tier, spec.density_figure)
    isfile(src) || error("Missing density comparison figure: $src")
    dest = joinpath(THEORY_ASSET_DIR, "$(spec.slug)_density.png")
    cp(src, dest; force=true)
    return dest
end

function build_curated_final_rows(specs)
    rows = NamedTuple[]
    citations = literature_citations()
    for spec in specs
        spec.tier === :final || continue
        meta = literature_meta(spec.case_id, spec.variant_family)
        include_in_final(meta) || continue
        row = summary_row(spec.case_id, spec.tier, spec.variant_id)
        ref_val = reference_value(meta, row)
        mean_energy = parse_float(row, "mean_energy")
        sem_energy = parse_float(row, "sem_energy")
        citation_text = join([citations[id].citation for id in meta.citation_ids], " | ")
        citation_urls = join([citations[id].url for id in meta.citation_ids], " | ")
        push!(rows, (
            case_id=spec.case_id,
            slug=spec.slug,
            display_label=meta.display_label,
            mean_energy=mean_energy,
            sem_energy=sem_energy,
            literature_reference=ref_val,
            energy_error=mean_energy - ref_val,
            absolute_error=abs(mean_energy - ref_val),
            theory_note=meta.theory_note,
            citations=citation_text,
            citation_urls=citation_urls,
        ))
    end
    return rows
end

function save_curated_final_dumbbell(rows)
    ordered = sort(collect(rows), by=row -> row.literature_reference)
    y_positions = collect(1:length(ordered))
    labels = [row.display_label for row in ordered]
    fig = plot(
        xlabel="energy",
        ylabel="benchmark case",
        title="GFMC final values vs paper-backed theory",
        yticks=(y_positions, labels),
        legend=:bottomright,
        size=(1300, 700),
    )
    for (y, row) in zip(y_positions, ordered)
        plot!(fig, [row.literature_reference, row.mean_energy], [y, y]; color=:gray60, linewidth=2.0, label="")
    end
    scatter!(fig, [row.literature_reference for row in ordered], y_positions; color=:black, marker=:diamond, markersize=7, label="theory")
    scatter!(fig, [row.mean_energy for row in ordered], y_positions; xerror=[row.sem_energy for row in ordered], color=:navy, marker=:circle, markersize=7, label="simulation")
    path = joinpath(THEORY_ASSET_DIR, "final_energy_vs_theory.png")
    savefig(fig, path)
    return path
end

function save_curated_final_error(rows)
    ordered = sort(collect(rows), by=row -> row.absolute_error, rev=true)
    fig = bar(
        [row.display_label for row in ordered],
        [row.energy_error for row in ordered];
        xlabel="benchmark case",
        ylabel="computed minus theory",
        title="GFMC final energy error vs paper-backed theory",
        legend=false,
        size=(1300, 650),
        xrotation=35,
        color=:firebrick,
    )
    hline!(fig, [0.0]; color=:black, linestyle=:dash, linewidth=1.5)
    path = joinpath(THEORY_ASSET_DIR, "final_energy_error_vs_theory.png")
    savefig(fig, path)
    return path
end

function write_theory_gallery_markdown(specs, final_rows, generated)
    citations = literature_citations()
    open(THEORY_DOC_PATH, "w") do io
        println(io, "# GFMC Benchmark Theory Comparisons")
        println(io)
        println(io, "This document collects the benchmark cases whose outputs can be compared directly to paper-backed theory references already cited in the repository.")
        println(io)
        println(io, "Regenerate the benchmark outputs and these curated figures from the repository root with:")
        println(io)
        println(io, "```bash")
        println(io, "julia --project=. Experiments/benchmarks/gfmc/generate_theory_comparison_gallery.jl")
        println(io, "```")
        println(io)
        println(io, "## Curated Aggregate Figures")
        println(io)
        println(io, "- [Final energies vs theory](assets/gfmc_theory_comparisons/final_energy_vs_theory.png)")
        println(io, "- [Final energy error vs theory](assets/gfmc_theory_comparisons/final_energy_error_vs_theory.png)")
        println(io)
        println(io, "## Included Cases")
        println(io)
        println(io, "| Case | History/Error Plot | Density Plot | Theory Reference |")
        println(io, "|---|---|---|---|")
        for spec in specs
            meta = literature_meta(spec.case_id, spec.variant_family)
            history_link = "assets/gfmc_theory_comparisons/$(spec.slug)_history.png"
            density_link = spec.density_figure === nothing ? "n/a" : "[density](assets/gfmc_theory_comparisons/$(spec.slug)_density.png)"
            citation_links = join(["[$(citations[id].short)]($(citations[id].url))" for id in meta.citation_ids], ", ")
            println(io, "| $(meta.display_label) | [history/error]($(history_link)) | $(density_link) | $(citation_links) |")
        end
        println(io)
        println(io, "## Final-Value Summary")
        println(io)
        println(io, "| Benchmark | Computed | Theory | Error |")
        println(io, "|---|---:|---:|---:|")
        for row in sort(collect(final_rows), by=row -> row.absolute_error, rev=true)
            println(io, @sprintf("| %s | %.8f +/- %.2e | %.8f | %.3e |",
                row.display_label,
                row.mean_energy,
                row.sem_energy,
                row.literature_reference,
                row.energy_error,
            ))
        end
        println(io)
        println(io, "## Notes")
        println(io)
        println(io, "- Only cases with direct paper-backed theory references are included here. The periodic-ion soft-Coulomb benchmarks remain exact within the repository model, but not against a matching published parameter set.")
        println(io, "- The cosine-lattice guided comparison uses the smallest time-step sweep run (`dt = 0.0010`) because the locked final benchmark for that case keeps only the unguided branch.")
        println(io, "- The unguided hydrogen case is included as a stress comparison against the odd-state 1D hydrogen reference; it is expected to fail badly, and that mismatch is the point of the plot.")
    end
    return THEORY_DOC_PATH
end

function main()
    specs = gallery_specs()

    println("Running paper-backed benchmark cases...")
    for (case_id, tier) in required_runs(specs)
        println("  - ", case_id, " [", tier, "]")
        run_benchmark_case(case_id, tier)
    end

    println("Generating curated comparison figures...")
    generated = NamedTuple[]
    for spec in specs
        meta = literature_meta(spec.case_id, spec.variant_family)
        row = summary_row(spec.case_id, spec.tier, spec.variant_id)
        history = history_rows(spec.case_id, spec.tier, spec.variant_id)
        history_path = save_history_comparison(spec, meta, row, history)
        density_path = copy_density_figure(spec)
        push!(generated, (slug=spec.slug, history_path=history_path, density_path=density_path))
    end

    final_rows = build_curated_final_rows(specs)
    final_plot = save_curated_final_dumbbell(final_rows)
    error_plot = save_curated_final_error(final_rows)
    doc_path = write_theory_gallery_markdown(specs, final_rows, generated)

    println("Curated theory gallery written.")
    println("Documentation: ", abspath(doc_path))
    println("Aggregate energy plot: ", abspath(final_plot))
    println("Aggregate error plot: ", abspath(error_plot))
    for item in generated
        println("Case history plot: ", abspath(item.history_path))
        item.density_path === nothing || println("Case density plot: ", abspath(item.density_path))
    end
end

main()
