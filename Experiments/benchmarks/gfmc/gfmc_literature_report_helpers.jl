using Printf
using Statistics
using Plots

const GFMC_BENCHMARK_DIR = @__DIR__
const GFMC_BENCHMARK_OUTPUTS = joinpath(GFMC_BENCHMARK_DIR, "outputs")
const LITERATURE_OUTPUT_ROOT = joinpath(GFMC_BENCHMARK_OUTPUTS, "literature")
const LITERATURE_TABLES_DIR = joinpath(LITERATURE_OUTPUT_ROOT, "tables")
const LITERATURE_FIGURES_DIR = joinpath(LITERATURE_OUTPUT_ROOT, "figures")

mkpath(LITERATURE_TABLES_DIR)
mkpath(LITERATURE_FIGURES_DIR)

default(; dpi=170)

function parse_csv_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    in_quotes = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            next_i = nextind(line, i)
            if in_quotes && next_i <= lastindex(line) && line[next_i] == '"'
                print(buf, '"')
                i = next_i
            else
                in_quotes = !in_quotes
            end
        elseif c == ',' && !in_quotes
            push!(fields, String(take!(buf)))
        else
            print(buf, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    return fields
end

function read_csv_rows(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return Dict{String,String}[]
    headers = parse_csv_line(lines[1])
    rows = Dict{String,String}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = parse_csv_line(line)
        row = Dict{String,String}()
        for (header, value) in zip(headers, values)
            row[header] = value
        end
        push!(rows, row)
    end
    return rows
end

parse_float(row::Dict{String,String}, key::AbstractString) = parse(Float64, row[key])
parse_int(row::Dict{String,String}, key::AbstractString) = parse(Int, row[key])

function output_summary_rows()
    return Dict(
        "smoke_all" => read_csv_rows(joinpath(GFMC_BENCHMARK_OUTPUTS, "suite", "all", "smoke", "tables", "summary.csv")),
        "sweep_accuracy" => read_csv_rows(joinpath(GFMC_BENCHMARK_OUTPUTS, "suite", "accuracy", "sweep", "tables", "summary.csv")),
        "final_accuracy" => read_csv_rows(joinpath(GFMC_BENCHMARK_OUTPUTS, "suite", "accuracy", "final", "tables", "summary.csv")),
        "final_stress" => read_csv_rows(joinpath(GFMC_BENCHMARK_OUTPUTS, "suite", "stress", "final", "tables", "summary.csv")),
    )
end

function literature_citations()
    return Dict(
        "schrodinger_1926_ii" => (
            short="Schrodinger 1926 (II)",
            citation="E. Schrodinger, \"Quantisierung als Eigenwertproblem (Zweite Mitteilung),\" Annalen der Physik 79, 489-527 (1926). DOI: 10.1002/andp.19263840602.",
            url="https://doi.org/10.1002/andp.19263840602",
        ),
        "vega_ojeda_guillen_mota_2024" => (
            short="Vega et al. 2024",
            citation="J. C. Vega, D. Ojeda-Guillen, and R. D. Mota, \"Exact solution of the isotropic and anisotropic Hamiltonian of two coupled harmonic oscillators,\" arXiv:2410.00021 (2024).",
            url="https://arxiv.org/abs/2410.00021",
        ),
        "loudon_1959" => (
            short="Loudon 1959",
            citation="R. Loudon, \"One-Dimensional Hydrogen Atom,\" American Journal of Physics 27, 649-655 (1959). DOI: 10.1119/1.1934950.",
            url="https://doi.org/10.1119/1.1934950",
        ),
        "palma_raff_2006" => (
            short="Palma and Raff 2006",
            citation="G. Palma and U. Raff, \"The one dimensional Hydrogen atom revisited,\" Canadian Journal of Physics 84, 787-800 (2006). DOI: 10.1139/P06-072.",
            url="https://doi.org/10.1139/P06-072",
        ),
        "bloch_1929" => (
            short="Bloch 1929",
            citation="F. Bloch, \"Uber die Quantenmechanik der Elektronen in Kristallgittern,\" Zeitschrift fur Physik 52, 555-600 (1929). DOI: 10.1007/BF01339455.",
            url="https://doi.org/10.1007/BF01339455",
        ),
        "slater_1937_periodic" => (
            short="Slater 1937",
            citation="J. C. Slater, \"Wave Functions in a Periodic Potential,\" Physical Review 51, 846-851 (1937). DOI: 10.1103/PhysRev.51.846.",
            url="https://doi.org/10.1103/PhysRev.51.846",
        ),
        "yin_2020_sinusoidal" => (
            short="Yin and Erwin 2020",
            citation="J. Yin and S. C. Erwin, \"Noninteracting Electrons in a Prototypical One-Dimensional Sinusoidal Potential,\" arXiv:2003.06647 (2020).",
            url="https://arxiv.org/abs/2003.06647",
        ),
        "slater_1929" => (
            short="Slater 1929",
            citation="J. C. Slater, \"The Theory of Complex Spectra,\" Physical Review 34, 1293-1322 (1929). DOI: 10.1103/PhysRev.34.1293.",
            url="https://doi.org/10.1103/PhysRev.34.1293",
        ),
        "girardeau_1960" => (
            short="Girardeau 1960",
            citation="M. Girardeau, \"Relationship between Systems of Impenetrable Bosons and Fermions in One Dimension,\" Journal of Mathematical Physics 1, 516-523 (1960). DOI: 10.1063/1.1703687.",
            url="https://doi.org/10.1063/1.1703687",
        ),
        "tonks_1936" => (
            short="Tonks 1936",
            citation="L. Tonks, \"The Complete Equation of State of One, Two and Three-Dimensional Gases of Hard Elastic Spheres,\" Physical Review 50, 955-963 (1936). DOI: 10.1103/PhysRev.50.955.",
            url="https://doi.org/10.1103/PhysRev.50.955",
        ),
        "lieb_liniger_1963" => (
            short="Lieb and Liniger 1963",
            citation="E. H. Lieb and W. Liniger, \"Exact Analysis of an Interacting Bose Gas. I. The General Solution and the Ground State,\" Physical Review 130, 1605-1616 (1963). DOI: 10.1103/PhysRev.130.1605.",
            url="https://doi.org/10.1103/PhysRev.130.1605",
        ),
    )
end

function literature_reference_catalog()
    return [
        (
            case_id="free_particle_ring",
            variant_family="unguided",
            display_label="Free Ring",
            reference_mode=:direct_value,
            reference_value=0.0,
            citation_ids=["schrodinger_1926_ii"],
            theory_note="Inference from Schrodinger's rotor/free-wave exact theory: the periodic free-particle ground state is uniform with E0 = 0 in the repository units.",
        ),
        (
            case_id="harmonic_oscillator_unguided",
            variant_family="unguided",
            display_label="HO Unguided",
            reference_mode=:direct_value,
            reference_value=0.5,
            citation_ids=["schrodinger_1926_ii"],
            theory_note="Direct exact harmonic-oscillator ground-state reference E0 = omega / 2 with omega = 1.",
        ),
        (
            case_id="harmonic_oscillator_guided",
            variant_family="guided",
            display_label="HO Guided",
            reference_mode=:direct_value,
            reference_value=0.5,
            citation_ids=["schrodinger_1926_ii"],
            theory_note="Same exact harmonic-oscillator reference as the unguided case, but with the exact Gaussian trial state.",
        ),
        (
            case_id="two_particle_ho_guided",
            variant_family="guided",
            display_label="Coupled 2-HO Guided",
            reference_mode=:direct_value,
            reference_value=0.5 * (1.0 + sqrt(1.0 + 2.0 * 0.7)),
            citation_ids=["vega_ojeda_guillen_mota_2024"],
            theory_note="Exact center-of-mass plus relative-mode energy for the benchmark parameters omega = 1 and kappa = 0.7.",
        ),
        (
            case_id="two_particle_ho_fixed_node",
            variant_family="fixed-node",
            display_label="Odd HO Fixed Node",
            reference_mode=:direct_value,
            reference_value=1.5,
            citation_ids=["schrodinger_1926_ii"],
            theory_note="Inference from the exact oscillator spectrum: the odd one-dimensional state has E = 3 omega / 2 for omega = 1.",
        ),
        (
            case_id="hydrogen_fixed_node",
            variant_family="fixed-node",
            display_label="1D Hydrogen Fixed Node",
            reference_mode=:direct_value,
            reference_value=-0.5,
            citation_ids=["loudon_1959", "palma_raff_2006"],
            theory_note="The benchmark targets the odd-parity hydrogenic branch; the n = 1 odd state has E = -1/2 in the convention used here. This is an inference from the cited odd-state treatments.",
        ),
        (
            case_id="hydrogen_unguided",
            variant_family="unguided",
            display_label="1D Hydrogen Unguided",
            reference_mode=:direct_value,
            reference_value=-0.5,
            citation_ids=["loudon_1959", "palma_raff_2006"],
            include_in_final=false,
            theory_note="Same odd-parity hydrogenic reference as the fixed-node case, used here only as a stress-path comparison because unguided propagation near the singular potential is intentionally unstable.",
        ),
        (
            case_id="cosine_lattice_ring",
            variant_family="unguided",
            display_label="Cosine Lattice Unguided",
            reference_mode=:row_reference_value,
            reference_value=NaN,
            citation_ids=["bloch_1929", "slater_1937_periodic", "yin_2020_sinusoidal"],
            theory_note="Paper-backed exact theory: a one-dimensional sinusoidal potential maps to the periodic Mathieu/Bloch problem. The plotted reference value is the benchmark's numerical evaluation of that exact theory at the repo parameters.",
        ),
        (
            case_id="cosine_lattice_ring",
            variant_family="guided",
            display_label="Cosine Lattice Guided",
            reference_mode=:row_reference_value,
            reference_value=NaN,
            citation_ids=["bloch_1929", "slater_1937_periodic", "yin_2020_sinusoidal"],
            theory_note="Same Mathieu/Bloch exact reference as the unguided cosine case, used here to stress-test the analytic guiding ansatz.",
        ),
    ]
end

function final_rows_by_key(rows::Vector{Dict{String,String}})
    lookup = Dict{Tuple{String,String},Dict{String,String}}()
    for row in rows
        lookup[(row["case_id"], row["variant_family"])] = row
    end
    return lookup
end

function reference_value(meta, row::Dict{String,String})
    if meta.reference_mode === :direct_value
        return meta.reference_value
    elseif meta.reference_mode === :row_reference_value
        return parse_float(row, "reference_value")
    end
    error("Unknown literature reference mode: $(meta.reference_mode)")
end

include_in_final(meta) = (:include_in_final in keys(meta) ? meta.include_in_final : true)

function write_rows_csv(path::AbstractString, rows::AbstractVector{<:NamedTuple})
    isempty(rows) && error("Cannot write CSV with no rows.")
    headers = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(headers), ","))
        for row in rows
            values = String[]
            for header in headers
                value = string(getproperty(row, header))
                value = replace(value, '"' => "\"\"")
                if occursin(',', value) || occursin('"', value) || occursin('\n', value)
                    push!(values, "\"" * value * "\"")
                else
                    push!(values, value)
                end
            end
            println(io, join(values, ","))
        end
    end
    return path
end

function build_final_comparison_rows(summary_rows::Dict{String,Vector{Dict{String,String}}})
    final_accuracy = final_rows_by_key(summary_rows["final_accuracy"])
    final_stress = final_rows_by_key(summary_rows["final_stress"])
    citations = literature_citations()

    rows = NamedTuple[]
    for meta in literature_reference_catalog()
        include_in_final(meta) || continue
        key = (meta.case_id, meta.variant_family)
        row = get(final_accuracy, key, nothing)
        row === nothing && (row = get(final_stress, key, nothing))
        row === nothing && continue

        ref_val = reference_value(meta, row)
        citation_text = join([citations[id].citation for id in meta.citation_ids], " | ")
        citation_urls = join([citations[id].url for id in meta.citation_ids], " | ")
        mean_energy = parse_float(row, "mean_energy")
        sem_energy = parse_float(row, "sem_energy")
        energy_error = mean_energy - ref_val
        push!(rows, (
            case_id=meta.case_id,
            variant_family=meta.variant_family,
            display_label=meta.display_label,
            tier=row["tier"],
            dt=parse_float(row, "dt"),
            mean_energy=mean_energy,
            sem_energy=sem_energy,
            literature_reference=ref_val,
            energy_error=energy_error,
            absolute_error=abs(energy_error),
            theory_note=meta.theory_note,
            citations=citation_text,
            citation_urls=citation_urls,
        ))
    end
    return rows
end

function build_series_rows(summary_rows::Dict{String,Vector{Dict{String,String}}})
    refs = literature_reference_catalog()
    all_rows = vcat(
        summary_rows["smoke_all"],
        summary_rows["sweep_accuracy"],
        summary_rows["final_accuracy"],
        summary_rows["final_stress"],
    )

    rows = NamedTuple[]
    for meta in refs
        for row in all_rows
            row["case_id"] == meta.case_id || continue
            row["variant_family"] == meta.variant_family || continue
            ref_val = reference_value(meta, row)
            mean_energy = parse_float(row, "mean_energy")
            sem_energy = parse_float(row, "sem_energy")
            push!(rows, (
                case_id=meta.case_id,
                variant_family=meta.variant_family,
                display_label=meta.display_label,
                tier=row["tier"],
                dt=parse_float(row, "dt"),
                mean_energy=mean_energy,
                sem_energy=sem_energy,
                literature_reference=ref_val,
                energy_error=mean_energy - ref_val,
                absolute_error=abs(mean_energy - ref_val),
            ))
        end
    end
    return rows
end

function save_final_dumbbell_plot(rows::AbstractVector{<:NamedTuple})
    ordered = sort(collect(rows), by=row -> row.literature_reference)
    y_positions = collect(1:length(ordered))
    labels = [row.display_label for row in ordered]

    fig = plot(
        xlabel="energy",
        ylabel="benchmark case",
        title="GFMC final results vs literature-backed references",
        yticks=(y_positions, labels),
        legend=:bottomright,
        size=(1300, 700),
    )

    for (y, row) in zip(y_positions, ordered)
        plot!(fig, [row.literature_reference, row.mean_energy], [y, y]; color=:gray60, linewidth=2.0, label="")
    end
    scatter!(fig, [row.literature_reference for row in ordered], y_positions; color=:black, marker=:diamond, markersize=7, label="literature/theory reference")
    scatter!(fig, [row.mean_energy for row in ordered], y_positions; xerror=[row.sem_energy for row in ordered], color=:navy, marker=:circle, markersize=7, label="computed final value")

    path = joinpath(LITERATURE_FIGURES_DIR, "final_vs_literature.png")
    savefig(fig, path)
    return path
end

function save_final_error_plot(rows::AbstractVector{<:NamedTuple})
    ordered = sort(collect(rows), by=row -> row.absolute_error, rev=true)
    fig = bar(
        [row.display_label for row in ordered],
        [row.energy_error for row in ordered];
        xlabel="benchmark case",
        ylabel="computed minus literature reference",
        title="GFMC final energy error vs literature",
        legend=false,
        size=(1300, 650),
        xrotation=35,
        color=:firebrick,
    )
    hline!(fig, [0.0]; color=:black, linestyle=:dash, linewidth=1.5)
    path = joinpath(LITERATURE_FIGURES_DIR, "final_error_vs_literature.png")
    savefig(fig, path)
    return path
end

function save_stress_series_plot(rows::AbstractVector{<:NamedTuple})
    selected_labels = Set([
        "HO Unguided",
        "HO Guided",
        "Cosine Lattice Unguided",
        "Cosine Lattice Guided",
        "1D Hydrogen Unguided",
        "1D Hydrogen Fixed Node",
    ])
    filtered = [row for row in rows if row.display_label in selected_labels]

    fig = plot(
        xscale=:log10,
        yscale=:log10,
        xlabel="GFMC time step",
        ylabel="absolute energy error",
        title="Paper-backed stress paths across benchmark tiers",
        legend=:topright,
        size=(1300, 700),
    )

    for label in unique(row.display_label for row in filtered)
        series = [row for row in filtered if row.display_label == label]
        sort!(series, by=row -> row.dt, rev=true)
        plot!(
            fig,
            [row.dt for row in series],
            [max(row.absolute_error, 1.0e-14) for row in series];
            marker=:circle,
            linewidth=2.2,
            label=label,
        )
    end

    path = joinpath(LITERATURE_FIGURES_DIR, "stress_paths_vs_literature.png")
    savefig(fig, path)
    return path
end

function write_literature_markdown(final_rows::AbstractVector{<:NamedTuple}, final_rows_path::AbstractString)
    citations = literature_citations()
    open(joinpath(LITERATURE_OUTPUT_ROOT, "report.md"), "w") do io
        println(io, "# GFMC Literature Comparison")
        println(io)
        println(io, "This report compares existing GFMC benchmark outputs against literature-backed exact references or literature-backed exact theories.")
        println(io)
        println(io, "Generated tables:")
        println(io, "- [final_comparison.csv](tables/final_comparison.csv)")
        println(io)
        println(io, "Generated figures:")
        println(io, "- `figures/final_vs_literature.png`")
        println(io, "- `figures/final_error_vs_literature.png`")
        println(io, "- `figures/stress_paths_vs_literature.png`")
        println(io)
        println(io, "## Final Comparison")
        println(io)
        println(io, "| Benchmark | Computed | Ref | Error | Note |")
        println(io, "|---|---:|---:|---:|---|")
        for row in sort(collect(final_rows), by=row -> row.absolute_error, rev=true)
            println(io, @sprintf("| %s | %.8f +/- %.2e | %.8f | %.3e | %s |",
                row.display_label,
                row.mean_energy,
                row.sem_energy,
                row.literature_reference,
                row.energy_error,
                row.theory_note,
            ))
        end
        println(io)
        println(io, "## Citations")
        println(io)
        for (id, citation) in sort(collect(citations), by=first)
            println(io, "- `", id, "`: ", citation.citation, " Link: ", citation.url)
        end
        println(io)
        println(io, "## Method Notes")
        println(io)
        println(io, "- For the canonical oscillator, rotor/ring, and odd-parity hydrogen benchmarks, the references are direct analytical values from the cited papers or straightforward inferences from those exact spectra.")
        println(io, "- The unguided hydrogen singular case is kept out of the final-value table on purpose because its extreme variance would dominate the final-energy scale; it is still included in `series_comparison.csv` and in the stress-path plot against the same odd-state reference.")
        println(io, "- For the cosine lattice benchmark, the cited papers establish the exact Mathieu/Bloch theory. The numeric reference value used here is the benchmark's own evaluation of that exact theory at the repository parameters.")
        println(io, "- The periodic-ion and bosonic scaffold benchmarks remain in the broader benchmark suite, but they are not included in this literature plot because their current references are exact within the repository's custom model rather than standalone values published for the exact same parameter set.")
    end
    return joinpath(LITERATURE_OUTPUT_ROOT, "report.md")
end
