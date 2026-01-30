# UseCases/dmc: Plotting utilities for DMC visualization

import Plots: histogram, scatter

"""
    plot_snapshot_1d_density(snap::AbstractVector{<:AbstractVector}; kwargs...)

Plot walker density as a histogram for a 1D snapshot.

# Arguments
- `snap`: Vector of walker positions (each position is a vector with first element x)
- `nbins::Int=200`: Number of histogram bins
- `xmin`, `xmax`: Custom axis limits (auto-detected if not provided)
- `normalize::Bool=true`: Whether to normalize as probability density
- `title::AbstractString`: Plot title
"""
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

"""
    plot_snapshot_1d_points(snap::AbstractVector{<:AbstractVector}; kwargs...)

Plot walker positions as a scatter plot for a 1D snapshot.

# Arguments
- `snap`: Vector of walker positions (each position is a vector with first element x)
- `title::AbstractString`: Plot title
"""
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
