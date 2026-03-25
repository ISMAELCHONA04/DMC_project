# UseCases/common: Step-level debug logging helpers

_debug_prefix(label::AbstractString) = isempty(label) ? "[debug]" : "[debug][$label]"

_debug_value_string(x::Integer) = string(x)
_debug_value_string(x::Bool) = x ? "true" : "false"

function _debug_value_string(x::AbstractFloat)
    xf = Float64(x)
    return isfinite(xf) ? string(round(xf; digits=6)) : string(xf)
end

_debug_value_string(x) = string(x)

function _emit_step_debug(
    io::IO,
    label::AbstractString,
    step::Integer,
    nsteps::Integer,
    step_elapsed::Real,
    elapsed::Real,
    fields::Pair...,
)
    frac = nsteps > 0 ? round((100.0 * step) / nsteps; digits=1) : 0.0
    print(
        io,
        _debug_prefix(label),
        " step ",
        step,
        "/",
        nsteps,
        " (",
        frac,
        "%, step ",
        round(Float64(step_elapsed); digits=3),
        "s, elapsed ",
        round(Float64(elapsed); digits=1),
        "s)",
    )
    for field in fields
        print(io, ", ", field.first, "=")
        print(io, _debug_value_string(field.second))
    end
    println(io)
    flush(io)
end
