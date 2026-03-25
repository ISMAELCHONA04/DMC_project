# Formatting helpers for significant figures

using Printf

function format_sigfig(x; sigfigs::Int=6)
    if !isfinite(x)
        return string(x)
    end
    decimals = _decimals_for_sigfigs(x, sigfigs)
    x_r = round(x; digits=decimals)
    if decimals > 0
        return @sprintf("%.*f", decimals, x_r)
    else
        return @sprintf("%.0f", x_r)
    end
end

function _decimals_for_sigfigs(x, sigfigs::Int)
    if x == 0
        return sigfigs
    end
    exp = floor(Int, log10(abs(x)))
    return sigfigs - 1 - exp
end

function _sigfigs_for_uncertainty(err)
    abs_err = abs(err)
    exp = floor(Int, log10(abs_err))
    leading = abs_err / 10.0^exp
    sigfigs = leading < 3 ? 2 : 1
    decimals = -exp + (sigfigs - 1)
    return sigfigs, decimals
end

function format_with_uncertainty(value, err)
    if !isfinite(err)
        return format_sigfig(value; sigfigs=6), string(err)
    elseif err == 0
        decimals = _decimals_for_sigfigs(value, 6)
        value_r = round(value; digits=decimals)
        if decimals > 0
            return @sprintf("%.*f", decimals, value_r), @sprintf("%.*f", decimals, 0.0)
        else
            return @sprintf("%.0f", value_r), @sprintf("%.0f", 0.0)
        end
    end
    _, decimals = _sigfigs_for_uncertainty(err)
    value_r = round(value; digits=decimals)
    err_r = round(abs(err); digits=decimals)
    if decimals > 0
        return @sprintf("%.*f", decimals, value_r), @sprintf("%.*f", decimals, err_r)
    else
        return @sprintf("%.0f", value_r), @sprintf("%.0f", err_r)
    end
end

function format_pair_sigfig(a, b; sigfigs::Int=6)
    if !isfinite(b)
        return format_sigfig(a; sigfigs=sigfigs), string(b)
    elseif b == 0
        decimals = _decimals_for_sigfigs(a, sigfigs)
        a_r = round(a; digits=decimals)
        if decimals > 0
            return @sprintf("%.*f", decimals, a_r), @sprintf("%.*f", decimals, 0.0)
        else
            return @sprintf("%.0f", a_r), @sprintf("%.0f", 0.0)
        end
    end
    decimals = _decimals_for_sigfigs(b, sigfigs)
    a_r = round(a; digits=decimals)
    b_r = round(b; digits=decimals)
    if decimals > 0
        return @sprintf("%.*f", decimals, a_r), @sprintf("%.*f", decimals, b_r)
    else
        return @sprintf("%.0f", a_r), @sprintf("%.0f", b_r)
    end
end
