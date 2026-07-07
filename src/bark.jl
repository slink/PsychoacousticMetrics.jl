# Bark <-> Hz and amplitude <-> dB conversions.
# Transcribed 2026-07-07 from MoSQITo commit d990c33f94f1:
#   mosqito/utils/conversion/freq2bark.py, bark2freq.py, amp2db.py, db2amp.py
# Tables: E. Zwicker, H. Fastl, Psychoacoustics (1990), table 6.1,
# linearly interpolated (numpy.interp semantics: clamped at both ends).

# numpy.interp equivalent: linear interpolation, clamped to endpoint values
# outside the table. xp must be sorted ascending.
function _lininterp(x::Real, xp::AbstractVector{<:Real}, yp::AbstractVector{<:Real})
    x <= xp[1] && return float(yp[1])
    x >= xp[end] && return float(yp[end])
    i = searchsortedlast(xp, x)
    x == xp[i] && return float(yp[i])
    t = (x - xp[i]) / (xp[i + 1] - xp[i])
    return yp[i] + t * (yp[i + 1] - yp[i])
end

const _BARK_TABLE_HZ = Float64[
    0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 510, 570, 630, 700,
    770, 840, 920, 1000, 1080, 1170, 1270, 1370, 1480, 1600, 1720, 1850,
    2000, 2150, 2320, 2500, 2700, 2900, 3150, 3400, 3700, 4000, 4400,
    4800, 5300, 5800, 6400, 7000, 7700, 8500, 9500, 10500, 12000, 13500,
    15500, 20000,
]
const _BARK_TABLE_Z = collect(0.0:0.5:24.5)   # numpy arange(0, 25, 0.5), 50 values

freq2bark(f::Real) = _lininterp(f, _BARK_TABLE_HZ, _BARK_TABLE_Z)
bark2freq(z::Real) = _lininterp(z, _BARK_TABLE_Z, _BARK_TABLE_HZ)

# MoSQITo amp2db: exact-zero amplitudes are replaced by 2e-12 before the log
# (upstream does this by mutating the input array; we keep the behaviour,
# not the mutation).
function _amp2db(amp::Real; ref::Real)
    ref == 0 && throw(ArgumentError("Reference must be different from 0"))
    a = amp == 0 ? 2e-12 : amp
    return 20 * log10(a / ref)
end

_db2amp(dB::Real; ref::Real) = 10.0^(0.05 * dB) * ref
