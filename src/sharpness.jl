# Sharpness [acum] from ISO 532-1 Method 1 specific loudness.
#
# All weightings evaluate
#     S = 0.11 * Σᵢ N'(zᵢ) g(zᵢ) zᵢ dz / N,    zᵢ = 0.1·i (i = 1..240), dz = 0.1
# following MoSQITo sharpness_din_from_loudness.py (transcribed 2026-07-06
# from https://github.com/Eomys/MoSQITo, master). The 0.11 constant and the
# rectangle-rule grid are part of that transcription.

const Z_BARK = collect(0.1:0.1:24.0)
const WEIGHTINGS = (:din, :aures, :bismarck, :fastl)

# DIN 45692:2009 weighting (MoSQITo sharpness_din_from_loudness.py:124-125)
_g_din(z) = z > 15.8 ? 0.15 * exp(0.42 * (z - 15.8)) + 0.85 : 1.0

# Aures 1985 — depends on total loudness N (line 127)
_g_aures(z, N) = 0.078 * (exp(0.171 * z) / z) * N / log(0.05 * N + 1)

# von Bismarck 1974 (lines 129-130)
_g_bismarck(z) = z > 15.0 ? 0.2 * exp(0.308 * (z - 15.0)) + 0.8 : 1.0

# Fastl & Zwicker 2007 — clamped linear interpolation, matching numpy.interp
function _g_fastl(z)
    z <= FASTL_X[1] && return FASTL_Y[1]
    z >= FASTL_X[end] && return FASTL_Y[end]
    i = searchsortedlast(FASTL_X, z)
    t = (z - FASTL_X[i]) / (FASTL_X[i + 1] - FASTL_X[i])
    return FASTL_Y[i] + t * (FASTL_Y[i + 1] - FASTL_Y[i])
end

"""
    sharpness(specific_loudness; weighting=:din) -> Float64
    sharpness(result::ZwickerResult; weighting=:din) -> Float64

Sharpness in acum from a 240-bin specific loudness N'(z) at 0.1 Bark
resolution (ISO 532-1:2017 Method 1), per DIN 45692:2009.

Total loudness is computed as `N = 0.1 * sum(specific_loudness)` — the
Riemann sum over the bins, NOT `ZwickerResult.loudness` (which the kernel
integrates analytically per spreading segment; the two differ by ~1 %).
Using the bin sum keeps numerator and denominator on the same
discretization. Returns 0.0 when N < 0.1 sone.

Four weightings are available, selected via `weighting`:
- `:din` — DIN 45692:2009 (the default).
- `:aures` — Aures 1985.
- `:bismarck` — von Bismarck 1974.
- `:fastl` — Fastl & Zwicker, *Psychoacoustics*, 2007.

Implementation transcribed from MoSQITo `sharpness_din_from_loudness.py`
(Apache-2.0).
"""
function sharpness(specific_loudness::AbstractVector{<:Real}; weighting::Symbol=:din)
    weighting in WEIGHTINGS ||
        throw(ArgumentError("weighting must be one of $(WEIGHTINGS), got $(repr(weighting))"))
    length(specific_loudness) == 240 ||
        throw(ArgumentError("specific_loudness must have 240 bins at 0.1 Bark resolution, got $(length(specific_loudness))"))
    all(>=(0), specific_loudness) ||
        throw(ArgumentError("specific_loudness values must be non-negative"))

    N = 0.1 * sum(specific_loudness)
    N < 0.1 && return 0.0

    g = _weighting_function(weighting, N)
    return 0.11 * sum(specific_loudness[i] * g(Z_BARK[i]) * Z_BARK[i] * 0.1 for i in 1:240) / N
end

sharpness(result::ZwickerResult; weighting::Symbol=:din) =
    sharpness(result.specific_loudness; weighting)

function _weighting_function(weighting::Symbol, N::Float64)
    weighting === :din && return _g_din
    weighting === :bismarck && return _g_bismarck
    weighting === :fastl && return _g_fastl
    return z -> _g_aures(z, N)
end
