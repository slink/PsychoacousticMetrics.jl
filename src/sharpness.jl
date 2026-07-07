# Sharpness [acum] from ISO 532-1 Method 1 specific loudness.
#
# All weightings evaluate
#     S = 0.11 * Σᵢ N'(zᵢ) g(zᵢ) zᵢ dz / N,    zᵢ = 0.1·i (i = 1..240), dz = 0.1
# following MoSQITo sharpness_din_from_loudness.py (transcribed 2026-07-06
# from https://github.com/Eomys/MoSQITo, master). The 0.11 constant and the
# rectangle-rule grid are part of that transcription.

const Z_BARK = collect(0.1:0.1:24.0)
const WEIGHTINGS = (:din,)

# DIN 45692:2009 weighting (MoSQITo sharpness_din_from_loudness.py:124-125)
_g_din(z) = z > 15.8 ? 0.15 * exp(0.42 * (z - 15.8)) + 0.85 : 1.0

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

Reference: DIN 45692:2009. Implementation transcribed from MoSQITo
`sharpness_din_from_loudness.py` (Apache-2.0).
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

_weighting_function(weighting::Symbol, N::Float64) = _g_din
