# Tables for the Daniel & Weber roughness model.
# Transcribed 2026-07-07 from MoSQITo commit d990c33f94f1:
#   mosqito/utils/LTQ.py (reference="roughness"),
#   .../roughness_dw/_ear_filter_coeff.py, _gzi_weighting.py, _H_weighting.py

# Threshold in quiet for roughness [dB SPL] over Bark (LTQ.py roughness table)
const _LTQ_R_X = Float64[0, 0.01, 0.17, 0.8, 1, 1.5, 2, 3.3, 4, 5, 6, 8, 10,
    12, 13.3, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 24.5, 25]
const _LTQ_R_Y = Float64[130, 70, 60, 30, 25, 20, 15, 10, 8.1, 6.3, 5, 3.5,
    2.5, 1.7, 0, -2.5, -4, -3.7, -1.5, 1.4, 3.8, 5, 7.5, 15, 48, 60, 130]
_ltq_roughness(z::Real) = _lininterp(z, _LTQ_R_X, _LTQ_R_Y)

# Zwicker a0 outer/middle ear transmission [dB] over Bark (fig. 8.18)
const _A0_X = Float64[0, 10, 12, 13, 14, 15, 16, 16.5, 17, 18, 18.5, 19, 20,
    21, 21.5, 22, 22.5, 23, 23.5, 24, 25, 26]
const _A0_Y = Float64[0, 0, 1.15, 2.31, 3.85, 5.62, 6.92, 7.38, 6.92, 4.23,
    2.31, 0, -1.43, -2.59, -3.57, -5.19, -7.41, -11.3, -20, -40, -130, -999]
_a0_ear_db(z::Real) = _lininterp(z, _A0_X, _A0_Y)

# Aures modulation-depth weighting gzi over channel center [Bark]
const _GZI_X = collect(0.0:1.0:24.0)
const _GZI_Y = Float64[0.15, 0.26, 0.38, 0.47, 0.54, 0.65, 0.76, 0.83, 0.90,
    0.98, 0.98, 0.90, 0.80, 0.70, 0.62, 0.54, 0.49, 0.43, 0.39, 0.35, 0.30,
    0.30, 0.30, 0.30, 0.30]
_gzi(zi::Real) = _lininterp(zi, _GZI_X, _GZI_Y)

# H weighting anchor curves (D&W 1997 figure 2, as digitized by MoSQITo):
const _H2_X = Float64[0, 17, 23, 25, 32, 37, 48, 67, 90, 114, 171, 206, 247, 294, 358]
const _H2_Y = Float64[0, 0.8, 0.95, 0.975, 1, 0.975, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0]
const _H5_X = Float64[0, 32, 43, 56, 69, 92, 120, 142, 165, 231, 277, 331, 397, 502]
const _H5_Y = Float64[0, 0.8, 0.95, 1, 0.975, 0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0]
const _H16_X = Float64[0, 23.5, 34, 47, 56, 63, 79, 100, 115, 135, 159, 172,
    194, 215, 244, 290, 348, 415, 500, 645]
const _H16_Y = Float64[0, 0.4, 0.6, 0.8, 0.9, 0.95, 1, 0.975, 0.95, 0.9,
    0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
const _H21_X = Float64[0, 19, 44, 52.5, 58, 75, 101.5, 114.5, 132.5, 143.5,
    165.5, 197.5, 241, 290, 348, 415, 500, 645]
const _H21_Y = Float64[0, 0.4, 0.8, 0.9, 0.95, 1, 0.95, 0.9, 0.85, 0.8, 0.7,
    0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
const _H42_X = Float64[0, 15, 41, 49, 53, 64, 71, 88, 94, 106, 115, 137, 180,
    238, 290, 348, 415, 500, 645]
const _H42_Y = Float64[0, 0.4, 0.8, 0.9, 0.965, 0.99, 1, 0.95, 0.9, 0.85,
    0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]

"""
    _h_weighting(n, fs) -> Matrix{Float64}  (47 × n)

Modulation-frequency weighting functions Hᵢ per channel, transcribed from
MoSQITo `_H_weighting.py`. Upstream quirks preserved deliberately:
- anchor curves are COPIED to neighbouring channels (the code, unlike its
  own docstring, does no interpolation between anchors);
- `last = floor(502/fs·n)` from the H5 bound is reused for H16/H21/H42,
  truncating their 645 Hz tables at 502 Hz;
- python fills rows at 0-based bins j ∈ [cut, last) with cut = 2, i.e.
  1-based columns j+1 ∈ [3, last+1); bins below `cut` stay zero.
"""
function _h_weighting(n::Int, fs::Real)
    H = zeros(47, n)
    cut = 2
    function fill_row!(row::Int, hx::Vector{Float64}, hy::Vector{Float64}, fmax::Real)
        last = floor(Int, (fmax / fs) * n)
        for j in cut:(last - 1)               # python arange(cut, last), 0-based
            freq = j * fs / n
            H[row, j + 1] = _lininterp(freq, hx, hy)
        end
    end
    fill_row!(2, _H2_X, _H2_Y, 358)           # python H[1, j]
    fill_row!(5, _H5_X, _H5_Y, 502)           # python H[4, j]
    fill_row!(16, _H16_X, _H16_Y, 502)        # python H[15, j] — quirk: 502, not 645
    fill_row!(21, _H21_X, _H21_Y, 502)        # python H[20, j] — quirk: 502, not 645
    fill_row!(42, _H42_X, _H42_Y, 502)        # python H[41, j] — quirk: 502, not 645
    # channel copies (python lines 212-220)
    for r in (1, 3, 4); H[r, :] = H[2, :]; end
    for r in 6:15;  H[r, :] = H[5, :];  end
    for r in 17:20; H[r, :] = H[16, :]; end
    for r in 22:41; H[r, :] = H[21, :]; end
    for r in 43:47; H[r, :] = H[42, :]; end
    return H
end
