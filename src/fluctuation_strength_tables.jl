# Parameter tables and filtering primitives for fluctuation_strength_osses.
#
# Licensing note: SQAT (github.com/ggrecow/SQAT) is CC BY-NC 4.0 and the TUe
# originals are GPL-2.0/unlicensed, so — unlike the roughness path — NO code
# here is transcribed or translated from the reference. The numeric values
# below are unprotectable facts, each cited to the file that defines it at
# SQAT commit 00b449e. Algorithm structure follows the published description:
# A. Osses Vecchi, PhD thesis, TU/e (2018), Appendix B; Osses, García &
# Kohlrausch, Proc. Mtgs. Acoust. 28, 050005 (2016), doi:10.1121/2.0000410.

# gzi centre-frequency weighting, FS variant. Knot values from
# il_Get_gzi_fluctuation (FluctuationStrength_Osses2016.m @ SQAT 00b449e);
# thesis App. B.2 documents only the 4 anchor gains (13→1, 15→0.9,
# 17.5→0.7, 23.5→0.5). Beyond 24 Bark the reference maps interp1's NaN to
# the last gain (0.5); clamped _lininterp is equivalent.
const _GZI_FS_BARK = Float64[0, 1, 2.5, 4.9, 6.5, 8, 9, 10, 11, 11.5, 13, 15, 17.5, 24]
const _GZI_FS_GAIN = Float64[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9, 0.7, 0.5]
_gzi_fluctuation(z::Real) = _lininterp(z, _GZI_FS_BARK, _GZI_FS_GAIN)

# Simplified a0 transmission curve ('fluctuationstrength_osses2016'), i.e.
# Fastl's free-field curve with the ear-canal resonance removed (flat to
# 19 Bark). Knots from calculate_a0.m @ SQAT 00b449e.
const _A0_FS_BARK = Float64[0, 10, 19, 20, 21, 21.5, 22, 22.5, 23, 23.5, 24, 25, 26]
const _A0_FS_DB = Float64[0, 0, 0, -1.43, -2.59, -3.57, -5.19, -7.41, -11.3, -20, -40, -130, -999]
_a0_fs_db(z::Real) = _lininterp(z, _A0_FS_BARK, _A0_FS_DB)

# Hearing threshold: the FS knots from TerhardtExcitationPatterns.m's
# il_calculate_MinExcdB (@ SQAT 00b449e) are knot-for-knot identical to
# _ltq_roughness's table (src/roughness_dw_tables.jl, transcribed from
# MoSQITo LTQ.py reference="roughness") — all 27 (Bark, dB) pairs match
# exactly, so this aliases rather than vendoring a duplicate table.
# Dual citation: TerhardtExcitationPatterns.m @ SQAT 00b449e (il_calculate_MinExcdB)
# and mosqito/utils/LTQ.py (reference="roughness") @ MoSQITo d990c33f94f1.
const _fs_threshold_db = _ltq_roughness

# Cascaded-biquad filtering (direct-form II transposed, zero initial state),
# semantics of MATLAB/Octave sosfilt as used by il_modulation_depths
# (@ SQAT 00b449e). SOS rows are [b0 b1 b2 a0 a1 a2].
function _sosfilt(sos::AbstractMatrix{Float64}, x::AbstractVector{Float64})
    y = Vector{Float64}(x)
    for s in axes(sos, 1)
        a0 = sos[s, 4]
        b0, b1, b2 = sos[s, 1] / a0, sos[s, 2] / a0, sos[s, 3] / a0
        a1, a2 = sos[s, 5] / a0, sos[s, 6] / a0
        z1 = 0.0
        z2 = 0.0
        for n in eachindex(y)
            xn = y[n]
            yn = b0 * xn + z1
            z1 = b1 * xn - a1 * yn + z2
            z2 = b2 * xn - a2 * yn
            y[n] = yn
        end
    end
    return y
end
