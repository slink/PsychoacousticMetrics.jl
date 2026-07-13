"""
    PsychoacousticAnnoyanceResult

Result of the signal-level `psychoacoustic_annoyance_widmann` wrapper:
`pa` [au] and the four component values it was computed from —
`loudness` [sone], `sharpness` [acum], `roughness` [asper],
`fluctuation_strength` [vacil] — plus `convention`, a `Symbol` describing
how those components relate to the signal. The canonical Widmann model
percentile-izes all four components (not just loudness); `convention =
:stationary` documents that here all four are whole-signal stationary
values rather than percentiles (a future `:percentile` convention will
describe all four as percentiles once a time-varying loudness path
exists).
"""
struct PsychoacousticAnnoyanceResult
    pa::Float64
    loudness::Float64
    sharpness::Float64
    roughness::Float64
    fluctuation_strength::Float64
    convention::Symbol
end

"""
    psychoacoustic_annoyance_widmann(N, S, R, FS) -> Float64

Psychoacoustic annoyance (PA) [au] per Widmann (1992, doctoral thesis, TU
München; formula p. 66, reference signal p. 65), widely known as "Zwicker
psychoacoustic annoyance" and commonly misattributed to Fastl & Zwicker's
book — see Lotinga & Torija (2025), JASA 157(5):3282–3285, which clarifies
that this model is the one popularized (but not originated) by Fastl &
Zwicker, *Psychoacoustics: Facts and Models*, §16. Numeric constants below
are cited facts read from the pinned oracle,
`PsychoacousticAnnoyance_Widmann1992_from_percentile.m` @ SQAT commit
`00b449e40599f1c1ef4abe0596094552213d57eb` (CC BY-NC,
github.com/ggrecow/SQAT) — no code transcribed, see
`.superpowers/sdd/pa-oracle-pins.md`.

Takes already-computed metric values: `N` loudness [sone], `S` sharpness
[acum], `R` roughness [asper], `FS` fluctuation strength [vacil]. The
canonical model takes the 5th-percentile (`N5`, `S5`, `R5`, `FS5`) of each
metric's time-varying course over a signal; this surface is
percentile-agnostic — it computes the formula on whatever values the
caller supplies. See `psychoacoustic_annoyance_widmann(signal, fs, ...)`
for the wrapper that composes this ecosystem's four stationary metrics
(a documented approximation of the percentile convention).

    PA = N · (1 + √(w_S² + w_FR²))
    w_S  = (S − 1.75) · log₁₀(N + 10) / 4    for S > 1.75, else 0
    w_FR = 2.18 / N^0.4 · (0.4·FS + 0.6·R)

The log base in `w_S` is unspecified in Widmann's thesis; log₁₀ follows
Fastl & Zwicker and subsequent literature (SQAT's implementation likewise
uses `log10`). The `S > 1.75` threshold is STRICT (pinned against the
oracle: `S == 1.75` takes the `w_S = 0` branch). `w_FR` is computed by
direct division (not guarded against `N = 0`) and then has any
non-finite result zeroed — replicating the reference's
`wfr(isinf(wfr)|isnan(wfr)) = 0` construct exactly, including the `N = 0,
R = FS = 0` sub-case where the naive value is `Inf * 0 = NaN` before
zeroing (pinned in `.superpowers/sdd/pa-oracle-pins.md` Step 4). The net
effect is `N = 0 => PA = 0` exactly, consistent with the analytic limit
`PA ~ 2.18·N^0.6·(0.4·FS + 0.6·R) -> 0` as `N -> 0`.

Published anchor (Widmann thesis p. 65; SQAT's header states its
implementation reproduces it): a 1 kHz, 40 dB SPL tone gives PA = 1 au.

Throws `ArgumentError` if any of `N`, `S`, `R`, `FS` is negative or
non-finite. The reference has no such guard — a negative `N` or `S` flows
into `log10`/`^0.4` and yields a complex result in MATLAB/Octave; this
package rejects that input instead of returning a complex or NaN value.
"""
function psychoacoustic_annoyance_widmann(N::Real, S::Real, R::Real, FS::Real)
    isfinite(N) && N >= 0 || throw(ArgumentError("N must be non-negative and finite, got $N"))
    isfinite(S) && S >= 0 || throw(ArgumentError("S must be non-negative and finite, got $S"))
    isfinite(R) && R >= 0 || throw(ArgumentError("R must be non-negative and finite, got $R"))
    isfinite(FS) && FS >= 0 || throw(ArgumentError("FS must be non-negative and finite, got $FS"))

    return _pa_widmann_arithmetic(Float64(N), Float64(S), Float64(R), Float64(FS))
end

# The Widmann arithmetic without the public surface's non-negativity guard.
# Shared by the guarded scalar formula above and by the signal wrapper below.
# The wrapper needs the unguarded path because it mirrors the reference's
# signal-level behavior of feeding RAW SIGNED R/FS into the w_FR sum:
# SQAT's PsychoacousticAnnoyance_Widmann1992.m does not clamp its percentile
# components — a tiny negative FS5/R5 (near-zero modulation artifact on
# steady tones) contributes its signed value to 0.4·FS + 0.6·R, and only
# w_FR's square washes the sign of that total. The guard stays on the public
# scalar surface unchanged (negative N/S would go complex; negative R/FS on
# the PUBLIC surface remains rejected as out-of-domain caller input).
function _pa_widmann_arithmetic(N::Float64, S::Float64, R::Float64, FS::Float64)
    wS = S > 1.75 ? (S - 1.75) * log10(N + 10) / 4 : 0.0
    (isinf(wS) || isnan(wS)) && (wS = 0.0)

    wFR = (2.18 / N^0.4) * (0.4 * FS + 0.6 * R)
    (isinf(wFR) || isnan(wFR)) && (wFR = 0.0)

    return N * (1 + sqrt(wS^2 + wFR^2))
end

"""
    psychoacoustic_annoyance_widmann(signal::AbstractVector{<:Real}, fs::Real,
                                      loudness::ZwickerResult; pa_per_unit::Real=1.0)
        -> PsychoacousticAnnoyanceResult

Signal-level wrapper composing this ecosystem's four stationary metrics into
Widmann (1992) psychoacoustic annoyance (see the scalar-formula method above
for full model documentation and attribution).

- `N` = `loudness.loudness` — the caller supplies `loudness`, typically
  computed on the SAME `signal`/`fs` via ZwickerLoudnessAudio's
  `loudness_zwst` (ISO 532-1:2017 Method 1). This surface does not depend on
  ZwickerLoudnessAudio at runtime; it only needs the resulting
  `ZwickerLoudness.ZwickerResult`.
- `S` = `sharpness(loudness)` (default `:din` weighting — matches SQAT PA's
  default and Widmann's original weighting).
- `R` = `roughness_dw(signal, fs; pa_per_unit).roughness`.
- `FS` = `fluctuation_strength_osses(signal, fs; pa_per_unit).fluctuation_strength`.

`pa_per_unit` is forwarded to `roughness_dw` and `fluctuation_strength_osses`
only (component pass-through, per the resolved design — `loudness` already
carries its own calibration from however it was computed). If `signal` is
not already in pascals, pass the SAME `pa_per_unit` to `loudness_zwst` (or
whatever produced `loudness`) and to this function, so all four components
see a consistently-calibrated signal.

`convention = :stationary` on the returned result: the canonical Widmann
model percentile-izes all four components' time-varying course (see the
scalar formula's docstring); here all four are whole-signal stationary
values instead — a documented approximation, pending a time-varying-loudness
percentile path in a future release.

`R` and `FS` enter the formula RAW, sign included — mirroring the
reference's signal-level arithmetic: SQAT's
`PsychoacousticAnnoyance_Widmann1992.m` feeds its raw signed percentile
components straight into `w_FR = 2.18/N^0.4 · (0.4·FS + 0.6·R)` with no
clamp. Near-stationary tones can drive `roughness_dw`/
`fluctuation_strength_osses` to a tiny negative value — model noise around a
true zero, not a bug: SQAT's own oracle exhibits the identical artifact on
the exact same stimulus (`steady_1k_40db`'s vendored `FS5 = -0.001666...`,
`steady_1k_40db_48k`'s `FS5 = -0.001599...`, see
`test/data/sqat_pa_crosscheck.jl`). Such a value contributes its signed
amount to the `0.4·FS + 0.6·R` sum (only `w_FR`'s square washes the sign of
the total), so this wrapper does the same rather than routing through the
public scalar formula's non-negativity guard — that guard is a policy on
CALLER-supplied values of the scalar surface and remains unchanged there.
The raw values are also what the returned
`PsychoacousticAnnoyanceResult.roughness`/`.fluctuation_strength` fields
carry.

Error handling is inherited from the underlying components: `roughness_dw`
requires `fs >= 44100`, `fluctuation_strength_osses` requires `fs` in
`(44100, 48000)`, and both require `pa_per_unit > 0` — this wrapper adds no
guards of its own beyond what those calls already raise.

# Example

```julia
using PsychoacousticMetrics, ZwickerLoudnessAudio

signal = ...      # time-domain samples, calibrated so 1.0 == 1 pascal
fs = 44100
loudness = loudness_zwst(signal, fs)
result = psychoacoustic_annoyance_widmann(signal, fs, loudness)
result.pa
```
"""
function psychoacoustic_annoyance_widmann(signal::AbstractVector{<:Real}, fs::Real,
                                          loudness::ZwickerResult; pa_per_unit::Real=1.0)
    N = loudness.loudness
    S = sharpness(loudness)
    R = roughness_dw(signal, fs; pa_per_unit=pa_per_unit).roughness
    FS = fluctuation_strength_osses(signal, fs; pa_per_unit=pa_per_unit).fluctuation_strength
    # Raw signed R/FS, mirroring the reference's signal-level arithmetic
    # (see docstring); intentionally bypasses the public scalar surface's
    # non-negativity guard via the shared internal arithmetic.
    pa = _pa_widmann_arithmetic(N, S, R, FS)
    return PsychoacousticAnnoyanceResult(pa, N, S, R, FS, :stationary)
end
