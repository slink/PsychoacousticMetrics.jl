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

    Nf, Sf, Rf, FSf = Float64(N), Float64(S), Float64(R), Float64(FS)

    wS = Sf > 1.75 ? (Sf - 1.75) * log10(Nf + 10) / 4 : 0.0
    (isinf(wS) || isnan(wS)) && (wS = 0.0)

    wFR = (2.18 / Nf^0.4) * (0.4 * FSf + 0.6 * Rf)
    (isinf(wFR) || isnan(wFR)) && (wFR = 0.0)

    return Nf * (1 + sqrt(wS^2 + wFR^2))
end
