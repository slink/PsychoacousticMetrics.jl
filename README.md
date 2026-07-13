# PsychoacousticMetrics.jl

[![CI](https://github.com/slink/PsychoacousticMetrics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/slink/PsychoacousticMetrics.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/slink/PsychoacousticMetrics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slink/PsychoacousticMetrics.jl)

Psychoacoustic sound-quality metrics for Julia, built on
[ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).

v0.3 implements **sharpness** [acum] per **DIN 45692:2009** (with Aures,
von Bismarck, and Fastl variants), **roughness** [asper] per
**Daniel & Weber (1997)** as implemented by MoSQITo, and **fluctuation
strength** [vacil] per **Osses, García & Kohlrausch (2016)**.

## Quick start

```julia
using ZwickerLoudness, PsychoacousticMetrics

spl = [60, 62, 65, 68, 70, 72, 74, 75, 73, 71,
       69, 67, 65, 63, 61, 59, 57, 55, 53, 50,
       47, 44, 41, 38, 35, 32, 29, 26.0]   # 28 third-octave bands [dB]

result = zwicker_loudness(spl)
sharpness(result)                      # DIN 45692 sharpness [acum]
sharpness(result; weighting=:aures)    # :din (default), :aures/:bismarck/:fastl
```

From an audio file, compose with
[ZwickerLoudnessAudio.jl](https://github.com/slink/ZwickerLoudnessAudio.jl):

```julia
using ZwickerLoudnessAudio
sharpness(loudness_zwst("recording.wav"))
```

### Roughness

```julia
using Statistics
fs = 48000
t = range(0, 1, length=fs)
signal = @. (1 + sin(2π * 70 * t)) * sin(2π * 1000 * t)   # 100% AM tone
signal .*= 2e-5 * 10^(60/20) / std(signal)                 # 60 dB SPL

r = roughness_dw(signal, fs)
r.roughness            # overall roughness [asper] (~1.0 for this anchor)
r.roughness_over_time  # per 200 ms frame
r.specific_roughness   # 47 half-Bark channels × frames
```

### Fluctuation strength

```julia
using Statistics
fs = 44100
t = range(0, 3, length=3fs)
signal = @. (1 + sin(2π * 4 * t)) * sin(2π * 1000 * t)  # 100% AM tone, 4 Hz
signal .*= 2e-5 * 10^(60/20) / std(signal)               # 60 dB SPL

r = fluctuation_strength_osses(signal, fs)
r.fluctuation_strength            # ~1.01 vacil for this anchor
r.fluctuation_strength_over_time  # per 2 s frame, 90% overlap
r.specific_fluctuation_strength   # 47 half-Bark channels × frames
```

Anchor: a 1 kHz tone at 60 dB SPL, 100% amplitude-modulated at 4 Hz,
is defined to be 1 vacil. `fs` must be 44100 or 48000 Hz (the
reference model ships precomputed envelope filters for those rates
only) — resample other rates first.

### Psychoacoustic annoyance

Psychoacoustic annoyance (PA) [au] per **Widmann (1992)**, doctoral
thesis, TU München (formula p. 66) — widely known as "Zwicker
psychoacoustic annoyance" and commonly misattributed to Fastl &
Zwicker's *Psychoacoustics: Facts and Models* (the model is the one
that book popularized, not originated). See Lotinga & Torija (2025),
JASA 157(5):3282-3285, for the correction, and Fastl & Zwicker Ch. 16
as a widely-cited secondary description. Anchor: a 1 kHz tone at
40 dB SPL is defined to be 1 au (Widmann thesis p. 65).

Two surfaces are exported:

```julia
psychoacoustic_annoyance_widmann(N, S, R, FS)          -> Float64
psychoacoustic_annoyance_widmann(signal, fs, loudness)
    -> PsychoacousticAnnoyanceResult
```

The first is the pure formula on already-computed loudness `N`
[sone], sharpness `S` [acum], roughness `R` [asper], and fluctuation
strength `FS` [vacil]. The second composes this package's own
metrics (`sharpness`, `roughness_dw`, `fluctuation_strength_osses`)
plus a caller-supplied `ZwickerResult` for `N` on a single signal —
typically produced by
[ZwickerLoudnessAudio.jl](https://github.com/slink/ZwickerLoudnessAudio.jl):

```julia
using Statistics
using ZwickerLoudnessAudio, PsychoacousticMetrics

fs = 44100
t = range(0, 5, length = 5fs)
signal = sin.(2π * 1000 .* t)
signal .*= 2e-5 * 10^(40 / 20) / std(signal)   # 1 kHz tone, 40 dB SPL

loudness = loudness_zwst(signal, fs)
result = psychoacoustic_annoyance_widmann(signal, fs, loudness)
```

Output of this exact run:

```
result.pa                   = 1.0142921749776441
result.loudness             = 1.014
result.sharpness            = 1.0324589374115618
result.roughness            = 0.00015961445927174095
result.fluctuation_strength = 9.285797426192075e-5
result.convention           = stationary
```

(~1.4% from the 1 au anchor, in the same direction and order of
magnitude as SQAT's own reference implementation misses it by.)

The canonical Widmann model takes the 5th-percentile (`N5`, `S5`,
`R5`, `FS5`) of each metric's time-varying course over a signal; the
wrapper above uses whole-signal *stationary* values instead (ISO
532-1 Method 1 loudness, and this package's stationary sharpness/
roughness/fluctuation-strength), so it is a documented
**approximation** of the canonical percentile convention, not a
reproduction of it (`result.convention == :stationary` records this
explicitly). A percentile-based path is planned now that
[ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl)
v0.3.0 has shipped time-varying loudness, giving this package an
`N(t)` course to take the 5th percentile of.

`R` and `FS` enter the wrapper's formula with their raw, signed
values, mirroring SQAT's signal-level reference implementation,
which does not clamp its percentile components before combining
them. On near-stationary tones `roughness_dw`/
`fluctuation_strength_osses` can return a tiny negative value (model
noise around a true zero, not a bug — SQAT's own reference output
shows the identical artifact on the identical stimulus); the wrapper
passes that signed value straight into the formula rather than
clamping it, so `PsychoacousticAnnoyanceResult`'s `roughness`/
`fluctuation_strength` fields can be (very slightly) negative.

## Conformance

Tested against all 41 DIN 45692:2009 chapter-6 reference signals
(21 narrowband + 20 broadband) within the standard's tolerance of
±max(5 %, 0.05 acum), and cross-checked against
[MoSQITo](https://github.com/Eomys/MoSQITo) for all four weightings.

Note: total loudness N is computed as the Riemann sum over the 240-bin
specific loudness (`0.1 * sum`), not `ZwickerResult.loudness` — see the
`sharpness` docstring.

Roughness is tested against the Zwicker & Fastl reference curves on
MoSQITo's validation grid (7 carrier × 11 modulation frequencies, ±0.1
asper) and cross-checked against MoSQITo's `roughness_dw` on identical
signals.

Fluctuation strength is tested against Fastl & Zwicker's AM/FM
reference curves (Osses 2018 thesis Table B.1, ±30% of reference) and
cross-checked against SQAT's `FluctuationStrength_Osses2016`, run
under Octave, on identical signals. Two known deviations: the
reference model overestimates FM tones with fmod > 4 Hz (thesis
§B.4.1), and on the AM-tone curve the fmod = 2 and 32 Hz shoulder/tail
points diverge from the published values in the reference model
itself — both confirmed against a direct Octave oracle run (see the
conformance test's header, not a bug in this package).

## Roadmap

v0.4 adds **psychoacoustic annoyance** (Widmann, 1992), composing
this package's own loudness/sharpness/roughness/fluctuation-strength
metrics into a stationary approximation of the model — see
[Psychoacoustic annoyance](#psychoacoustic-annoyance) above.

Next: a **percentile-based** psychoacoustic annoyance path (`N5`,
`S5`, `R5`, `FS5`, matching the canonical Widmann convention exactly
rather than approximating it), now that ZwickerLoudness.jl v0.3.0 has
shipped time-varying loudness (ISO 532-1 Method 2) for this package
to draw a 5th-percentile `N(t)` from.

## License

MIT. Test data, the roughness implementation, and the Fastl/Bark/roughness
weighting tables are derived from MoSQITo (Apache-2.0). Fluctuation
strength's parameter tables and psychoacoustic annoyance's formula
constants, plus both metrics' cross-check fixtures, reference SQAT
(CC BY-NC 4.0) as factual data only — no code is reused. See
`test/data/NOTICE`.
