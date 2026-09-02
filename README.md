# PsychoacousticMetrics.jl

[![CI](https://github.com/slink/PsychoacousticMetrics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/slink/PsychoacousticMetrics.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/slink/PsychoacousticMetrics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slink/PsychoacousticMetrics.jl)
[![DOI](https://zenodo.org/badge/1292553110.svg)](https://zenodo.org/badge/latestdoi/1292553110)

Psychoacoustic sound-quality metrics for Julia, built on
[ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).

v0.5 implements **sharpness** [acum] per **DIN 45692:2009** (with Aures,
von Bismarck, and Fastl variants), **roughness** [asper] per
**Daniel & Weber (1997)** as implemented by MoSQITo, **fluctuation
strength** [vacil] per **Osses, García & Kohlrausch (2016)**, and
**psychoacoustic annoyance** [au] per **Widmann (1992)** in both its
canonical percentile convention and a stationary approximation.

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

Three surfaces are exported:

```julia
psychoacoustic_annoyance_widmann(N, S, R, FS)          -> Float64
psychoacoustic_annoyance_widmann(signal, fs, tv::ZwickerTimeVaryingResult)
    -> PsychoacousticAnnoyanceResult   # convention = :percentile
psychoacoustic_annoyance_widmann(signal, fs, loudness::ZwickerResult)
    -> PsychoacousticAnnoyanceResult   # convention = :stationary
```

The first is the pure formula on already-computed loudness `N`
[sone], sharpness `S` [acum], roughness `R` [asper], and fluctuation
strength `FS` [vacil].

#### Percentile convention (canonical)

The canonical Widmann model takes the value exceeded 5 % of the time
(`N5`, `S5`, `R5`, `FS5`) of each metric's time-varying course over a
signal. The `ZwickerTimeVaryingResult` method implements exactly that:
`N5` comes from the caller-supplied time-varying loudness (ISO 532-1
Method 2, from
[ZwickerLoudnessAudio.jl](https://github.com/slink/ZwickerLoudnessAudio.jl)'s
`loudness_zwtv`), `S5` from DIN 45692 sharpness applied to every 2 ms
specific-loudness frame of that same result, and `R5`/`FS5` from
`roughness_dw`/`fluctuation_strength_osses`'s per-frame tracks. Each
percentile is taken over the metric's own native frame rate, as in the
reference implementation. `fs` must be 48000 Hz (the time-varying
loudness front end is 48 kHz-only) and `tv` must have been computed
from the same signal — a duration mismatch beyond one 2 ms loudness
block throws `ArgumentError`.

```julia
using Statistics
using ZwickerLoudnessAudio, PsychoacousticMetrics

fs = 48000
t = range(0, 5, length = 5fs)
signal = sin.(2π * 1000 .* t)
signal .*= 2e-5 * 10^(40 / 20) / std(signal)   # 1 kHz tone, 40 dB SPL

tv = loudness_zwtv(signal, fs)
result = psychoacoustic_annoyance_widmann(signal, fs, tv)
```

Output of this exact run:

```
result.pa                   = 1.003160292953811
result.loudness             = 1.000999999999856
result.sharpness            = 1.0357392751279124
result.roughness            = 0.00016138026780951667
result.fluctuation_strength = 0.002233844522856296
result.convention           = percentile
```

(0.3 % from the 1 au anchor on this 5 s tone; the suite's 4 s
cross-check case lands 0.1 % from it.)

Percentiles are Hyndman–Fan type 7 (`Statistics.quantile(track,
0.95)`, the same definition ZwickerLoudness.jl uses for its own `N5`).
SQAT's reference uses a nearest-rank selection with no interpolation
instead; on the 2 ms loudness/sharpness tracks the two agree to ~1e-9,
but on the coarse ~200 ms fluctuation-strength track (a dozen frames
on a 4 s signal) the gap is a full order-statistic step. The
cross-check tests measure and attribute that, the roughness-lineage
difference (MoSQITo vs. canonical MATLAB Daniel & Weber), and the
loudness-chain difference per component before any tolerance is set.

#### Stationary approximation

The `ZwickerResult` method composes this package's own metrics
(`sharpness`, `roughness_dw`, `fluctuation_strength_osses`) plus a
caller-supplied stationary `ZwickerResult` for `N` on a single signal —
typically produced by `loudness_zwst`:

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

This wrapper uses whole-signal *stationary* values (ISO 532-1 Method 1
loudness, and this package's stationary sharpness/roughness/
fluctuation-strength) in place of the canonical percentiles, so it is
a documented **approximation** of the percentile convention, not a
reproduction of it (`result.convention == :stationary` records this
explicitly). It works at 44.1 kHz as well as 48 kHz. On steady tones
the two conventions nearly coincide (~1.4 % PA difference on the
anchor tone, measured); on modulated signals they diverge by design,
since `N5` of a modulated signal sits near the modulation crests
rather than at the mean.

In both wrappers `R` and `FS` enter the formula with their raw,
signed values, mirroring SQAT's signal-level reference
implementation, which does not clamp its percentile components before
combining them. On near-stationary tones `roughness_dw`/
`fluctuation_strength_osses` can return a tiny negative value (model
noise around a true zero, not a bug — SQAT's own reference output
shows the identical artifact on the identical stimulus); the wrappers
pass that signed value straight into the formula rather than clamping
it, so `PsychoacousticAnnoyanceResult`'s `roughness`/
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
asper; 5 points at fc = 2 kHz are `@test_broken`, an upstream model
behavior that MoSQITo itself shows) and cross-checked against
MoSQITo's `roughness_dw` on identical signals. Against the Daniel &
Weber curves themselves the suite reports, informationally, 63 of 77
grid points within 30 % of the reference.

Fluctuation strength is tested against Fastl & Zwicker's AM/FM
reference curves (Osses 2018 thesis Table B.1) with a within-30 %-of-
reference gate, and cross-checked against SQAT's
`FluctuationStrength_Osses2016`, run under Octave, on identical
signals. The honest tally the suite prints is 8 of 18 Table B.1
points within 30 %: the reference model overestimates FM tones with
fmod > 4 Hz (thesis §B.4.1), on the AM-tone curve the fmod = 2 and
32 Hz shoulder/tail points diverge from the published values in the
reference model itself — both confirmed against a direct Octave
oracle run — and all six AM-broadband-noise points are skipped because
this package's noise stimulus does not reproduce SQAT's band-limited
one (see the conformance test's header). The Julia-vs-SQAT
cross-check on identical signals holds to rtol 1e-6, so these are
model/stimulus deviations, not bugs in this package.

Psychoacoustic annoyance is cross-checked against SQAT's
`PsychoacousticAnnoyance_Widmann1992` (an 840-case formula grid,
reproduced bit-exactly, plus six signal-level cases with SQAT's own
`N5`/`S5`/`R5`/`FS5` vendored so every deviation is attributed per
component) and gated at the published 1 au anchor: on the 4 s
cross-check case the percentile convention lands 0.1 % from it and
the stationary approximation 0.8 %; on the README's 5 s examples the
figures are 0.3 % and 1.4 %.

## Roadmap

v0.4 added **psychoacoustic annoyance** (Widmann, 1992) as a
stationary approximation; v0.5 adds the canonical **percentile**
convention (`N5`, `S5`, `R5`, `FS5`) on top of ZwickerLoudness.jl
v0.3's time-varying loudness — see
[Psychoacoustic annoyance](#psychoacoustic-annoyance) above. With
that, the core Zwicker-family metric set (loudness, sharpness,
roughness, fluctuation strength, annoyance) is complete.

Candidates for what comes next: test-suite hardening (Aqua.jl,
CompatHelper), and tonality.

## License

MIT. Test data, the roughness implementation, and the Fastl/Bark/roughness
weighting tables are derived from MoSQITo (Apache-2.0). Fluctuation
strength's parameter tables and psychoacoustic annoyance's formula
constants, plus both metrics' cross-check fixtures, reference SQAT
(CC BY-NC 4.0) as factual data only — no code is reused. See
`test/data/NOTICE`.
