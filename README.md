# PsychoacousticMetrics.jl

[![CI](https://github.com/slink/PsychoacousticMetrics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/slink/PsychoacousticMetrics.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/slink/PsychoacousticMetrics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/slink/PsychoacousticMetrics.jl)

Psychoacoustic sound-quality metrics for Julia, built on
[ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).

v0.2 implements **sharpness** [acum] per **DIN 45692:2009** (with Aures,
von Bismarck, and Fastl variants) and **roughness** [asper] per
**Daniel & Weber (1997)** as implemented by MoSQITo.

## Quick start

```julia
using ZwickerLoudness, PsychoacousticMetrics

spl = [60, 62, 65, 68, 70, 72, 74, 75, 73, 71,
       69, 67, 65, 63, 61, 59, 57, 55, 53, 50,
       47, 44, 41, 38, 35, 32, 29, 26.0]   # 28 third-octave bands [dB]

result = zwicker_loudness(spl)
sharpness(result)                      # DIN 45692 sharpness [acum]
sharpness(result; weighting=:aures)    # :din (default), :aures, :bismarck, :fastl
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

## Roadmap

Fluctuation strength and psychoacoustic annoyance are planned for this
package. Time-varying loudness (ISO 532-1 Method 2) belongs in
ZwickerLoudness.jl.

## License

MIT. Test data, the roughness implementation, and the Fastl/Bark/roughness
weighting tables are derived from MoSQITo (Apache-2.0) — see
`test/data/NOTICE`.
