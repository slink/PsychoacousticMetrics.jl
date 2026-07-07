# PsychoacousticMetrics.jl

Psychoacoustic sound-quality metrics for Julia, built on
[ZwickerLoudness.jl](https://github.com/slink/ZwickerLoudness.jl).

v0.1 implements **sharpness** [acum] per **DIN 45692:2009**, with the
Aures, von Bismarck, and Fastl weighting variants.

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

## Conformance

Tested against all 41 DIN 45692:2009 chapter-6 reference signals
(21 narrowband + 20 broadband) within the standard's tolerance of
±max(5 %, 0.05 acum), and cross-checked against
[MoSQITo](https://github.com/Eomys/MoSQITo) for all four weightings.

Note: total loudness N is computed as the Riemann sum over the 240-bin
specific loudness (`0.1 * sum`), not `ZwickerResult.loudness` — see the
`sharpness` docstring.

## Roadmap

Roughness, fluctuation strength, and psychoacoustic annoyance are planned
for this package. Time-varying loudness (ISO 532-1 Method 2) belongs in
ZwickerLoudness.jl.

## License

MIT. Test data and the Fastl weighting table are derived from MoSQITo
(Apache-2.0) — see `test/data/NOTICE`.
