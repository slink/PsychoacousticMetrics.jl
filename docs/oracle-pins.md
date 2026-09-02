# Oracle pins and measured facts

Tracked record of the reference implementations this package's fixtures
and tolerances were derived from, and of the behavioral facts pinned
against them. Source comments and test headers point here by section.

Everything below was measured once, on the pinned checkouts and tool
versions stated, and is quoted from the working notes taken at the time.
Where a note did not survive (see §3), that is stated rather than
reconstructed.

## 1. Pinned references

| Reference | Pin | License | Used for |
|---|---|---|---|
| SQAT (github.com/ggrecow/SQAT) | commit `00b449e40599f1c1ef4abe0596094552213d57eb` (commit date 2026-06-01), local checkout `/tmp/sqat-pinned` (override with `SQAT_DIR`) | CC BY-NC 4.0 | Fluctuation strength and psychoacoustic annoyance oracles: numeric outputs and parameter values only, no code transcribed (see `test/data/NOTICE`) |
| GNU Octave | 11.3.0 (arm64, darwin 25.4.0, `brew install octave`) with octave-signal 1.4.7 (octave-control 4.2.2 as dependency) | GPL | Runs the SQAT oracle |
| MoSQITo (github.com/Eomys/MoSQITo) | commit `d990c33f94f1` | Apache-2.0 | Sharpness and roughness references and transcriptions (see `test/data/NOTICE`) |

No shims are needed to run the SQAT checkout under Octave: every path,
including each model's `private/` directory, is added with plain
`addpath`. Octave does not enforce MATLAB's private-function scoping
rule, so `private/cos_ramp.m`, `private/Get_Hweight_fluctuation.m` and
`private/Get_gzi_fluctuation.m` are directly callable without modifying
the checkout.

## 2. Fluctuation strength (Osses, García & Kohlrausch 2016) oracle

### 2.1 Smoke run

Anchor stimulus (1 kHz carrier, 4 Hz modulation, m = 1, 60 dB SPL, 4 s
at 44100 Hz, synthesized with the convention in §2.2):
`FluctuationStrength_Osses2016(x, 44100, 1, 0)` ran without error and
gave `mean(OUT.InstantaneousFluctuationStrength) = 1.005363` vacil
(`FSmean` agrees).

### 2.2 Stimulus synthesis convention (binds `test/support/fs_stimuli.jl`)

SQAT's own validation scripts do not synthesize their stimuli; they read
pre-rendered `.wav` files from a Zenodo dataset
(doi:10.5281/zenodo.7933206, `validation_SQAT_v1_0.zip`, 390 MB, not
downloaded). The scripts reveal the calibration convention:
`insig = insig * 10^((dBFS_in-dBFS_out)/20)` with `dBFS_out = 94`, after
which `SPL = 20*log10(rms(insig)/2e-5)` is described as "verify final SPL
of the signal". SPL calibration therefore binds the RMS of the whole,
already-modulated signal.

Confirmed by sample inspection of the one stimulus SQAT ships in-repo,
`sound_files/reference_signals/RefSignal_FluctuationStrength_Osses2016.wav`
(60 dB SPL, 1 kHz, 100 % AM at 4 Hz, documented to give about 1 vacil):

- `fs = 44100`, float32 PCM, 220501 samples.
- `rms(wav) = 0.020000000015800…`, matching `2e-5*10^(60/20) = 0.02` Pa
  to 9 significant figures (post-modulation RMS; a carrier-based
  calibration would give `0.02/sqrt(0.75) = 0.02309`).
- `wav[1:10]` reversed equals `wav[end-9:end]` exactly, and
  `wav[1] = 0.0461876` is within 1e-7 of the global maximum. Only a cosine
  envelope times a cosine carrier puts both factors at their extrema at
  `t = 0`. Predicted peak `2k` for `k = p/sqrt(0.5+m^2/4)`, `m = 1`:
  `0.046188022` vs measured `0.046187568`.
- No onset/offset ramp: the global maximum is at sample 1, and local RMS
  is periodic with period 0.25 s across the whole file.

Pinned convention (`fs_am_tone`):
`x(t) = (1 + m*cos(2*pi*fmod*t)) * cos(2*pi*fc*t)`, scaled by
`k = p / sqrt(0.5 + m^2/4)` with `p = 2e-5 * 10^(SPL/20)`, so
`rms(k*x) == p`.

Cross-checked against two published points baked into SQAT's validation
scripts (Fastl & Zwicker reference values, not wav-derived):

| case | oracle [vacil] | published reference |
|---|---|---|
| `am_8hz_70db` (fc 1000, fmod 8, m 1, 70 dB) | 1.311230 | 1.30 ± 0.13 |
| `fm_tone_4hz` (fc 1500, fdev 700, fmod 4, 70 dB) | 2.212498 | 2.0 ± 0.2 |

FM and broadband-noise conventions are inferred by analogy (cosine
carrier; FM phase term `(fdev/fmod)*sin(2*pi*fmod*t)`; constant envelope
so `k = p*sqrt(2)`); no FM or noise reference wav ships with SQAT. The
`am_bbn_4hz` stimulus is deliberately a lowpass-to-16 kHz noise, not
SQAT's bandpass-around-8 kHz noise, so its oracle value (3.86 vacil vs
the published 1.80) is a self-consistency fixture only.

### 2.3 Buffer semantics (`buffer(x, N, round(0.9N), 'nodelay')`)

Provenance: octave-signal's `buffer` as called by the reference, measured
07-Jul-2026. For `N = 2*fs = 88200` (fs 44100), `V = round(0.9*88200) =
79380`, `hop = N - V = 8820`:

| L [s] | len = round(L·fs) | frames | `ceil((len-V)/hop)` | `floor((len-V)/hop)+1` |
|---|---|---|---|---|
| 2.0 | 88200 | 2 | 1 (wrong) | 2 |
| 3.0 | 132300 | 7 | 6 (wrong) | 7 |
| 3.05 | 134505 | 7 | 7 | 7 |
| 5.0 | 220500 | 17 | 16 (wrong) | 17 |

Pins:

- `nframes = floor((len - V)/hop) + 1`, not `ceil((len-V)/hop)`, which is
  wrong whenever `(len - V)` is an exact multiple of `hop`.
- Frame start times `t_start(k) = ((k-1)*hop + 1)/fs`, an arithmetic
  progression with step `hop/fs = 0.2` s and offset `1/fs =
  2.2676e-5` s (verified against every recorded row).
- Trailing partial frames are zero-padded, not dropped.
- Consequence for `method = 0` ("stationary"): `N = length(insig)`, so
  `len - V == hop` exactly and the stationary path always yields two
  frames (90 % overlapping), never one. `stationary_anchor` and
  `short_fallback` fixtures both have `length(fluct) == 2`.

### 2.4 Two-point spline extrapolation under Octave

`interp1([0 0.5], [0.1 0.2], 1, 'spline')` returns `NA` silently under
Octave (not the linear value 0.3). In `il_cross_correlation` (Chno = 47)
this cascades deterministically:

1. `ki(1,Chno-1) = interp1([0 0.5], ki(1,Chno-3:Chno-2), 1, 'spline')`
   returns `NA`, assignment succeeds.
2. `ki(1,Chno) = interp1([0 0.5], ki(1,Chno-2:Chno-1), 1, 'spline')` now
   has `NA` in its input and errors (`spline: requires at least 2
   points`).
3. The reference's `try/catch` runs the fallback for all four boundary
   assignments: `ki(1,Chno-1) = ki(1,Chno) = ki(1,Chno-2)`;
   `ki(2,1) = ki(2,2) = ki(1,3)`.

Pin: under Octave all four boundary corrections take the nearest-neighbor
copy, never the spline or linear-extrapolation value. Verified by direct
trace with random and fixed data; never observed to differ. The
instrumented reproduction is `scripts/dump_fs_stage.m`.

### 2.5 Content at or above ~24 Bark

15.8 kHz, 70 dB pure tone, 4 s at 44100 Hz, method 1: runs without error,
`FSmean = 0.0` exactly, `TimeAveragedSpecificFluctuationStrength` all zero
across the 47 channels. Pin: finite zero, silently, not an error.

### 2.6 Anchor at 48 kHz

Same anchor stimulus synthesized at 48000 Hz: `FSmean = 1.004505` vacil,
confirming the 48 kHz `Hweight-48000-Hz-{HP,LP}.mat` filters load.

### 2.7 `cos_ramp`

`cos_ramp(2*44100, 44100, 50, 50)` (`a = round(44100*50/1000) = 2205`):
first five samples `5.0748e-07, 2.0299e-06, 4.5674e-06, 8.1197e-06,
1.2687e-05`; last five `8.1197e-06, 4.5674e-06, 2.0299e-06, 5.0748e-07,
0.0`. Pins: `w(1) > 0` (equals `0.5*(1-cos(pi*1/2205))`), `w(a) = 1`
exactly, `w(a+1) = 1`, `w(end) = 0` exactly.

### 2.8 Filter coefficients (`src/fluctuation_strength_fir.jl`)

Extracted 2026-07-07 by `scripts/extract_fs_filters.jl` (values only,
`%.17g`):

- `Hweight-{fs}-Hz-HP.mat`: one SOS section; `Hweight-{fs}-Hz-LP.mat`:
  two SOS sections; `a0` column is 1 in every row, both sample rates.
  These are plain MATLAB v5 numeric data files read via `scipy.io.loadmat`.
- `calculate_a0(fs, 4096, 'fluctuationstrength_osses2016')` returns 4097
  taps for both sample rates, symmetric to machine precision (max
  `|B - fliplr(B)|` 5.6e-17 at 44100, 1.4e-17 at 48000).

The vendored constants in `src/fluctuation_strength_fir.jl` are the
authoritative copy. The intermediate JSON dumps the script writes
(`hweight_sos.json`, ~200 KB `a0_fir_b.json`) are regenerable scratch
output under `scratch/` and are not tracked.

### 2.9 Stage dump (anchor frame 1)

`scripts/dump_fs_stage.m` runs the pipeline's frame 1 for the
`anchor_44k` stimulus by calling the reference's own functions
(`cos_ramp`, `calculate_a0`, `TerhardtExcitationPatterns`,
`Get_Hweight_fluctuation`, `Get_gzi_fluctuation`) and dumps `h0`,
`mdepth_pre`, `mdepth_post`, `ki1`, `ki2`, `fi` (47 channels each), the
windowed frame's first and last 10 samples, and `post_a0_rms`
(`0.019247376652625438`). Measured 07-Jul-2026. The 47-channel literals
are vendored directly in `test/test_fluctuation_strength.jl` ("stage
dumps" testset) and `test/test_fluctuation_strength_stages.jl`; the JSON
itself is not tracked.

### 2.10 Generator run: `fs_mean` per fixture case

| case | fs | method | frames | fs_mean [vacil] |
|---|---|---|---|---|
| anchor_44k | 44100 | 1 | 12 | 1.005363 |
| anchor_48k | 48000 | 1 | 12 | 1.004505 |
| am_8hz_70db | 44100 | 1 | 12 | 1.311231 |
| fm_tone_4hz | 44100 | 1 | 12 | 2.212498 |
| am_bbn_4hz | 44100 | 1 | 12 | 3.862075 |
| stationary_anchor | 44100 | 0 | 2 (§2.3) | 1.004052 |
| short_fallback | 44100 | 1 → 0 (internal fallback, signal < 2 s) | 2 | 0.995734 |
| tone_25bark | 44100 | 1 | 12 | 0.0 (exact) |

### 2.11 Tolerance derivations

The per-case measured deviations behind every FS tolerance (anchor gate,
`rtol = 1e-6` cross-check, thesis-curve gate) are quoted inline in the
header comments of `test/test_fluctuation_strength.jl`,
`test/test_crosscheck_sqat_fs.jl` and `test/test_conformance_fs_thesis.jl`.
The full working log they were taken from was scratch and was not
retained.

## 3. Psychoacoustic annoyance, stationary convention (v0.4.0)

The working notes for the stationary-era rig (`pa-oracle-pins.md` and
its task reports) were scratch in a worktree that has since been removed
and are not recoverable. The facts the code and tests cite from them are
listed here; the negative-input pin was re-verified for §4.

- Formula constants (`S > 1.75` threshold, 2.18, the 0.4/0.6 weights,
  `N^0.4`, `log10(N+10)/4`) read from
  `PsychoacousticAnnoyance_Widmann1992_from_percentile.m`; the `log10`
  base follows SQAT and Fastl & Zwicker.
- The `S > 1.75` threshold is strict: `S == 1.75` takes the `w_S = 0`
  branch.
- The reference zeroes non-finite `w_S` and `w_FR`
  (`x(isinf(x)|isnan(x)) = 0`). With `N = 0, R = FS = 0` the naive `w_FR`
  is `Inf * 0 = NaN` before zeroing, so `N = 0` gives `PA = 0` exactly;
  all 120 formula-grid rows with `N == 0.0` have `pa == 0.0`.
- No epsilon guard is needed beyond that construct.
- The 840-case formula grid (`SQAT_PA_FORMULA_CASES`) is non-negative
  because the reference returns complex values for negative `N` (§4.3).
- Measured deviation of the stationary wrapper's components from SQAT's
  vendored percentiles, per case, is tabulated in
  `test/test_psychoacoustic_annoyance.jl` ("CONVENTION GAP"); the script
  log behind it was not retained.

## 4. Psychoacoustic annoyance, percentile convention (v0.5.0)

`/tmp/sqat-pinned` was already at the pinned SHA; no re-clone was needed.

### 4.1 SQAT signal-level PA source path

File `psychoacoustic_metrics/PsychoacousticAnnoyance_Widmann1992/PsychoacousticAnnoyance_Widmann1992.m`,
component tracks in call order (lines 151–185):

- **N(t)**: `Loudness_ISO532_1(insig, fs, LoudnessField, 2, time_skip, 0)`
  (method 2 = time-varying ISO 532-1:2017). `L.N5` is the percentile of
  `Total_Loudness`. Output rate 2 ms (`SR_LOUDNESS = 500`,
  `Loudness_ISO532_1.m:113`; `OUT.time = (0:…)' * 2e-3`, line 726).
- **S(t)**: `Sharpness_DIN45692_from_loudness(L.InstantaneousSpecificLoudness,
  'DIN45692', L.time, time_skip, 0)`, derived from the same time-varying
  specific loudness as N(t) (`ns_dec`, `[nTimeSteps x 240]`,
  `Loudness_ISO532_1.m:705`; Bark axis `(1:240)/10`, line 120), on
  N(t)'s own time axis (`Sharpness_DIN45692_from_loudness.m:108`).
  Weighting `g(z<15.8)=1; g(z>=15.8)=0.15*exp(0.42*(z-15.8))+0.85`
  (lines 167–168); `s = k*sum(SpecificLoudness.*g.*z.*0.10, 2)./loudness_sones`
  with `k = 0.11` (lines 85–86). The PA source's own comment (lines
  30–31) states the original PA sharpness weighting equals DIN 45692.
- **R(t)**: `Roughness_Daniel1997(insig, fs, time_skip, 0)` on the raw
  signal. Window `time_resolution = 0.2` s, `hopsize = N/2`
  (`Roughness_Daniel1997.m:69,78`), so ~100 ms hop.
- **FS(t)**: `FluctuationStrength_Osses2016(insig, fs, method_FS,
  time_skip, 0)` on the raw signal. Window 2 s
  (`FluctuationStrength_Osses2016.m:131`), `overlap = round(0.9*N)`
  (line 153), so ~200 ms hop.

The time-varying PA composition (lines 219–298, only when
`time_insig(end) >= 2`) spline-resamples the R and FS tracks onto
`L.time` (lines 254, 256) solely to form the optional `InstantaneousPA`
trace; the percentile scalars are each metric's own native-resolution
`X5`. The scalar composition (lines 282–298) is the same arithmetic as
`PsychoacousticAnnoyance_Widmann1992_from_percentile.m:75–91`, compared
line against line at this pin: strict `S > 1.75`,
`ws = (S-1.75)*log10(N+10)/4`, Inf/NaN zeroing,
`wfr = 2.18/N^0.4*(0.4*FS+0.6*R)`, `PA = N*(1+sqrt(ws^2+wfr^2))`.

### 4.2 Percentile mechanics

`N.N5`, `S.S5`, `R.R5`, `FS.FS5` (and `OUT.PA5`) each come from
`get_statistics(<track>(idx:end), '<MetricName>')`
(`Loudness_ISO532_1.m:740`, `Sharpness_DIN45692_from_loudness.m:116`,
`Roughness_Daniel1997.m:591`, `FluctuationStrength_Osses2016.m:235`,
`PsychoacousticAnnoyance_Widmann1992.m:316`). `get_statistics`
(`utilities/get_statistics.m:36–86`) maps `'5'` to
`get_exceeded_value(input, 5)`, which does not call `prctile`: it
computes `X_index = floor((100-PercentValue)/100 * size(input,1))`
(clamped to 1) and returns `sort(input)(X_index)`, a nearest-rank
selection with no interpolation (Hyndman–Fan type 1 style). No `prctile(`
call exists on this code path (the only ones in the checkout are in an
unrelated publication script and in `Roughness_ECMA418_2.m:677`).

Percentiles are taken on each metric's native track before the `interp1`
resampling used only for the `PA(t)` plot trace. This package uses
`Statistics.quantile(track, 0.95)` (type 7, linear interpolation); on
the 2 ms N(t)/S(t) tracks the two definitions agree to ~1e-9, on the
12–17-frame FS(t) track the gap can be a full order-statistic step.

### 4.3 Negative-input pin (re-verified, Octave 11.3.0)

```
PsychoacousticAnnoyance_Widmann1992_from_percentile(-1, 1, 0.3, 0.3)
  -> -1.2021+0.62199i   (iscomplex == 1)
```

`N^0.4` is complex for `N < 0` and propagates through
`sqrt(ws^2 + wfr^2)`. Negative `N` with `R = FS = 0` stays real (`-1`),
and negative `R`/`FS` with positive `N` stays real
(`(4, 1, -0.3, -0.3) -> 5.5025`), so the grid's non-negativity constraint
binds on `N`; signed `R`/`FS` feed through safely.

### 4.4 48 kHz signal cases

Added to `SIGNAL_CASE_DEFS` in `scripts/generate_sqat_pa_crosscheck.jl`:
`steady_1k_60db_48k` (1000 Hz, 0 Hz mod, 60 dB, 4 s, mdepth 0) and
`am_4hz_60db_48k` (1000 Hz, 4 Hz mod, 60 dB, 4 s, mdepth 1), 4 s to match
`fs_stimuli.jl`'s convention (the other cases are 5 s). SQAT's vendored
values:

| case | pa | N5 | S5 | R5 | FS5 |
|---|---|---|---|---|---|
| steady_1k_60db_48k | 4.039729273183468 | 4.034881749794696 | 1.0306524227322555 | 0.0006155141096380355 | 0.0014838887431786161 |
| am_4hz_60db_48k | 7.0104524211345876 | 4.706023832732632 | 1.035585574479048 | 0.01299911916394573 | 1.0239096429240029 |

Sanity vs the 44.1 kHz siblings: pa 4.0397 vs 4.0393 (steady), 7.0105
vs 7.0144 and FS5 1.0239 vs 1.0247 (AM). Regeneration byte-stability:
`test/data/sqat_pa_crosscheck.jl` pre-regeneration SHA-256
`e8df99132f03d0268325f165d098c969b4ba5bac8e102231dfad0b8121184eee`
(869 lines); the diff after regeneration showed only comment-line changes
and the two appended rows, with all 840 formula rows and 4 pre-existing
signal rows byte-identical.

### 4.5 Duration-guard bound (`ZwickerLoudnessAudio.loudness_zwtv` v0.3.0)

| signal duration | `length(signal)/fs` | `time_axis[end]` | gap [s] |
|---|---|---|---|
| 4.0 s (192000 samples) | 4.0 | 3.9984998124765596 | 0.0015001875234404 |
| 5.0 s (240000 samples) | 5.0 | 4.998499849984999 | 0.0015001500150014 |

Root cause: the front end's raw 0.5 ms axis is endpoint-inclusive
(`range(0.0, n/fs; length=n_time)`, `n_time = length(1:24:n)`), then
decimated by 4, dropping the last `(n_time − 1) mod 4` raw steps, so the
gap is quantized to {0, ~0.5, ~1.0, ~1.5} ms:

| n (samples) | n_time raw | (n_time−1) mod 4 | gap [s] |
|---|---|---|---|
| 192000 | 8000 | 3 | 0.0015001875234403883 |
| 192024 | 8001 | 0 | 0.0 |
| 192048 | 8002 | 1 | 0.0005000624921889241 |
| 192072 | 8003 | 2 | 0.0010001249687574543 |
| 192096 | 8004 | 3 | 0.0015001874297135842 |

Whole-second durations at 48 kHz always land on the ~1.5 ms worst case,
strictly below one 2 ms block, hence the wrapper's `> 0.002` s guard.

### 4.6 Measured deviations and tolerance derivation

Track lengths (native, per metric) on the three percentile fixture
cases: N(t)/S(t) 2500/2000/2000 frames at 2 ms; R(t) 49/39/39 frames;
FS(t) 17/12/12 frames. Mechanisms isolated by direct measurement:

- **N5**: loudness-implementation lineage (MoSQITo-lineage `loudness_zwtv`
  vs SQAT's ISO 532-1:2017 reference code). Percentile-definition share
  exactly 0.0. Substituting only our N5 into SQAT's components gives PA
  reldev 5.80e-3 / 4.43e-3 / 7.35e-3, i.e. N5 explains essentially the
  whole PA deviation on every case.
- **S5**: same DIN weighting, different loudness chain;
  percentile-definition share ≤ 1.8e-9. `S5 ≈ 1.03 < 1.75` on all cases,
  so S contributes exactly 0 to PA.
- **R5**: model lineage (canonical-MATLAB Daniel1997 vs MoSQITo
  `roughness_dw`); percentile-definition share ≤ 7e-12. The AM probe
  measured 4.54e-4 asper absolute (3.5 % relative).
- **FS5**: nearest-rank vs type-7 discretization on the 12–17-frame
  track. Applying SQAT's own selection rule to our track reproduces
  SQAT's vendored FS5 to 3.2e-10 (40 dB) and 2.0e-9 (60 dB); the 1.31e-2
  vacil deviation on the 60 dB case is percentile definition, not model
  drift. AM case: 6.08e-10.

Gates (each tolerance set after the attribution above):

- Substitution gate (SQAT's vendored N5/S5/R5/FS5 through our
  arithmetic vs SQAT's pa): reldev 0.0 exactly on all cases; test
  rtol 1e-12.
- Attribution gate (SQAT N5 + our S5/R5/FS5 vs SQAT pa): reldev
  4.68e-5 / 6.49e-3 / 2.15e-4; tolerances 1e-4 / 1.5e-2 / 5e-4.
- PA-level: measured 5.75e-3 / 2.04e-3 / 7.57e-3; tolerances
  1.2e-2 / 4.5e-3 / 1.5e-2.
- Anchor gate: `ANCHOR_ATOL = 2 × |SQAT pa(steady_1k_40db_48k) − 1| =
  2 × 6.866e-3 = 1.373e-2`; measured `|ours − 1| = 1.077e-3`.
- Percentile vs kernel N5: `quantile(tv.loudness_over_time, 0.95) − tv.N5
  = 0.0` exactly on all cases; test atol `4·eps(N5)`.
- Percentile vs stationary (steady 40 dB): `|pct PA − stat PA| =
  1.409e-2`; tolerance 3e-2.
- Duration gaps on matched pairs: 1.50e-3 s on all three cases.

Suite at the time: 1532 pass / 16 broken before the percentile work,
1577 pass / 16 broken after (+45 tests).
