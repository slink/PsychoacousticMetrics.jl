# Thesis-conformance check against Fastl & Zwicker (2007) experimental
# fluctuation-strength curves, transcribed from Osses, García & Kohlrausch
# (2018 PhD thesis), Table B.1 (AM tones / AM broadband noise / FM tones,
# fmod in {1,2,4,8,16,32} Hz). Mirrors the D&W fig-3 conformance layer
# (test_conformance_dw_fig3.jl): most points gate at a relative tolerance
# against the published curve — `abs(computed - ref)/ref <= 0.30`, i.e. 30%
# OF THE REFERENCE value. Deliberately NOT `isapprox(...; rtol=0.30)`:
# isapprox scales by max(|computed|, |ref|), which would grant
# overestimating points an effective ~43%-of-reference bound; the
# reference-relative formula matches the informational tally below and the
# roughness layer's convention. Known model deviations are `@test_broken`,
# the AM-BBN row (stimulus mismatch, see below) is `@test_skip`, and an
# informational 30%-of-reference count over all 18 points is emitted
# regardless of gate outcome.
#
# NOTE on suite totals: Julia's Test stdlib buckets BOTH `@test_broken` and
# `@test_skip` under the "Broken" column, so this file contributes 11 to
# the suite's Broken count: 5 `@test_broken` (AM fmod=2 and fmod=32; FM
# fmod=8/16/32) plus 6 `@test_skip` (the AM-BBN row). Together with the 5
# pre-existing roughness-layer `@test_broken` points
# (test_conformance_dw_fig3.jl) the suite total is 16 Broken.
#
# --- AM tones: two adjudicated deviations (fmod=2 and fmod=32 Hz) ---
#
# Both points fail the 30%-of-reference gate and were checked, per the
# binding adjudication rule, against direct one-off runs of the real Octave
# oracle (SQAT @ 00b449e, same machine / Octave 11.3.0 / octave-signal 1.4.7
# as fs-oracle-pins.md) on the IDENTICAL `fs_am_tone(1000,fmod,70,4,44100)`
# samples this test uses (2026-07-08):
#
#   fmod=2:  computed 1.1096 vs. published 0.84 (rel dev 32.1%). Oracle
#     method=1 (time_varying, this test's default): FSmean = 1.109603 —
#     agrees with this package's 1.1096027 to displayed precision. (Oracle
#     method=0, stationary: 1.156835, even further from the published
#     value.) The reference model itself overestimates the published curve
#     by the same ~32% on this stimulus.
#   fmod=32: computed 0.0161 vs. published 0.06 (a real ~3.7x ratio).
#     Oracle method=1: FSmean = 0.016144 — agrees with this package's
#     0.016144 to displayed precision. (Oracle method=0: 0.015576.) The
#     reference model itself undershoots the published 0.06 by 3.7-3.8x.
#
# In both cases the real reference MODEL ITSELF fails the gate on the same
# signal — upstream model behavior at the curve's shoulder/tail, not a
# transcription bug in this package — so both are `@test_broken` (not a
# gate-widening). fmod in {1,4,8,16} pass the gate as measured (rel devs
# 1.9%, 5.6%, 0.9%, 6.1%), no adjudication needed.
#
# --- FM known deviation (thesis §B.4.1) ---
#
# The reference model overestimates fluctuation strength for FM tones with
# fmod > 4 Hz, with the band-pass peak shifted to 8 Hz rather than 4 Hz
# (thesis §B.4.1, cited directly — not re-derived). Task 1's oracle rig
# measured fm_tone_4hz (fmod = 4, the boundary case) at 2.212 vacil vs. the
# published 2.0 +/- 0.2 (JND) reference: 6% over, just outside the strict
# band, consistent with the overestimation already starting at the fmod=4
# boundary. fmod in {8,16,32} are therefore `@test_broken` per §B.4.1
# (documented upstream/model behavior); fmod in {1,2,4} pass the
# 30%-of-reference gate as measured (computed 0.8441, 1.4801, 2.2125 vacil;
# rel devs 0.7%, 26.5%, 10.6%).
#
# --- AM broadband-noise caveat (read before touching the BBN rows) ---
#
# Task 1's oracle rig measured the ORACLE ITSELF (SQAT's own
# FluctuationStrength_Osses2016, not this package) giving 3.86 vacil on the
# `am_bbn_4hz` fixture's noise stimulus vs. the published 1.80 vacil
# (fs-oracle-pins.md, "Step 7" / "Cross-checked against two more entries").
# The mechanism is stimulus mismatch, not model or transcription error:
# SQAT's own published validation curve was measured against
# `Fc-8010_BW-15980` (SQAT's original band-pass-around-8kHz noise, loaded
# from a 390 MB Zenodo wav dataset out of scope for this repo), while this
# package's AM BBN test stimulus (`fs_am_noise`, defined in fs_stimuli.jl) is
# a brick-wall 20 Hz-16 kHz low-pass-to-16kHz noise on a plain xorshift64
# stream — a materially broader-band, spectrally different signal. Because
# the oracle itself disagrees with the published reference by this much on
# OUR stimulus, the published Table B.1 AM-BBN row cannot gate our
# implementation via our stimulus: a pass would not confirm correctness and
# a fail would not indicate a bug.
#
# Before writing off all 6 BBN points, two cheap synthesis variants were
# tried (fmod in {1,2,4,8,16,32}, 60 dB SPL, 4 s @ 44100 Hz) to see whether a
# plausible modulation-depth/level convention change would bring OUR
# implementation's output (oracle-parity to 1e-9, per Task 4) near the
# published curve. Measured 2026-07-08:
#
#   fmod   as-is FS (vacil)   pre-mod x sqrt(1.5) FS   published ref
#    1      0.9670             1.0249                  1.12
#    2      3.1719             3.3122                  1.58
#    4      3.8621             4.0326                  1.80
#    8      3.8765             4.0471                  1.57
#   16      0.7950             0.8403                  0.48
#   32      0.0339             0.0374                  0.14
#
# ("pre-mod" = scaling the as-is (post-modulation-RMS-calibrated) signal by
# sqrt(1.5) =~ +1.76 dB, matching the AM-tone m=1 pre- vs. post-modulation
# calibration convention gap.) Neither variant tracks the published curve's
# SHAPE (fmod=2,4,8 overestimate by ~2x; fmod=32 undershoots by ~4x; fmod=1
# happens to land within 30% by coincidence, fmod=16 nearly so) — a uniform
# +-1.76 dB level shift cannot repair a shape mismatch. This confirms the
# deviation is dominated by the noise BANDWIDTH mismatch (16 kHz low-pass vs.
# SQAT's ~8-16 kHz band-pass, a materially different spectral content), not
# a modulation-depth or calibration convention, exactly as the mechanism
# above predicts. All 6 AM-BBN points are therefore marked `@test_skip`
# (not `@test_broken`): `@test_skip` never evaluates its expression and
# reports "broken" regardless of outcome, which is required here because the
# stimulus-mismatch reasoning means a pass or fail at any individual point is
# not meaningful either way (indeed fmod=1 passes the 30%-of-reference gate
# by coincidence — rel dev 13.7% — which would make `@test_broken` error
# with "Unexpected Pass"). The within-30%-informational tally below still
# counts these points. The crosscheck fixtures are unaffected (they compare
# oracle vs. Julia on the IDENTICAL as-is stimulus, a self-consistency
# check, not a validation one).

isdefined(@__MODULE__, :fs_am_tone) ||
    include(joinpath(@__DIR__, "support", "fs_stimuli.jl"))

# (fmod [Hz], FS [vacil]) — experimental references, Osses 2018 thesis Table
# B.1 (transcribed from Fastl & Zwicker 2007 via the thesis).
const FS_REF_AM_70DB = [(1.0, 0.39), (2.0, 0.84), (4.0, 1.25), (8.0, 1.30), (16.0, 0.36), (32.0, 0.06)]
const FS_REF_FM_70DB = [(1.0, 0.85), (2.0, 1.17), (4.0, 2.00), (8.0, 0.70), (16.0, 0.27), (32.0, 0.02)]
const FS_REF_AMBBN_60DB = [(1.0, 1.12), (2.0, 1.58), (4.0, 1.80), (8.0, 1.57), (16.0, 0.48), (32.0, 0.14)]

# AM tones: fmod=2 (overestimate, ~32%) and fmod=32 Hz (undershoot, ~3.7x)
# reproducibly deviate from the published curve in the real oracle itself
# too (one-off Octave runs, see header) -- upstream model behavior, not a
# transcription bug.
const _AM_KNOWN_BROKEN = Set([2.0, 32.0])
# FM fmod > 4 Hz: overestimation, band-pass peak shifted to 8 Hz (§B.4.1).
const _FM_KNOWN_BROKEN = Set([8.0, 16.0, 32.0])

# Gate (roughness/D&W convention): within 30% OF THE REFERENCE value.
# See header for why this is not isapprox(; rtol = 0.30).
const _GATE_RTOL = 0.30
_fs_gate(computed, ref) = abs(computed - ref) / ref <= _GATE_RTOL

@testset "thesis Table B.1 conformance (Fastl & Zwicker via Osses 2018)" begin
    within30 = 0
    total = 0

    @testset "AM tones, 70 dB SPL" begin
        for (fmod, fs_ref) in FS_REF_AM_70DB
            sig = fs_am_tone(1000.0, fmod, 70.0, 4.0, 44100.0)
            r = fluctuation_strength_osses(sig, 44100)
            total += 1
            within30 += _fs_gate(r.fluctuation_strength, fs_ref)
            @testset "fmod=$(fmod)Hz" begin
                if fmod in _AM_KNOWN_BROKEN
                    @test_broken _fs_gate(r.fluctuation_strength, fs_ref)
                else
                    @test _fs_gate(r.fluctuation_strength, fs_ref)
                end
            end
        end
    end

    @testset "FM tones, 70 dB SPL, fdev=700 Hz" begin
        for (fmod, fs_ref) in FS_REF_FM_70DB
            sig = fs_fm_tone(1500.0, fmod, 700.0, 70.0, 4.0, 44100.0)
            r = fluctuation_strength_osses(sig, 44100)
            total += 1
            within30 += _fs_gate(r.fluctuation_strength, fs_ref)
            @testset "fmod=$(fmod)Hz" begin
                if fmod in _FM_KNOWN_BROKEN
                    @test_broken _fs_gate(r.fluctuation_strength, fs_ref)
                else
                    @test _fs_gate(r.fluctuation_strength, fs_ref)
                end
            end
        end
    end

    @testset "AM broadband noise, 60 dB SPL" begin
        for (fmod, fs_ref) in FS_REF_AMBBN_60DB
            sig = fs_am_noise(16000.0, fmod, 60.0, 4.0, 44100.0)
            r = fluctuation_strength_osses(sig, 44100)
            total += 1
            within30 += _fs_gate(r.fluctuation_strength, fs_ref)
            @testset "fmod=$(fmod)Hz" begin
                # @test_skip, not @test_broken: the stimulus-mismatch
                # reasoning above means individual points may coincidentally
                # pass (fmod=1 does) or fail; @test_skip never evaluates the
                # expression, so it can't error on "Unexpected Pass".
                @test_skip _fs_gate(r.fluctuation_strength, fs_ref)
            end
        end
    end

    @info "thesis Table B.1: $within30/$total points within 30% of the Fastl & Zwicker reference curves (informational)"
end
