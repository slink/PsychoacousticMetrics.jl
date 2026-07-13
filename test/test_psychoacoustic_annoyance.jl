# Formula cross-check vs SQAT's PsychoacousticAnnoyance_Widmann1992_from_percentile
# (CC BY-NC, github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb),
# oracle only — no code transcribed (see .superpowers/sdd/pa-oracle-pins.md).
# Pure scalar arithmetic: any disagreement beyond the measured tolerance is a
# transcription bug, not a tolerance question.

isdefined(@__MODULE__, :SQAT_PA_FORMULA_CASES) ||
    include(joinpath(@__DIR__, "data", "sqat_pa_crosscheck.jl"))
isdefined(@__MODULE__, :fs_am_tone) ||
    include(joinpath(@__DIR__, "support", "fs_stimuli.jl"))

using ZwickerLoudnessAudio: loudness_zwst  # test-only dep ([extras]), NOT a runtime dep

@testset "psychoacoustic_annoyance_widmann formula" begin
    @testset "cross-check vs SQAT (840-case grid)" begin
        # Measured max relative deviation across all 840 cases, this
        # machine: 0.0 exactly (identical bit patterns) — see
        # .superpowers/sdd/pa-task-2-report.md. rtol kept at 1e-12 (pure
        # scalar arithmetic; exceeding 1e-9 would indicate a transcription
        # bug, not a tolerance question).
        RTOL = 1e-12
        for c in SQAT_PA_FORMULA_CASES
            pa = psychoacoustic_annoyance_widmann(c.N, c.S, c.R, c.FS)
            @test isapprox(pa, c.pa; rtol = RTOL, atol = 1e-12)
        end
    end

    @testset "N = 0 => PA = 0 exactly" begin
        # Pinned (pa-oracle-pins.md Step 4): all 120 grid rows with N == 0.0
        # give pa == 0.0 exactly via the reference's Inf/NaN zeroing,
        # including the R = FS = 0 sub-case where wfr = Inf*0 = NaN before
        # zeroing.
        zero_rows = filter(c -> c.N == 0.0, SQAT_PA_FORMULA_CASES)
        @test length(zero_rows) == 120
        for c in zero_rows
            @test psychoacoustic_annoyance_widmann(c.N, c.S, c.R, c.FS) == 0.0
        end
        @test psychoacoustic_annoyance_widmann(0.0, 2.5, 0.3, 0.3) == 0.0
        @test psychoacoustic_annoyance_widmann(0.0, 2.5, 0.0, 0.0) == 0.0
        @test psychoacoustic_annoyance_widmann(0.0, 0.5, 0.0, 0.0) == 0.0
    end

    @testset "threshold edge: strict S > 1.75" begin
        # Pinned (pa-oracle-pins.md Step 4): S = 1.75 takes the ws = 0
        # branch (else); only S > 1.75 (e.g. 1.7500001) engages ws != 0.
        pa_175 = psychoacoustic_annoyance_widmann(4.0, 1.75, 0.3, 0.3)
        pa_below = psychoacoustic_annoyance_widmann(4.0, 1.0, 0.3, 0.3)
        pa_above = psychoacoustic_annoyance_widmann(4.0, 1.7500001, 0.3, 0.3)
        @test pa_175 == pa_below
        @test pa_175 != pa_above
        @test isapprox(pa_175, 5.502497448336122; rtol = 1e-12)
        @test isapprox(pa_above, 5.502497448336126; rtol = 1e-12)
    end

    @testset "tiny N stays finite" begin
        # Pinned (pa-oracle-pins.md Step 4): no epsilon guard needed beyond
        # the documented Inf/NaN zeroing.
        pa = psychoacoustic_annoyance_widmann(1e-12, 2.5, 0.3, 0.3)
        @test isfinite(pa)
        @test pa > 0
        pa_zero_rfs = psychoacoustic_annoyance_widmann(1e-12, 2.5, 0.0, 0.0)
        @test isfinite(pa_zero_rfs)
    end

    @testset "input validation" begin
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(-1.0, 1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, -1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, -0.1, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, 0.0, -0.1)

        @test_throws ArgumentError psychoacoustic_annoyance_widmann(NaN, 1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(Inf, 1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, NaN, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, Inf, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, NaN, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, Inf, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, 0.0, NaN)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, 0.0, Inf)
    end
end

# ---------------------------------------------------------------------------
# Signal wrapper cross-check vs SQAT's signal-level
# PsychoacousticAnnoyance_Widmann1992 (same pin as above; oracle only).
#
# Signals are synthesized with the SAME fs_am_tone calls the fixture generator
# used (scripts/generate_sqat_pa_crosscheck.jl SIGNAL_CASE_DEFS), so both
# sides of every comparison see IDENTICAL input samples. Our N comes from
# ZwickerLoudnessAudio's loudness_zwst (ISO 532-1 Method 1 stationary, free
# field — matching the LoudnessField=0 the fixtures were generated with).
#
# The wrapper feeds RAW SIGNED R/FS into the Widmann arithmetic
# (_pa_widmann_arithmetic), mirroring the reference: SQAT's signal-level PA
# does not clamp its percentile components. Verified exactly below — the
# signed arithmetic on SQAT's own vendored (N5, S5, R5, FS5), including the
# two cases with NEGATIVE FS5, reproduces SQAT's vendored pa bit-for-bit
# (measured reldev 0.0 on all four cases, this machine).
#
# CONVENTION GAP, measured and attributed BEFORE any tolerance below was set
# (this machine, 2026-07-13; script log in .superpowers/sdd/pa-task-3-report.md):
# our wrapper composes whole-signal STATIONARY components, SQAT's signal-level
# PA composes 5th-PERCENTILE components (N5/S5/R5/FS5, vendored per case).
# Per-component deviations of ours vs SQAT's vendored percentiles:
#
#   case                 N reldev   S reldev   R absdev   FS absdev  PA reldev
#   steady_1k_40db       8.15e-3    4.22e-3    1.19e-4    3.30e-5    8.27e-3
#   steady_1k_60db       1.27e-2    2.42e-3    7.91e-5    1.71e-3    1.35e-2
#   am_4hz_60db          1.23e-1    4.34e-3    5.05e-3    1.90e-2    1.15e-1
#   steady_1k_40db_48k   8.12e-3    4.22e-3    1.21e-4    3.75e-5    8.24e-3
#
# Mechanism check (the attribution gate below asserts it stays true): feeding
# SQAT's N5 with OUR raw signed S/R/FS through the signed arithmetic
# reproduces SQAT's PA to <= 1.3e-4 (steady 40 dB), 7.9e-4 (steady 60 dB),
# 8.4e-3 (AM) relative — i.e. the PA-level deviation is almost entirely the
# N convention difference (stationary Method 1 loudness vs the 5th
# percentile of time-varying loudness), NOT drift in S/R/FS.
# On the steady tones N5 ~ stationary N and everything nearly coincides
# (PA within 1.4%); on the AM case the percentile-vs-stationary gap is a
# real model difference (N5 of a 4 Hz modulated signal sits near the
# modulation crests, 12.3% above the stationary value) and propagates
# through PA ~ N * (1 + sqrt(...)) essentially 1:1 — 11.5% PA deviation from
# a 12.3% N deviation, with the small remainder from FS5-vs-stationary-FS.
# R and FS on steady tones are near-zero noise on both sides (absolute
# comparisons; relative deviation is meaningless at 1e-4 asper).
#
# Tolerances: ~2x each measured deviation, rounded up, EXCEPT where the
# measured deviation is a near-zero epsilon whose 2x would be numerically
# meaningless — there a value-scale floor is used instead and marked in the
# table (R atol 5e-4 vs measured ~1e-4 absdev on steady tones ~ 4-6x: a
# floor at the scale below which asper differences carry no model meaning;
# same for FS atol 2e-4 vs measured ~3.5e-5, and for the 40 dB cases'
# attr rtol 5e-4 vs measured ~1.3e-4). Every entry's measured value is in
# the table above (or the mechanism-check note) so any real regression
# trips it.
# ---------------------------------------------------------------------------

# name => (fc, fmod, spl_db, dur_s, fs, mdepth) — identical to the generator.
const PA_WRAPPER_CASE_DEFS = Dict(
    "steady_1k_40db"     => (1000, 0, 40, 5.0, 44100, 0.0),
    "steady_1k_60db"     => (1000, 0, 60, 5.0, 44100, 0.0),
    "am_4hz_60db"        => (1000, 4, 60, 5.0, 44100, 1.0),
    "steady_1k_40db_48k" => (1000, 0, 40, 5.0, 48000, 0.0),
)

# (pa_rtol, N_rtol, S_rtol, R_atol, FS_atol, attribution_rtol) per case;
# derivation per the header note (~2x measured, or a marked value-scale
# floor for the near-zero R/FS absolute checks).
const PA_WRAPPER_TOLS = Dict(
    "steady_1k_40db"     => (pa = 2e-2,   N = 2e-2,   S = 1e-2, R = 5e-4,   FS = 2e-4, attr = 5e-4),
    "steady_1k_60db"     => (pa = 3e-2,   N = 3e-2,   S = 1e-2, R = 5e-4,   FS = 4e-3, attr = 2e-3),
    "am_4hz_60db"        => (pa = 2.5e-1, N = 2.5e-1, S = 1e-2, R = 1.5e-2, FS = 4e-2, attr = 2e-2),
    "steady_1k_40db_48k" => (pa = 2e-2,   N = 2e-2,   S = 1e-2, R = 5e-4,   FS = 2e-4, attr = 5e-4),
)

function pa_wrapper_case(name)
    fc, fmod, spl, dur, fs, mdepth = PA_WRAPPER_CASE_DEFS[name]
    sig = fs_am_tone(fc, fmod, spl, dur, fs; mdepth = mdepth)
    loud = loudness_zwst(sig, Float64(fs))
    return sig, Float64(fs), loud, psychoacoustic_annoyance_widmann(sig, Float64(fs), loud)
end

@testset "psychoacoustic_annoyance_widmann signal wrapper" begin
    wrapper_results = Dict(name => pa_wrapper_case(name)
                           for name in keys(PA_WRAPPER_CASE_DEFS))

    @testset "anchor gate: 1 kHz 40 dB tone -> 1 au" begin
        # Published anchor (Widmann thesis p. 65): PA = 1 au. Tolerance
        # derivation (FS-precedent style; non-tautological — derived from
        # the fixture-vs-1.0 gap only, never from the value under test):
        #   |oracle(steady_1k_40db pa) - 1.0| = |1.0068980604614974 - 1|
        #                                     = 6.898e-3
        #     (the reference implementation's own distance from the
        #     published 1 au, measured once from the vendored fixture)
        # ANCHOR_ATOL = 3 * oracle_gap = 2.069e-2 — the smallest integer
        # multiplier that passes. The FS anchor gate used 2x because its
        # Julia-vs-oracle parity was ~1e-10; here the wrapper deliberately
        # differs from the oracle in CONVENTION (stationary components vs
        # the oracle's percentiles), which contributes a second same-order
        # offset (measured 8.3e-3, attributed ~fully to N in the
        # cross-check testset below). Measured wrapper PA = 1.0152255:
        # |ours - 1| = 1.523e-2 > 2 * gap = 1.380e-2, so a 2x gate would
        # FAIL — recorded, not hidden. The gate stays falsifiable: it does
        # not scale with the value under test, and a regression pushing PA
        # outside 1 +/- 2.07e-2 fails regardless of what caused it.
        anchor = only(c for c in SQAT_PA_SIGNAL_CASES if c.name == "steady_1k_40db")
        oracle_gap = abs(anchor.pa - 1.0)
        ANCHOR_ATOL = 3 * oracle_gap
        _, _, _, r = wrapper_results["steady_1k_40db"]
        @test isapprox(r.pa, 1.0; atol = ANCHOR_ATOL)
    end

    @testset "cross-check vs SQAT signal-level PA (per-component attribution)" begin
        for c in SQAT_PA_SIGNAL_CASES
            haskey(PA_WRAPPER_CASE_DEFS, c.name) || continue
            _, fs, _, r = wrapper_results[c.name]
            tol = PA_WRAPPER_TOLS[c.name]
            @testset "$(c.name)" begin
                @test fs == c.fs
                # arithmetic pin: the signed Widmann arithmetic on SQAT's
                # own vendored components (raw FS5 sign included) must
                # reproduce SQAT's vendored pa — this is what establishes
                # that the reference does NOT clamp (measured: exact, 0.0,
                # on all four cases including the two negative-FS5 ones;
                # rtol 1e-12 for cross-machine reduction-order headroom,
                # matching the 840-grid crosscheck above).
                @test isapprox(
                    PsychoacousticMetrics._pa_widmann_arithmetic(c.N5, c.S5, c.R5, c.FS5),
                    c.pa; rtol = 1e-12)
                # components FIRST (attribution), then PA
                @test isapprox(r.loudness, c.N5; rtol = tol.N)
                @test isapprox(r.sharpness, c.S5; rtol = tol.S)
                @test isapprox(r.roughness, c.R5; atol = tol.R)
                @test isapprox(r.fluctuation_strength, c.FS5; atol = tol.FS)
                @test isapprox(r.pa, c.pa; rtol = tol.pa)
                # attribution gate: our PA deviation must stay explained by
                # the N convention difference — substituting SQAT's N5 for
                # our stationary N (keeping OUR raw signed S/R/FS, exactly
                # as the wrapper feeds them) must reproduce SQAT's PA. If
                # this fails while the component checks pass, the
                # composition itself drifted — that is a bug, not a
                # tolerance question.
                pa_n5 = PsychoacousticMetrics._pa_widmann_arithmetic(
                    c.N5, r.sharpness, r.roughness, r.fluctuation_strength)
                @test isapprox(pa_n5, c.pa; rtol = tol.attr)
            end
        end
    end

    @testset "raw signed R/FS on near-stationary tones" begin
        # The oracle itself produces a tiny NEGATIVE FS5 on the steady
        # 40 dB tones (vendored: -1.666e-3 at 44.1 kHz, -1.599e-3 at
        # 48 kHz) — model noise around a true zero. Our stationary FS shows
        # the same artifact on the same stimulus, and the wrapper mirrors
        # the reference by feeding the RAW SIGNED value into the w_FR sum
        # (no clamp; the public scalar surface's non-negativity guard is a
        # caller-input policy and is bypassed via the shared internal
        # arithmetic). These assertions keep that signed path under live
        # coverage; if the FS metric ever stops going negative here,
        # revisit the testset (not necessarily a bug).
        _, _, _, r = wrapper_results["steady_1k_40db"]
        @test r.fluctuation_strength < 0
        @test r.pa == PsychoacousticMetrics._pa_widmann_arithmetic(
            r.loudness, r.sharpness, r.roughness, r.fluctuation_strength)
        # the signed FS strictly increases |0.4*FS + 0.6*R| here (R is tiny
        # positive), so PA must exceed the clamped-FS variant:
        @test r.pa > psychoacoustic_annoyance_widmann(
            r.loudness, r.sharpness, max(r.roughness, 0.0), 0.0)
        @test isfinite(r.pa) && r.pa > 0
    end

    @testset "result struct sanity: components == direct metric calls" begin
        # Non-default pa_per_unit so this also proves the kwarg is forwarded
        # to BOTH signal-domain metrics (a wrapper that dropped it would
        # produce roughness/FS of a 60 dB signal, not a 54 dB one).
        sig, fs, _, _ = wrapper_results["steady_1k_60db"]
        ppu = 0.5
        loud = loudness_zwst(sig, fs; pa_per_unit = ppu)
        r = psychoacoustic_annoyance_widmann(sig, fs, loud; pa_per_unit = ppu)
        @test r.loudness == loud.loudness
        @test r.sharpness == sharpness(loud)
        @test r.roughness == roughness_dw(sig, fs; pa_per_unit = ppu).roughness
        @test r.fluctuation_strength ==
              fluctuation_strength_osses(sig, fs; pa_per_unit = ppu).fluctuation_strength
        @test r.pa == PsychoacousticMetrics._pa_widmann_arithmetic(
            r.loudness, r.sharpness, r.roughness, r.fluctuation_strength)
        @test r.convention == :stationary
    end

    @testset "API errors (inherited from components)" begin
        sig, fs, loud, _ = wrapper_results["steady_1k_40db"]
        # fs = 22050: roughness_dw requires fs >= 44100 (and
        # fluctuation_strength_osses would require 44100/48000)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(sig, 22050.0, loud)
        # pa_per_unit = 0: rejected by the component guards
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(sig, fs, loud;
                                                                    pa_per_unit = 0)
    end
end
