# Thesis-conformance check against Fastl & Zwicker (2007) experimental
# fluctuation-strength curves, transcribed from Osses, García & Kohlrausch
# (2018 PhD thesis), Table B.1 (AM tones / AM broadband noise / FM tones,
# fmod in {1,2,4,8,16,32} Hz).
#
# The stimuli are synthesized with the validation-dataset laws
# (`fs_*_osses2016` in test/support/fs_stimuli.jl), reverse-engineered from
# the 18 files SQAT's own validation scripts read (docs/oracle-pins.md
# §2.13). Three assertions per grid point:
#
# 1. ORACLE PARITY on identical samples: `fluctuation_strength_osses` must
#    reproduce SQAT's `FluctuationStrength_Osses2016` (CC BY-NC,
#    github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb,
#    Octave 11.3.0 / octave-signal 1.4.7) per frame and as the mean, values
#    vendored in test/data/sqat_fs_crosscheck.jl (`thesis_*`). Measured
#    2026-09-02: max relative deviation 4.3e-10 (mean), 2.6e-9 (per frame);
#    rtol 1e-6 leaves > 2.5 orders of margin. Any failure is a transcription
#    bug.
#
# 2. DATASET PARITY: our synthetic stimulus must give the value SQAT gives
#    on the REAL dataset file for that grid point
#    (test/data/sqat_fs_osses2016_dataset.jl). For the deterministic tones
#    this pins the stimulus law itself: measured max 1.0e-3 (AM 32 Hz and
#    FM 1 Hz), asserted at rtol 5e-3. For the noise row our random draw
#    differs from the file's, measured max 11.4 % (32 Hz, a 0.035 vacil
#    value) and <= 5.6 % elsewhere; asserted at rtol 0.25.
#
# 3. PUBLISHED-CURVE GATE where the reference model meets it:
#    `abs(computed - ref)/ref <= 0.30`, i.e. 30 % OF THE REFERENCE value
#    (deliberately not `isapprox(; rtol)`, which scales by max(|a|,|b|) and
#    would grant overestimating points an effective ~43 % bound). Asserted
#    at the 11 points where SQAT itself, on the real dataset files, lands
#    within 30 % of the published value: AM 1/4/8/16 Hz, FM 1/2/4 Hz,
#    noise 2/4/8/16 Hz. An informational 30 %-of-reference tally over all
#    18 points is emitted regardless.
#
# --- Where the reference model misses the published curve ---
#
# Measured on the real dataset files (SQAT = ours to the digits shown):
#
#   point          SQAT     published  ratio
#   AM   2 Hz      1.1062   0.84       1.32x   (shoulder)
#   AM  32 Hz      0.0164   0.06       0.27x   (tail)
#   FM   8 Hz      2.5988   0.70       3.71x   thesis §B.4.1: the model
#   FM  16 Hz      1.2123   0.27       4.49x   overestimates FM tones for
#   FM  32 Hz      0.0545   0.02       2.72x   fmod > 4 Hz, peak shifted to 8 Hz
#   BBN  1 Hz      0.4565   1.12       0.41x
#   BBN 32 Hz      0.0349   0.14       0.25x
#
# These are the reference model's own misses (SQAT's validation figures
# show the same), so a published-curve gate there would assert a property
# the model does not have; assertions 1 and 2 cover those points.
#
# History: until 2026-09-02 this file synthesized the grid with the
# cosine-carrier / cosine-envelope convention of SQAT's shipped reference
# signal, which the dataset files do NOT use (AM: (1 - cos) x sin carrier;
# FM: sin carrier with the opposite deviation sign; noise: (1 - cos) applied
# to POWER, i.e. a |sin(pi fmod t)| amplitude envelope). The noise
# difference alone doubled the fluctuation strength and was misattributed
# to a bandwidth mismatch. See docs/oracle-pins.md §2.13.
#
# Suite totals: this file contributes 0 to the Broken column.

isdefined(@__MODULE__, :fs_am_tone) ||
    include(joinpath(@__DIR__, "support", "fs_stimuli.jl"))
isdefined(@__MODULE__, :SQAT_FS_CASES) ||
    include(joinpath(@__DIR__, "data", "sqat_fs_crosscheck.jl"))
isdefined(@__MODULE__, :SQAT_FS_OSSES2016_DATASET) ||
    include(joinpath(@__DIR__, "data", "sqat_fs_osses2016_dataset.jl"))

# (fmod [Hz], FS [vacil]) — experimental references, Osses 2018 thesis Table
# B.1 (transcribed from Fastl & Zwicker 2007 via the thesis).
const FS_REF_AM_70DB = [(1.0, 0.39), (2.0, 0.84), (4.0, 1.25), (8.0, 1.30), (16.0, 0.36), (32.0, 0.06)]
const FS_REF_FM_70DB = [(1.0, 0.85), (2.0, 1.17), (4.0, 2.00), (8.0, 0.70), (16.0, 0.27), (32.0, 0.02)]
const FS_REF_AMBBN_60DB = [(1.0, 1.12), (2.0, 1.58), (4.0, 1.80), (8.0, 1.57), (16.0, 0.48), (32.0, 0.14)]

# Points where SQAT on the real dataset files lands within 30 % of the
# published value (header); the published-curve gate is asserted only there.
const _FS_PUBLISHED_GATED = Set([
    ("am", 1.0), ("am", 4.0), ("am", 8.0), ("am", 16.0),
    ("fm", 1.0), ("fm", 2.0), ("fm", 4.0),
    ("bbn", 2.0), ("bbn", 4.0), ("bbn", 8.0), ("bbn", 16.0),
])

const _FS_ORACLE_RTOL = 1e-6                     # measured 2.6e-9 per frame, 4.3e-10 mean
const _FS_DATASET_RTOL = Dict("am" => 5e-3, "fm" => 5e-3, "bbn" => 0.25)  # measured 1.0e-3 / 1.0e-3 / 0.114

# Gate (roughness/D&W convention): within 30% OF THE REFERENCE value.
const _GATE_RTOL = 0.30
_fs_gate(computed, ref) = abs(computed - ref) / ref <= _GATE_RTOL

_fs_oracle_case(kind, fmod) =
    only(filter(c -> c.name == "thesis_$(kind)_$(Int(fmod))hz", SQAT_FS_CASES))

@testset "thesis Table B.1 conformance (Fastl & Zwicker via Osses 2018)" begin
    within30 = 0
    total = 0

    for (label, kind, refs) in (
            ("AM tones, 70 dB SPL", "am", FS_REF_AM_70DB),
            ("FM tones, 70 dB SPL, fdev=700 Hz", "fm", FS_REF_FM_70DB),
            ("AM broadband noise, 60 dB SPL", "bbn", FS_REF_AMBBN_60DB),
        )
        @testset "$label" begin
            for (fmod, fs_ref) in refs
                oracle = _fs_oracle_case(kind, fmod)
                sig, fs, _ = synthesize_case(oracle.name)   # the samples the oracle saw
                r = fluctuation_strength_osses(sig, fs)
                total += 1
                within30 += _fs_gate(r.fluctuation_strength, fs_ref)
                @testset "fmod=$(fmod)Hz" begin
                    # 1. oracle parity on identical samples
                    @test length(r.fluctuation_strength_over_time) == length(oracle.fluct)
                    @test isapprox(r.fluctuation_strength_over_time, oracle.fluct; rtol = _FS_ORACLE_RTOL)
                    @test isapprox(r.fluctuation_strength, oracle.fs_mean; rtol = _FS_ORACLE_RTOL)
                    # 2. our stimulus reproduces SQAT's value on the real dataset file
                    @test isapprox(r.fluctuation_strength, SQAT_FS_OSSES2016_DATASET[(kind, fmod)];
                                   rtol = _FS_DATASET_RTOL[kind])
                    # 3. published curve, where the reference model meets it
                    if (kind, fmod) in _FS_PUBLISHED_GATED
                        @test _fs_gate(r.fluctuation_strength, fs_ref)
                    end
                end
            end
        end
    end

    @info "thesis Table B.1: $within30/$total points within 30% of the Fastl & Zwicker reference curves (informational)"
end
