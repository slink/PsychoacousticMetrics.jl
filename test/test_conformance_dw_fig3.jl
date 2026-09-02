# Daniel & Weber figure-3 conformance against the Zwicker & Fastl reference
# curves — exactly MoSQITo's validation gate
# (validation_roughness_danielweber.py: |R - ref_zf| <= 0.1 asper, AM tones
# 1.5 s / 48 kHz / 60 dB / m=1 / overlap 0, first frame). The script's
# 30 %-relative comparison against the D&W curves is informational upstream
# and is not asserted here. Nothing outside this grid is asserted.
#
# KNOWN MODEL DEVIATION (documented, gated on the reference implementation
# instead): the following 5/77 grid points, all at fc=2000 Hz across the
# mid/high fmod band, exceed the ±0.1 asper Zwicker & Fastl gate. MoSQITo
# itself (the model this package transcribes) misses the same gate at the
# same five points by the same amount, and also overshoots Daniel & Weber's
# own published curve there by 15–31 % — so this is upstream model behavior
# at 2 kHz, not a transcription error. Values computed with this exact test
# setup (fs=48000, 1.5 s AM tone, 60 dB SPL, m=1, overlap=0, first frame),
# MoSQITo @ d990c33f94f1 / mosqito==1.2.1 via
# scripts/generate_mosqito_roughness_crosscheck.jl:
#
#   fc [Hz]  fmod [Hz]  R (ours)     R (MoSQITo)  R_ref (Z&F)  |ours-ZF|  |MoSQITo-ZF|  |ours-MoSQITo|
#   2000      80        0.913193     0.912519     0.810968     0.102      0.102         6.7e-4
#   2000      90        0.867606     0.867606     0.740062     0.128      0.128         1.8e-8
#   2000     100        0.800337     0.800179     0.632556     0.168      0.168         1.6e-4
#   2000     120        0.642240     0.643844     0.469837     0.172      0.174         1.6e-3
#   2000     140        0.505632     0.505441     0.361517     0.144      0.144         1.9e-4
#
# For these five points the assertion is therefore agreement with MoSQITo's
# first-frame value (vendored in test/data/mosqito_roughness_crosscheck.jl)
# within 5e-3 asper (~3x the measured maximum, and 20x below the Z&F gate),
# plus a bound that the Z&F miss stays within 0.2 asper so a drift in
# either direction is still caught. MoSQITo's own test suite never asserts
# fig. 3 at all (it gates only the 1 kHz / 70 Hz anchor at ±17 %); its
# validation script computes the ±0.1 compliance flag but only plots it.
#
# Nearest passing point to the gate: fc=500 Hz, fmod=40 Hz, margin +0.00931
# (R=0.6742020504829194 vs R_ref=0.5835151601969246). No other passing point
# has margin < 0.01. All five 2 kHz margins and the (500, 40) margin were
# re-measured across two independent `julia --project=.` processes, both with
# FFTW pinned to 1 thread (see below) and with FFTW's default thread count on
# this machine (also 1) — every value was bit-identical across runs
# (measured max cross-run |ΔR| over all 77 grid points = 0). This is expected:
# `fft`/`ifft` here use FFTW's default ESTIMATE-mode planning, which does not
# do runtime auto-tuning, so it should not vary run to run regardless of
# thread count. An earlier report of ~2e-4 cross-process variation at
# fc=2000/fmod=80 was NOT reproduced in this measurement — it is recorded
# here only as a hypothesis (possibly environment-specific FFT planning
# nondeterminism, e.g. a different BLAS/FFTW build or thread count on another
# machine), not a confirmed cause. `FFTW.set_num_threads(1)` below is kept as
# a cheap, harmless determinism pin for gate-boundary points regardless.
#
# All other 72/77 grid points, and all fc != 2000 Hz points, pass. The ±0.1
# gate itself is never adjusted.

import FFTW
FFTW.set_num_threads(1)  # pin run-to-run determinism for gate-boundary points

# This file uses `am_sine`, normally defined by test_roughness_dw.jl earlier
# in the runtests.jl include order; guard against running this file alone.
isdefined(@__MODULE__, :am_sine) ||
    include(joinpath(@__DIR__, "support", "am_generator.jl"))

include(joinpath(@__DIR__, "data", "dw_fig3_references.jl"))
isdefined(@__MODULE__, :MOSQITO_ROUGHNESS_CROSSCHECK) ||
    include(joinpath(@__DIR__, "data", "mosqito_roughness_crosscheck.jl"))

# (fc, fmod) => MoSQITo first-frame roughness for the five known 2 kHz
# deviation points (overlap 0, same synthesis as this test).
const _DW_FIG3_MOSQITO_2K = Dict(
    (2000, fm) => only(filter(c -> c[1] == "fc2000_fm$(fm)", MOSQITO_ROUGHNESS_CROSSCHECK))[4][1]
    for fm in (80, 90, 100, 120, 140)
)

@testset "D&W fig.3 conformance (Zwicker & Fastl, ±0.1 asper)" begin
    fs = 48000
    t = range(0, 1.5; length = Int(1.5 * fs))    # python linspace endpoint-inclusive
    within30 = 0                                 # informational D&W comparison
    for (i, (fc, fmod, R_ref)) in enumerate(DW_FIG3_REF_ZF)
        stim = am_sine(sin.(2π * fmod .* t), fs, fc, 60)
        r = roughness_dw(stim, fs; overlap = 0.0)
        R = r.roughness_over_time[1]
        @testset "fc=$(fc)Hz fmod=$(fmod)Hz" begin
            if haskey(_DW_FIG3_MOSQITO_2K, (fc, fmod))
                # Known upstream deviation (header): gate on the reference
                # implementation, and bound the Z&F miss.
                @test abs(R - _DW_FIG3_MOSQITO_2K[(fc, fmod)]) <= 5e-3
                @test abs(R - R_ref) <= 0.2
            else
                @test abs(R - R_ref) <= 0.1
            end
        end
        R_dw_ref = DW_FIG3_REF_DW[i][3]
        within30 += (abs(R - R_dw_ref) / R_dw_ref) <= 0.30
    end
    # Upstream reports (but never gates) agreement with the Daniel & Weber
    # curves within 30 % relative; we do the same.
    @info "D&W fig.3: $within30/$(length(DW_FIG3_REF_ZF)) points within 30% of the Daniel & Weber curves (informational)"
end
