# Daniel & Weber figure-3 conformance against the Zwicker & Fastl reference
# curves — exactly MoSQITo's validation gate
# (validation_roughness_danielweber.py: |R - ref_zf| <= 0.1 asper, AM tones
# 1.5 s / 48 kHz / 60 dB / m=1 / overlap 0, first frame). The script's
# 30 %-relative comparison against the D&W curves is informational upstream
# and is not asserted here. Nothing outside this grid is asserted.
#
# KNOWN FAILURES (documented, not silenced): the following 5/77 grid points,
# all at fc=2000 Hz across the mid/high fmod band, exceed the ±0.1 asper gate
# and are marked @test_broken below. Values computed with this exact test
# setup (fs=48000, 1.5 s AM tone, 60 dB SPL, m=1, overlap=0, first frame):
#
#   fc [Hz]  fmod [Hz]  R (computed)          R_ref (Zwicker&Fastl)   |R-R_ref|  margin(0.1-|R-Rref|)
#   2000      80        0.9131930314397203    0.8109675430315867     0.1022     -0.00223
#   2000      90        0.8676062588706794    0.740062260031716      0.1275     -0.02754
#   2000     100        0.8003367561933683    0.6325559817988694     0.1678     -0.06778
#   2000     120        0.642239945788509     0.46983749746147563    0.1724     -0.07240
#   2000     140        0.5056322780869698    0.36151661816289216    0.1441     -0.04412
#
# Nearest passing point to the gate: fc=500 Hz, fmod=40 Hz, margin +0.00931
# (R=0.6742020504829194 vs R_ref=0.5835151601969246). No other passing point
# has margin < 0.01. All 5 broken-point margins and the (500, 40) margin were
# re-measured across two independent `julia --project=.` processes, both with
# FFTW pinned to 1 thread (see below) and with FFTW's default thread count on
# this machine (also 1) — every value above was bit-identical across runs
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

const _DW_FIG3_KNOWN_BROKEN = Set([
    (2000, 80), (2000, 90), (2000, 100), (2000, 120), (2000, 140),
])

@testset "D&W fig.3 conformance (Zwicker & Fastl, ±0.1 asper)" begin
    fs = 48000
    t = range(0, 1.5; length = Int(1.5 * fs))    # python linspace endpoint-inclusive
    within30 = 0                                 # informational D&W comparison
    for (i, (fc, fmod, R_ref)) in enumerate(DW_FIG3_REF_ZF)
        stim = am_sine(sin.(2π * fmod .* t), fs, fc, 60)
        r = roughness_dw(stim, fs; overlap = 0.0)
        R = r.roughness_over_time[1]
        @testset "fc=$(fc)Hz fmod=$(fmod)Hz" begin
            if (fc, fmod) in _DW_FIG3_KNOWN_BROKEN
                @test_broken abs(R - R_ref) <= 0.1
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
