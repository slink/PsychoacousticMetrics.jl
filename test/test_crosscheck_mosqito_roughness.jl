# Transcription cross-check: MoSQITo's roughness_dw run on identical signals.
#
# Tolerance procedure (from the task brief): start rtol = 1e-6; if that fails
# but the deviation is < 1e-3, set rtol to 10x the measured max relative
# deviation and record the number; > 1e-3 indicates a transcription bug.
# Measured max relative deviations (2026-07-07, this machine, FFTW vs
# numpy.fft):
#
#   anchor_1k_70    2.237e-6   -> rtol 2.3e-5 (10x measured, per procedure)
#   mislabeled_fs   1.115e-14  -> rtol 2.3e-5 (well inside)
#   fc250_fm40      1.449e-3   -> rtol 1e-2 (adjudicated, see below)
#   fc2000_fm100    1.646e-3   -> rtol 1e-2 (adjudicated, see below)
#
# The two > 1e-3 cases were probed and are NOT transcription bugs (controller
# adjudication of the > 1e-3 rule). Mechanism, pinned per-band on the worst
# frame of each case: one hBP channel is analytically zero but numerically
# pure FFT round-off dust (fc2000 frame 2 band 37: maxabs(hBP) = 2.5e-22
# FFTW vs 5.3e-22 numpy, ~17 orders below that channel's h0 = 1.4e-5;
# fc250 frame 4 band 14: 9.8e-22 vs 1.6e-21). The ki guard's exact-zero
# predicate passes the dust on BOTH sides (exact zeros in hBP: 0/9600 in
# both implementations — no guard flip), so ki becomes the correlation of a
# real envelope with implementation-dependent rounding noise: fc2000
# ki[35] = -0.2415 (Julia/FFTW) vs -0.0017 (MoSQITo/numpy); fc250
# ki[12] = -0.0195 vs +0.0292. Every deterministic quantity in those frames
# (h0, mdepth, neighboring ki) matches to ~1e-15. The deviation is inherent
# FP-dust sensitivity of upstream's guard design, not a transcription error;
# rtol 1e-2 (~6x the measured deviation) covers it while still catching any
# real transcription regression.
#
# fc2000_fm100 verdict for the fig-3 @test_broken cluster: MoSQITo itself
# computes R = 0.80018 asper here vs the Zwicker & Fastl reference 0.63256
# (|diff| = 0.1676, far beyond the 0.1 gate) — it fails its own validation
# gate at this point just like we do. Our max absolute difference vs MoSQITo
# across all frames is 0.00132 asper, smaller than the smallest fig-3 failure
# overshoot (0.00223 asper at fc=2000/fmod=80). The fig-3 fc=2000 failures
# are upstream model behavior, not our transcription.

include(joinpath(@__DIR__, "data", "mosqito_roughness_crosscheck.jl"))

@testset "MoSQITo roughness cross-check" begin
    fs_synth = 48000
    params = Dict(
        "anchor_1k_70"  => (1000.0, 70.0, 1.0),
        "mislabeled_fs" => (1000.0, 70.0, 1.0),
        "fc250_fm40"    => (250.0, 40.0, 1.5),
        "fc2000_fm100"  => (2000.0, 100.0, 1.5),
    )
    # Per-case rtol per the header: 2.3e-5 by the brief's 10x-measured rule;
    # 1e-2 for the two adjudicated FP-dust cases.
    rtols = Dict(
        "anchor_1k_70"  => 2.3e-5,
        "mislabeled_fs" => 2.3e-5,
        "fc250_fm40"    => 1e-2,
        "fc2000_fm100"  => 1e-2,
    )
    for (name, fs_analyze, overlap, R_expected) in MOSQITO_ROUGHNESS_CROSSCHECK
        fc, fm, dur = params[name]
        t = range(0, dur; length = Int(round(dur * fs_synth)))
        stim = am_sine(sin.(2π * fm .* t), fs_synth, fc, 60)
        r = roughness_dw(stim, fs_analyze; overlap)
        @testset "$name" begin
            @test length(r.roughness_over_time) == length(R_expected)
            @test all(isapprox.(r.roughness_over_time, R_expected; rtol = rtols[name]))
        end
    end
end
