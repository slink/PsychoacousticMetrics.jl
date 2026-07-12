isdefined(@__MODULE__, :fs_am_tone) ||
    include(joinpath(@__DIR__, "support", "fs_stimuli.jl"))
isdefined(@__MODULE__, :SQAT_FS_CASES) ||
    include(joinpath(@__DIR__, "data", "sqat_fs_crosscheck.jl"))

using Statistics: mean

@testset "fluctuation_strength_osses" begin
    P = PsychoacousticMetrics

    @testset "anchor gate" begin
        # 1 kHz, 60 dB, 100% AM at 4 Hz -> 1 vacil (thesis App. B.3; SQAT
        # docstring). Tolerance derivation (measured 2026-07-08, see
        # .superpowers/sdd/fs-task-4-report.md for the full log):
        #   |oracle(anchor_44k fs_mean) - 1.0|  = |1.0053626323366578 - 1| = 5.3626e-3
        #     (this is the reference model's own distance from the ideal 1
        #     vacil, measured once from the fixture — not derived from the
        #     value under test, so the gate below is falsifiable)
        #   |Julia FS(anchor) - oracle fs_mean|  = 1.0474e-10 (measured
        #     separately; this package's own oracle rig gives
        #     near-machine-precision parity, consistent with the 1.3e-13
        #     parity measured at the excitation stage in Task 3)
        # ANCHOR_ATOL is set to 2x oracle_dev (2 * 5.3626e-3 = 1.0725e-2):
        # documented, not a silently-widened magic number, and independent
        # of r.fluctuation_strength so a regression in the code under test
        # cannot inflate its own tolerance. The tight rtol=1e-6 check below
        # (against oracle_fs_mean, not 1.0) is the gate that actually
        # catches Julia-vs-oracle drift.
        sig = fs_am_tone(1000.0, 4.0, 60.0, 4.0, 44100.0)
        r = fluctuation_strength_osses(sig, 44100)
        oracle_fs_mean = 1.0053626323366578
        oracle_dev = abs(oracle_fs_mean - 1.0)
        ANCHOR_ATOL = 2 * oracle_dev
        @test isapprox(r.fluctuation_strength, 1.0; atol = ANCHOR_ATOL)
        # tight self-consistency check vs the oracle fixture (not vs 1.0)
        @test isapprox(r.fluctuation_strength, oracle_fs_mean; rtol = 1e-6)
    end

    @testset "stage dumps: mdepth, ki, specific fi (anchor frame 1)" begin
        # Literals from .superpowers/sdd/stage_dumps.json (provenance: SQAT @
        # 00b449e40599f1c1ef4abe0596094552213d57eb, Octave 11.3.0,
        # octave-signal 1.4.7, dump_fs_stage.m, untracked scratch instrumentation).
        # All 47 channels are compared — no channel exclusions.
        #
        # Tolerances (measured 2026-07-08, both with and without
        # --check-bounds=yes, which is Pkg.test's default and changes float
        # reduction/SIMD order): mdepth/ki1/ki2 worst relative deviation is
        # 6.1e-7 (ki1, channel 25 — its correlation partner channel 27 has
        # h0 = 4.8e-6, ~3 orders below peak, so summation-order changes
        # shift it at the 1e-7 level) -> rtol 1e-6 holds. fi amplifies the
        # ki deviation by the model's ^1.7 exponent (1.7 x 6.1e-7 =
        # 1.04e-6, measured exactly), so an rtol of 1e-6 on fi is
        # mathematically inconsistent with allowing 1e-6 on ki; fi uses
        # rtol 3e-6 (~3x headroom over the measured 1.04e-6 worst case).
        sig = fs_am_tone(1000.0, 4.0, 60.0, 4.0, 44100.0)
        N = 2 * 44100
        frames, _ = P._fs_frames(sig, N)
        w = P._cos_ramp_window(N, 44100.0, 50.0)
        frame = P._apply_a0_fir(w .* frames[:, 1], P._A0_FIR_B_44100)
        ei = P._terhardt_excitation_fs(frame, 44100.0)
        mdepth, hBP = P._fs_modulation_depths(ei, P._HWEIGHT_SOS_44100)
        ki1, ki2 = P._fs_cross_covariance(hBP)
        gzi = P._gzi_fluctuation.((1:47) ./ 2)
        fi = 0.4980 .* P._fs_specific(mdepth, ki1, ki2, gzi)   # cal @ 00b449e (matches main function)

        # NOTE: P._fs_modulation_depths returns the PRE-compression modulation
        # depth (compression to slope 0.3 above 0.7 happens inside
        # _fs_specific), so this compares against stage_dumps.json's
        # `mdepth_pre` field, not `mdepth_post`. The compression formula
        # itself is exercised indirectly by the `fi` comparison below, since
        # `_fs_specific` applies it before raising to the 1.7 power.
        ORACLE_MDEPTH_PRE = Float64[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0.6607249573264759, 0.674427409998579, 0.706708135045457,
            0.6828180146947915, 0.710343508747022, 0.6956590774617057,
            0.6479293931750263, 0.593921226344826, 0.5409168114104557,
            0.493107811632327, 0.45820840027065909, 0.42577374256203645,
            0.3956312306547457, 0.000022674477188234647, 0.00002267447718813283,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ]
        ORACLE_KI1 = Float64[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0.9777848431681683, 0.9896992696797889, 0.9999553360574617,
            0.9999970529899658, 0.9991090174050349, 0.9972178482863879,
            0.9937315248598249, 0.9903387090779031, 0.9940710236917404,
            0.999999998908596, 0.9999999991320973, 0.101240501817624,
            0.10124527951540761,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ]
        ORACLE_KI2 = Float64[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0.9777848431681683, 0.9896992696797889, 0.9999553360574617,
            0.9999970529899658, 0.9991090174050349, 0.9972178482863879,
            0.9937315248598249, 0.9903387090779031, 0.9940710236917404,
            0.999999998908596, 0.9999999991320973, 0.101240501817624,
            0.10124527951540761,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ]
        ORACLE_FI = Float64[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0.2626613381323144, 0.25580035420900257, 0.27319360706836529,
            0.2674509938085147, 0.23524787588286809, 0.20106840750509528,
            0.1716003531240847, 0.14725258412347864, 0.13081156051485433,
            0.0023764885218376826, 0.002097782403948298,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ]
        for i in 1:47
            @test isapprox(mdepth[i], ORACLE_MDEPTH_PRE[i]; rtol = 1e-6, atol = 1e-9)
            @test isapprox(ki1[i], ORACLE_KI1[i]; rtol = 1e-6, atol = 1e-9)
            @test isapprox(ki2[i], ORACLE_KI2[i]; rtol = 1e-6, atol = 1e-9)
            @test isapprox(fi[i], ORACLE_FI[i]; rtol = 3e-6, atol = 1e-9)  # see tolerance note above
        end
    end

    @testset "API behavior" begin
        sig = fs_am_tone(1000.0, 4.0, 60.0, 4.0, 44100.0)
        @test_throws ArgumentError fluctuation_strength_osses(sig, 22050)
        @test_throws ArgumentError fluctuation_strength_osses(sig, 44100; pa_per_unit = 0.0)
        @test_throws ArgumentError fluctuation_strength_osses(sig, 44100; method = :nope)
        r = fluctuation_strength_osses(zeros(3 * 44100), 44100)
        @test r.fluctuation_strength == 0.0                    # silence, never NaN
        @test !any(isnan, r.specific_fluctuation_strength)

        # method = :stationary on the 4 s anchor: N == length(sig), so
        # L - V == hop EXACTLY (pinned quirk, fs-oracle-pins.md Step 3.1 /
        # SQAT_FS_CASES "stationary_anchor") -> 2 frames, NOT 1.
        rs = fluctuation_strength_osses(sig, 44100; method = :stationary)
        @test length(rs.fluctuation_strength_over_time) == 2

        # signal shorter than 2 s falls back to :stationary with a warning;
        # same pinned quirk applies (SQAT_FS_CASES "short_fallback" -> 2 frames)
        short = sig[1:round(Int, 1.5 * 44100)]
        rw = @test_logs (:warn, r"stationary") match_mode=:any fluctuation_strength_osses(short, 44100)
        @test length(rw.fluctuation_strength_over_time) == 2

        @test r.bark_axis == collect(0.5:0.5:23.5)
        @test size(r.specific_fluctuation_strength, 1) == 47
        # time axis = frame START times (SQAT t_b(1,:) — differs from
        # roughness's frame centers)
        hop = 2 * 44100 - round(Int, 0.9 * 2 * 44100)
        @test r.time_axis[1] == 1 / 44100
        @test r.time_axis[2] - r.time_axis[1] ≈ hop / 44100
    end

    @testset "band-pass character" begin
        # FS(4 Hz) must exceed FS(0.5 Hz) and FS(32 Hz) for AM tones
        f(fm) = fluctuation_strength_osses(
            fs_am_tone(1000.0, fm, 70.0, 4.0, 44100.0), 44100).fluctuation_strength
        @test f(4.0) > f(0.5)
        @test f(4.0) > f(32.0)
    end

    @testset "24 Bark behavior" begin
        # Pinned fact (fs-oracle-pins.md Step 3.3): a 15.8 kHz, 70 dB tone
        # (>= 24 Bark) yields FSmean = 0.0 EXACTLY, silently — not an error,
        # not NaN. Fixture: SQAT_FS_CASES "tone_25bark" (fs_mean = 0.0).
        # This is the INAUDIBLE regime: the a0 roll-off pushes the level
        # below the hearing threshold before the excitation stage even sees
        # it, so `audible` is empty and the short-circuit-to-zero path
        # (never the >= 24 Bark guard below) is what fires.
        sig, fs, _ = synthesize_case("tone_25bark")
        r = fluctuation_strength_osses(sig, fs)
        @test r.fluctuation_strength == 0.0
        @test !any(isnan, r.specific_fluctuation_strength)

        # AUDIBLE regime: raise the level so the component survives the a0
        # roll-off and is still above threshold at >= 24 Bark. Measured
        # (2026-07-11): a 15.5 kHz, 100 dB SPL, 4 s pure tone at 44100 Hz
        # has its highest audible bin at exactly bark = 24.0 (verified via
        # the same _terhardt_excitation_fs audible-bin computation this
        # guard uses) — i.e. it DOES trigger the audible->=24-Bark path,
        # unlike tone_25bark above. Upstream (TerhardtExcitationPatterns.m
        # @ SQAT 00b449e) crashes on this exact input: its kk1 expansion
        # (line 68) indexes past the 47-column array once
        # floor(2*bark) >= 48, i.e. bark >= 24 (verified against the
        # pinned Octave oracle). We fail loudly instead of returning an
        # unvalidated number (the unguarded code silently returned
        # 0.2070074124444851 here).
        t = (0:round(Int, 4.0 * 44100)-1) ./ 44100
        p = 2e-5 * 10.0^(100 / 20)
        loud_sig = (p * sqrt(2)) .* cos.(2π * 15500 .* t)
        @test_throws DomainError fluctuation_strength_osses(loud_sig, 44100)

        # Boundary: bark in (23.5, 24), audible, must NOT throw (upstream's
        # `if i ~= 47` branch protects this range). Measured: a 15.4 kHz,
        # 100 dB SPL tone has its highest audible bin at bark = 23.975 —
        # audible, just under 24 — and completes normally.
        boundary_sig = (p * sqrt(2)) .* cos.(2π * 15400 .* t)
        rb = fluctuation_strength_osses(boundary_sig, 44100)
        @test isfinite(rb.fluctuation_strength)
    end
end
