isdefined(@__MODULE__, :fs_am_tone) ||
    include(joinpath(@__DIR__, "support", "fs_stimuli.jl"))

using Statistics: mean

@testset "fluctuation strength pipeline stages" begin
    P = PsychoacousticMetrics

    @testset "_fs_frames matches MATLAB buffer('nodelay') pins" begin
        # Frame counts pinned against Octave (Task 1,
        # .superpowers/sdd/buffer_pins.json, SQAT @ 00b449e, Octave 11.3.0,
        # octave-signal 1.4.7): nframes = floor((L-V)/hop) + 1, NOT
        # ceil((L-V)/hop) — the two formulas disagree whenever L-V is an
        # exact multiple of hop, which is 3 of the 4 pinned cases below.
        fs = 44100
        N = 2 * fs
        V = round(Int, 0.9 * N)
        hop = N - V
        PIN_2s = 2
        PIN_3s = 7
        PIN_3_05s = 7
        PIN_5s = 17
        for (dur, want) in ((2.0, PIN_2s), (3.0, PIN_3s), (3.05, PIN_3_05s), (5.0, PIN_5s))
            L = round(Int, dur * fs)
            frames, starts = P._fs_frames(zeros(L), N)
            @test size(frames) == (N, want)
            @test starts == [1 + (k - 1) * hop for k in 1:want]
        end
        # zero-padded tail: last frame of the 3.05 s case ends in zeros
        L = round(Int, 3.05 * 44100)
        frames, starts = P._fs_frames(ones(L), N)
        @test frames[end, end] == 0.0
        @test frames[1, 1] == 1.0
        # L < V: outside the oracle's reachable domain (octave-signal's
        # buffer returns an empty N×0 matrix or errors internally there —
        # verified against the pinned checkout, see the _fs_frames comment
        # block). Our clamp yields one zero-padded frame instead of an
        # opaque allocation error; deliberate graceful extension.
        frames, starts = P._fs_frames(ones(4410), N)
        @test size(frames) == (N, 1)
        @test starts == [1]
        @test all(frames[1:4410, 1] .== 1.0)
        @test all(frames[4411:end, 1] .== 0.0)
    end

    @testset "_cos_ramp_window quirks (cos_ramp @ SQAT 00b449e)" begin
        fs = 44100
        N = 2 * fs
        w = P._cos_ramp_window(N, fs, 50.0)
        a = round(Int, fs * 50 / 1000)
        @test w[1] == 0.5 * (1 - cos(π * 1 / a))   # starts ABOVE zero (1-based)
        @test w[a] == 1.0                           # attack ends exactly at 1
        @test w[N] == 0.0                           # release ends exactly at 0
        @test w[N - a] == 1.0                       # release starts at 1 (cos(-π))
        @test all(w[a+1:N-a-1] .== 1.0)             # flat middle
    end

    @testset "_apply_a0_fir delay compensation and gain" begin
        fs = 44100.0
        t = (0:2*44100-1) ./ fs
        x = sin.(2π * 1000.0 .* t)
        y = P._apply_a0_fir(x, P._A0_FIR_B_44100)
        @test length(y) == length(x)
        seg = 22050:66150                            # avoid edge transients
        @test isapprox(sqrt(mean(abs2, y[seg])) / sqrt(mean(abs2, x[seg])), 1.0; atol = 0.02)
        # zero-delay: cross-correlation peak at lag 0 (linear phase compensated)
        @test abs(sum(x[seg] .* y[seg])) > 0.98 * sum(abs2, x[seg])
    end

    @testset "_terhardt_excitation_fs anchor frame vs oracle dump" begin
        # h0 per channel for anchor_44k frame 1, dumped from the Octave
        # oracle (Task 1, .superpowers/sdd/stage_dumps.json, provenance:
        # SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb, Octave 11.3.0,
        # octave-signal 1.4.7, dump_fs_stage.m). h0 = mean(abs(ei)) is
        # computed in Task 4's stage, so compare it here from the raw ei.
        sig = fs_am_tone(1000.0, 4.0, 60.0, 4.0, 44100.0)
        N = 2 * 44100
        frames, _ = P._fs_frames(sig, N)
        w = P._cos_ramp_window(N, 44100.0, 50.0)
        frame = P._apply_a0_fir(w .* frames[:, 1], P._A0_FIR_B_44100)
        ei = P._terhardt_excitation_fs(frame, 44100.0)
        @test size(ei) == (47, N)
        h0 = vec(mean(abs.(ei); dims = 2))
        ORACLE_H0 = Float64[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0.000032159114674729847, 0.00015326981202008495,
            0.0007217895976721765, 0.0034014301169236298,
            0.0033997404716480026, 0.0033999596234046264,
            0.0016406547878079087, 0.0007912871893871107,
            0.0003817410929765959, 0.00018407520060818755,
            0.00008878934636662715, 0.00004282815719867357,
            0.000020658573525663808, 0.000009963665296158618,
            0.000004806195621909302,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ]
        @test isapprox(h0, ORACLE_H0; rtol = 1e-6, atol = 1e-12)
        # silence in, silence out
        @test all(P._terhardt_excitation_fs(zeros(N), 44100.0) .== 0.0)
    end
end
