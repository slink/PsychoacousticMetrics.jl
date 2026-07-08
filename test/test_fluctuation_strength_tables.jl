function roots_of_quadratic(a, b, c)
    d = sqrt(complex(b^2 - 4a*c))
    return ((-b + d) / 2a, (-b - d) / 2a)
end

@testset "fluctuation strength tables" begin
    P = PsychoacousticMetrics

    @testset "gzi weighting (FS variant)" begin
        # knots from il_Get_gzi_fluctuation @ SQAT 00b449e
        @test P._gzi_fluctuation(0.5) == 1.0
        @test P._gzi_fluctuation(13.0) == 1.0
        @test P._gzi_fluctuation(15.0) == 0.9
        @test P._gzi_fluctuation(17.5) == 0.7
        @test P._gzi_fluctuation(23.5) ≈ 0.7 + (23.5 - 17.5) / (24 - 17.5) * (0.5 - 0.7)
        @test P._gzi_fluctuation(24.5) == 0.5     # beyond last knot: reference maps NaN→0.5
    end

    @testset "a0 simplified curve (FS variant)" begin
        # knots from calculate_a0 @ SQAT 00b449e, a0_type='fluctuationstrength_osses2016'
        @test P._a0_fs_db(0.0) == 0.0
        @test P._a0_fs_db(19.0) == 0.0            # flat through 19 Bark (no ear-canal resonance)
        @test P._a0_fs_db(20.0) == -1.43
        @test P._a0_fs_db(23.5) == -20.0
        @test P._a0_fs_db(25.0) == -130.0
    end

    @testset "Hweight SOS shape and stability" begin
        for sos in (P._HWEIGHT_SOS_44100, P._HWEIGHT_SOS_48000)
            @test size(sos, 2) == 6
            @test size(sos, 1) == 3                    # 1 HP section + 2 LP sections
            @test all(sos[:, 4] .== 1.0)               # normalized a0
            # poles inside unit circle
            for s in axes(sos, 1)
                r = roots_of_quadratic(sos[s, 4], sos[s, 5], sos[s, 6])
                @test all(abs.(r) .< 1.0)
            end
        end
    end

    @testset "_sosfilt matches oracle dump" begin
        # Task 1 dumped sosfilt(Hweight_44100, probe) for a fixed probe vector;
        # vendored alongside the SOS constants (small, 64 samples).
        y = P._sosfilt(P._HWEIGHT_SOS_44100, P._FS_SOSFILT_PROBE_IN)
        @test isapprox(y, P._FS_SOSFILT_PROBE_OUT; rtol = 1e-12, atol = 1e-15)
    end

    @testset "_sosfilt passband/stopband physics" begin
        fs = 44100.0
        t = (0:4*44100-1) ./ fs
        pass = P._sosfilt(P._HWEIGHT_SOS_44100, sin.(2π * 4.0 .* t))
        stop = P._sosfilt(P._HWEIGHT_SOS_44100, sin.(2π * 0.1 .* t))
        seg = 2*44100:4*44100-1                     # skip transients
        @test maximum(abs, pass[seg]) > 0.7          # 4 Hz in passband
        @test maximum(abs, stop[seg]) < 0.1          # 0.1 Hz rejected
    end

    @testset "a0 FIR vendored coefficients" begin
        @test length(P._A0_FIR_B_44100) == 4097
        @test length(P._A0_FIR_B_48000) == 4097
        # linear phase: symmetric taps
        @test isapprox(P._A0_FIR_B_44100, reverse(P._A0_FIR_B_44100); atol = 1e-12)
        # DC region: |H| ≈ 1 below 19 Bark. Evaluate at 1 kHz via DFT of taps.
        H(B, f, fs) = abs(sum(B .* cis.(-2π * f / fs .* (0:length(B)-1))))
        @test isapprox(H(P._A0_FIR_B_44100, 1000.0, 44100.0), 1.0; atol = 0.02)
        @test H(P._A0_FIR_B_44100, 15000.0, 44100.0) < 0.2   # deep in the roll-off
    end
end
