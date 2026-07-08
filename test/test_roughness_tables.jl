@testset "roughness tables" begin
    P = PsychoacousticMetrics

    @testset "threshold in quiet (roughness variant)" begin
        # exact knots from MoSQITo utils/LTQ.py, reference="roughness"
        @test P._ltq_roughness(0.0) == 130.0
        @test P._ltq_roughness(13.3) == 0.0
        @test P._ltq_roughness(16.0) == -4.0
        @test P._ltq_roughness(25.0) == 130.0
    end

    @testset "a0 ear transmission [dB]" begin
        # knots from _ear_filter_coeff.py
        @test P._a0_ear_db(0.0) == 0.0
        @test P._a0_ear_db(10.0) == 0.0
        @test P._a0_ear_db(16.5) == 7.38
        @test P._a0_ear_db(24.0) == -40.0
    end

    @testset "gzi Aures weighting" begin
        # knots from _gzi_weighting.py (x = 0:24)
        @test P._gzi(9.0) == 0.98
        @test P._gzi(0.0) == 0.15
        @test P._gzi(24.0) == 0.30
        @test P._gzi(0.5) ≈ (0.15 + 0.26) / 2
    end

    @testset "H weighting matrix" begin
        n, fs = 9600, 48000
        H = P._h_weighting(n, fs)
        @test size(H) == (47, n)
        # channel groups are exact copies (upstream derivation):
        @test H[1, :] == H[2, :] == H[3, :] == H[4, :]
        @test all(H[i, :] == H[5, :] for i in 6:15)
        @test all(H[i, :] == H[16, :] for i in 17:20)
        @test all(H[i, :] == H[21, :] for i in 22:41)
        @test all(H[i, :] == H[42, :] for i in 43:47)
        # H2 peaks at the 30 Hz bin (5 Hz grid straddles the 32 Hz knot):
        # _lininterp between the (25, 0.975) and (32, 1.0) knots, exactly
        @test maximum(H[2, :]) == 0.975 + (30 - 25) / (32 - 25) * (1.0 - 0.975)
        @test argmax(H[2, :]) == 7  # 1-based bin 7 = 30 Hz
        # zero at DC bins (j < cut)
        @test H[2, 1] == 0.0 && H[2, 2] == 0.0
        # upstream quirk: H16/H21/H42 are truncated at the H5 bound
        # (last = floor(502/fs*n)), so bins above 502 Hz are zero
        last_bin = floor(Int, 502 / fs * n)
        @test all(H[16, last_bin + 2:end] .== 0.0)
    end
end
