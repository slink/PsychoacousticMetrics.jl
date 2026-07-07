# Unit tests for the sharpness API surface. DIN 45692 conformance is covered
# separately in test_conformance_din45692.jl.

@testset "sharpness API" begin
    # README example spectrum of ZwickerLoudness.jl
    spl = [60, 62, 65, 68, 70, 72, 74, 75, 73, 71,
           69, 67, 65, 63, 61, 59, 57, 55, 53, 50,
           47, 44, 41, 38, 35, 32, 29, 26.0]
    r = zwicker_loudness(spl)

    @testset "struct and vector methods agree exactly" begin
        @test sharpness(r) == sharpness(r.specific_loudness)
        @test sharpness(r) > 0.0
    end

    @testset "input validation" begin
        @test_throws ArgumentError sharpness(zeros(239))
        @test_throws ArgumentError sharpness(zeros(241))
        @test_throws ArgumentError sharpness(r.specific_loudness; weighting=:zwicker)
        # weighting is validated even for silent input
        @test_throws ArgumentError sharpness(zeros(240); weighting=:zwicker)
        bad = copy(r.specific_loudness)
        bad[1] = -0.1
        @test_throws ArgumentError sharpness(bad)
    end

    @testset "near-silence returns zero" begin
        @test sharpness(zeros(240)) == 0.0
        # N = 0.1 * sum = 0.09 sone < 0.1 sone threshold
        @test sharpness(fill(0.9 / 240, 240)) == 0.0
    end

    @testset "DIN weighting function transcription" begin
        # DIN 45692 via MoSQITo sharpness_din_from_loudness.py:124-125:
        # g = 1 for z <= 15.8; 0.15*exp(0.42*(z-15.8)) + 0.85 above.
        @test PsychoacousticMetrics._g_din(0.1) == 1.0
        @test PsychoacousticMetrics._g_din(15.8) == 1.0
        @test PsychoacousticMetrics._g_din(20.0) ≈ 0.15 * exp(0.42 * (20.0 - 15.8)) + 0.85
    end

    @testset "generic real input" begin
        @test sharpness(big.(r.specific_loudness)) isa Real
        @test sharpness(big.(r.specific_loudness)) ≈ sharpness(r.specific_loudness) rtol = 1e-12
        @test sharpness(Float32.(r.specific_loudness)) isa Real
    end
end
