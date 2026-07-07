# Bark/level conversion unit tests. Tables from MoSQITo (Zwicker & Fastl
# table 6.1); expected values are exact table entries or hand-computed
# linear interpolations between adjacent entries.

@testset "bark and level conversions" begin
    P = PsychoacousticMetrics

    @testset "freq2bark table anchors" begin
        @test P.freq2bark(0.0) == 0.0
        @test P.freq2bark(1000.0) == 8.5     # table: 1000 Hz ↔ 8.5 Bark
        @test P.freq2bark(20000.0) == 24.5
        @test P.freq2bark(75.0) ≈ 0.75       # midpoint of 50→100 Hz = 0.5→1.0 Bark
        @test P.freq2bark(25000.0) == 24.5   # clamped above
    end

    @testset "bark2freq table anchors" begin
        @test P.bark2freq(8.5) == 1000.0
        @test P.bark2freq(0.0) == 0.0
        @test P.bark2freq(24.5) == 20000.0
        @test P.bark2freq(0.75) ≈ 75.0
        @test P.bark2freq(30.0) == 20000.0   # clamped above
    end

    @testset "amp2db / db2amp" begin
        @test P._amp2db(2e-5; ref=2e-5) == 0.0
        @test P._amp2db(2e-4; ref=2e-5) ≈ 20.0
        # upstream quirk: exact zero amplitude is replaced by 2e-12
        @test P._amp2db(0.0; ref=2e-5) == 20 * log10(2e-12 / 2e-5)
        @test P._db2amp(20.0; ref=2e-5) ≈ 2e-4
        @test P._db2amp(0.0; ref=1.0) == 1.0
    end

    @testset "_lininterp clamps like numpy.interp" begin
        xp = [1.0, 2.0, 4.0]; yp = [10.0, 20.0, 40.0]
        @test P._lininterp(3.0, xp, yp) == 30.0
        @test P._lininterp(0.0, xp, yp) == 10.0   # below → first
        @test P._lininterp(9.0, xp, yp) == 40.0   # above → last
        @test P._lininterp(2.0, xp, yp) == 20.0   # exact knot
    end
end
