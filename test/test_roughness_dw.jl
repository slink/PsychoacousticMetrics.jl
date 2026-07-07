include(joinpath(@__DIR__, "support", "am_generator.jl"))

@testset "roughness_dw" begin
    P = PsychoacousticMetrics

    @testset "anchor: 1 kHz, 60 dB, fmod 70 Hz, m=1 → 1 asper ±17%" begin
        # Gate adopted from MoSQITo test_roughness_dw.py:122-123. Their CI
        # analyzes a mislabeled-fs variant of this signal (see spec §Testing);
        # the clean anchor is expected to pass — a failure here is a finding
        # to investigate, never a tolerance to loosen.
        t = range(0, 1; length = 48000)          # python linspace(0, 1, 48000)
        xmod = sin.(2π * 70 .* t)
        stim = am_sine(xmod, 48000, 1000, 60)
        r = roughness_dw(stim, 48000)
        @test all(abs.(r.roughness_over_time .- 1.0) .<= 0.17)
        @test abs(r.roughness - 1.0) <= 0.17
    end

    @testset "result shape" begin
        t = range(0, 1; length = 48000)
        stim = am_sine(sin.(2π * 70 .* t), 48000, 1000, 60)
        r = roughness_dw(stim, 48000)
        nseg = length(r.roughness_over_time)
        @test size(r.specific_roughness) == (47, nseg)
        @test r.bark_axis == collect(0.5:0.5:23.5)
        @test length(r.time_axis) == nseg
        # 200 ms frames, default overlap 0.5 → hop 4800 samples: frame
        # centers step by 0.1 s (time_segmentation returns centers)
        nseg > 1 && @test r.time_axis[2] - r.time_axis[1] ≈ 4800 / 48000
    end

    @testset "silence → 0.0, no NaN" begin
        r = roughness_dw(zeros(48000), 48000)
        @test r.roughness == 0.0
        @test all(r.roughness_over_time .== 0.0)
        @test !any(isnan, r.specific_roughness)
    end

    @testset "calibration equivalence" begin
        t = range(0, 1; length = 48000)
        stim = am_sine(sin.(2π * 70 .* t), 48000, 1000, 60)
        r1 = roughness_dw(stim, 48000)
        r2 = roughness_dw(stim ./ 2, 48000; pa_per_unit = 2.0)
        @test r1.roughness ≈ r2.roughness
    end

    @testset "input validation" begin
        @test_throws ArgumentError roughness_dw(zeros(48000), 22050)      # fs too low
        @test_throws ArgumentError roughness_dw(zeros(100), 48000)        # < one frame
        @test_throws ArgumentError roughness_dw(zeros(48000), 48000; overlap = 1.0)
        @test_throws ArgumentError roughness_dw(zeros(48000), 48000; overlap = -0.1)
        @test_throws ArgumentError roughness_dw(zeros(48000), 48000; pa_per_unit = 0.0)

        # loud >15.5 kHz content exceeds the D&W 47-channel range;
        # upstream MoSQITo crashes on the same input (IndexError)
        t_hf = range(0, 1; length = 48000)
        hf = sin.(2π * 15600 .* t_hf)
        hf .*= 2e-5 * 10^(90 / 20) / std(hf; corrected = false)
        @test_throws DomainError roughness_dw(hf, 48000)
    end
end
