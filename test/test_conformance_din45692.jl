# DIN 45692:2009 chapter-6 conformance, both reference tables.
# Tolerance per the standard as implemented by MoSQITo
# validation_sharpness_din.py:257-263: within ±max(5 %, 0.05 acum).

include(joinpath(@__DIR__, "data", "din45692_band_levels.jl"))

din45692_tolerance(S_ref) = max(0.05 * S_ref, 0.05)

@testset "DIN 45692 conformance — narrowband (21 signals)" begin
    for (fc, S_ref, spl) in DIN45692_NARROWBAND
        @testset "$(fc) Hz" begin
            S = sharpness(zwicker_loudness(spl))
            @test abs(S - S_ref) <= din45692_tolerance(S_ref)
        end
    end
end

@testset "DIN 45692 conformance — broadband (20 signals)" begin
    for (fc, S_ref, spl) in DIN45692_BROADBAND
        @testset "$(fc) Hz" begin
            S = sharpness(zwicker_loudness(spl))
            @test abs(S - S_ref) <= din45692_tolerance(S_ref)
        end
    end
end
