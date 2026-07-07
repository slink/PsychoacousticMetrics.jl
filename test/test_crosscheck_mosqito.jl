# Weighting-transcription cross-check against MoSQITo. Expected values were
# computed by MoSQITo's own sharpness_din_from_loudness from THIS package's
# kernel output, with N = 0.1*sum(N') passed explicitly — so any deviation
# is a transcription bug, not a filterbank difference. Tight tolerance is
# intentional; do not loosen it.

include(joinpath(@__DIR__, "data", "crosscheck_cases.jl"))
include(joinpath(@__DIR__, "data", "mosqito_crosscheck.jl"))

@testset "MoSQITo weighting cross-check" begin
    results = Dict(name => zwicker_loudness(spl) for (name, spl) in CROSSCHECK_CASES)
    for (name, weighting, S_expected) in MOSQITO_CROSSCHECK
        S = sharpness(results[name]; weighting)
        @test isapprox(S, S_expected; rtol=1e-9)
    end
end
