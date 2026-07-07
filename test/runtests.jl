using PsychoacousticMetrics
using ZwickerLoudness
using Test

@testset "PsychoacousticMetrics" begin
    include("test_sharpness.jl")
    include("test_conformance_din45692.jl")
end
