using PsychoacousticMetrics
using ZwickerLoudness
using Test

@testset "PsychoacousticMetrics" begin
    include("test_bark.jl")
    include("test_sharpness.jl")
    include("test_conformance_din45692.jl")
    include("test_crosscheck_mosqito.jl")
end
