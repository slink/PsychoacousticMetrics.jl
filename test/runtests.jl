using PsychoacousticMetrics
using ZwickerLoudness
using Test

@testset "PsychoacousticMetrics" begin
    @testset "package loads" begin
        @test PsychoacousticMetrics isa Module
    end
end
