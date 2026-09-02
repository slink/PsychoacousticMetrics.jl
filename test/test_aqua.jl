using Aqua

@testset "Aqua.jl quality checks" begin
    Aqua.test_all(PsychoacousticMetrics)
end
