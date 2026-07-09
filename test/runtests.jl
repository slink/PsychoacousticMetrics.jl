using PsychoacousticMetrics
using ZwickerLoudness
using Test

@testset "PsychoacousticMetrics" begin
    include("test_bark.jl")
    include("test_roughness_tables.jl")
    include("test_fluctuation_strength_tables.jl")
    include("test_fluctuation_strength_stages.jl")
    include("test_fluctuation_strength.jl")
    include("test_conformance_fs_thesis.jl")
    include("test_crosscheck_sqat_fs.jl")
    include("test_sharpness.jl")
    include("test_conformance_din45692.jl")
    include("test_crosscheck_mosqito.jl")
    include("test_roughness_dw.jl")
    include("test_conformance_dw_fig3.jl")
    include("test_crosscheck_mosqito_roughness.jl")
end
