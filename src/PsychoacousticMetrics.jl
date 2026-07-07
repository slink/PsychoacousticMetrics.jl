module PsychoacousticMetrics

using ZwickerLoudness: ZwickerResult

export sharpness

include("weighting_fastl.jl")
include("sharpness.jl")

end # module
