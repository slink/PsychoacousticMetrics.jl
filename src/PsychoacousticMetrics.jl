module PsychoacousticMetrics

using ZwickerLoudness: ZwickerResult
using Statistics: mean, std, cor
using FFTW: fft, ifft

export sharpness, roughness_dw, RoughnessResult

include("bark.jl")
include("roughness_dw_tables.jl")
include("fluctuation_strength_tables.jl")
include("fluctuation_strength_fir.jl")
include("fluctuation_strength_osses.jl")
include("weighting_fastl.jl")
include("sharpness.jl")
include("roughness_dw.jl")

end # module
