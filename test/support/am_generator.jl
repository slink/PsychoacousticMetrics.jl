# AM sine stimulus, transcribed from MoSQITo utils/am_sine_generator.py
# (commit d990c33f94f1). NOTE: level normalization uses the standard
# deviation (numpy std, ddof=0), not the RMS — transcribed as-is.
using Statistics: std

function am_sine(xmod::AbstractVector{<:Real}, fs::Real, fc::Real, spl::Real)
    Nt = length(xmod)
    T = Nt / fs
    dt = 1 / fs
    t = range(0, T - dt; length = floor(Int, T * fs))   # python linspace(0, T-dt, int(T*fs)); int() truncates
    xc = sin.(2π * fc .* t)
    y = (1 .+ xmod) .* xc
    A = 2e-5 * 10.0^(spl / 20)
    return y .* (A / std(y; corrected = false))
end
