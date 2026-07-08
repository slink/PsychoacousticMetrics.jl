# Reproducer for the "FFT dust" adjudication in
# test/test_crosscheck_mosqito_roughness.jl (see that header for the full
# mechanism and the matching MoSQITo/numpy-side numbers). For each affected
# cross-check case this script re-runs the frame pipeline with per-channel
# instrumentation and, in every frame, flags channels whose bandpass
# envelope hBP is analytically zero but numerically pure FFT round-off dust:
# maxabs(hBP) many orders below that channel's mean envelope h0, with zero
# exact zeros — so the ki exact-zero guard never fires and ki becomes a
# correlation of a real envelope against rounding noise.
#
# Stages 1-3 below are an instrumented transcription of
# _roughness_dw_frame (src/roughness_dw.jl); segmentation, spectrum, and
# tables reuse the package internals directly. To catch drift, each frame's
# recomputed R is asserted equal to roughness_dw's own output — if the
# source function changes, this script fails loudly instead of probing a
# stale algorithm.
#
# Usage: julia --project=. scripts/probe_fft_dust.jl
#
# Reference numbers (2026-07-07, this machine, FFTW; see test header):
#   anchor_1k_70  frame 8, channel 27: h0 = 1.616e-5, maxabs(hBP) ~ 4.2e-22,
#                 ki[25] = 7.911e-3, R_spec[25] = 9.120e-6
#   fc250_fm40    frame 4, channel 14: maxabs(hBP) ~ 9.8e-22, ki[12] = -0.0195
#   fc2000_fm100  frame 2, channel 37: maxabs(hBP) ~ 2.5e-22, ki[35] = -0.2415

using PsychoacousticMetrics
using FFTW: fft, ifft
using Statistics: mean, cor
using Printf: @printf

const P = PsychoacousticMetrics

include(joinpath(@__DIR__, "..", "test", "support", "am_generator.jl"))

# A channel is "dust" when its envelope is nonzero (guard can't fire) yet
# sits >= 12 orders of magnitude below its own mean envelope: pure
# round-off, no analytic content.
const DUST_RATIO = 1e-12

# Instrumented replica of _roughness_dw_frame: returns per-channel
# h0, maxabs(hBP), exact-zero count of hBP, mdepth, ki, R_spec, and R.
function probe_frame(spec::Vector{ComplexF64}, freq_axis, fs, gzi47, H)
    L = length(spec)
    spec2 = vcat(spec, reverse(spec))
    n = 2L
    bark_axis = P.freq2bark.(freq_axis)

    a0 = zeros(n)
    a0[1:L] = P._db2amp.(P._a0_ear_db.(bark_axis); ref = 1.0)
    spec2 = a0 .* spec2

    mod_ = abs.(spec2[1:L])
    spec_dB = P._amp2db.(mod_; ref = 2e-5)
    threshold = P._ltq_roughness.(bark_axis)
    audible = findall(spec_dB .> threshold)
    n_aud = length(audible)

    s1 = -27.0
    s2 = [min(-24.0 - 230.0 / freq_axis[audible[k]] + 0.2 * spec_dB[audible[k]], 0.0)
          for k in 1:n_aud]
    zi = (1:47) ./ 2
    zb = P.bark2freq.(zi) .* n ./ fs
    binaxis = collect(1.0:L)
    minexcit = [P._lininterp(z, binaxis, threshold) for z in zb]

    ch_low = [floor(Int, 2 * bark_axis[audible[i]]) - 1 for i in 1:n_aud]
    ch_high = [ceil(Int, 2 * bark_axis[audible[i]]) - 1 for i in 1:n_aud]

    slopes = zeros(n_aud, 47)
    for k in 1:n_aud
        lev = spec_dB[audible[k]]
        b = bark_axis[audible[k]]
        for j in 0:ch_low[k]
            sl = s1 * (b - (j + 1) * 0.5) + lev
            sl > minexcit[j + 1] && (slopes[k, j + 1] = P._db2amp(sl; ref = 2e-5))
        end
        for j in ch_high[k]:46
            sl = s2[k] * ((j + 1) * 0.5 - b) + lev
            sl > minexcit[j + 1] && (slopes[k, j + 1] = P._db2amp(sl; ref = 2e-5))
        end
    end

    hBP = zeros(47, n)
    h0s = zeros(47)
    mdepth = zeros(47)
    for i in 0:46
        exc = zeros(ComplexF64, n)
        for j in 1:n_aud
            ind = audible[j]
            ampl = if ch_low[j] == i || ch_high[j] == i
                1.0
            elseif ch_high[j] > i
                slopes[j, i + 2] / mod_[ind]
            else
                slopes[j, i] / mod_[ind]
            end
            exc[ind] = ampl * spec2[ind]
        end
        temporal = abs.(n .* real.(ifft(exc)))
        h0 = mean(temporal)
        h0s[i + 1] = h0
        env_spec = fft(temporal .- h0) .* @view(H[i + 1, :])
        hBP[i + 1, :] = 2 .* real.(ifft(env_spec))
        if h0 > 0
            mdepth[i + 1] = min(sqrt(mean(abs2, @view(hBP[i + 1, :]))) / h0, 1.0)
        end
    end

    ki = zeros(47)
    for i in 1:45
        if all(!=(0.0), @view(hBP[i, :])) && all(!=(0.0), @view(hBP[i + 2, :]))
            ki[i] = cor(@view(hBP[i, :]), @view(hBP[i + 2, :]))
        end
    end

    R_spec = Vector{Float64}(undef, 47)
    R_spec[1] = gzi47[1] * (mdepth[1] * ki[1])^2
    R_spec[2] = gzi47[2] * (mdepth[2] * ki[2])^2
    for i in 3:45
        R_spec[i] = gzi47[i] * (mdepth[i] * ki[i] * ki[i - 2])^2
    end
    R_spec[46] = gzi47[46] * (mdepth[46] * ki[44])^2
    R_spec[47] = gzi47[47] * (mdepth[47] * ki[45])^2

    maxabs_hBP = [maximum(abs, @view(hBP[i, :])) for i in 1:47]
    nzeros_hBP = [count(==(0.0), @view(hBP[i, :])) for i in 1:47]
    return 0.25 * sum(R_spec), h0s, maxabs_hBP, nzeros_hBP, mdepth, ki, R_spec
end

# (name, fc, fmod, duration, overlap) — the three dust-adjudicated cases,
# synthesized exactly as in scripts/generate_mosqito_roughness_crosscheck.jl
cases = [
    ("anchor_1k_70", 1000.0, 70.0, 1.0, 0.5),
    ("fc250_fm40", 250.0, 40.0, 1.5, 0.0),
    ("fc2000_fm100", 2000.0, 100.0, 1.5, 0.0),
]

fs = 48000
for (name, fc, fm, dur, overlap) in cases
    t = range(0, dur; length = Int(round(dur * fs)))
    sig = am_sine(sin.(2π * fm .* t), fs, fc, 60)

    ref = roughness_dw(sig, fs; overlap = overlap)

    nperseg = floor(Int, 0.2 * fs)
    noverlap = floor(Int, overlap * nperseg)
    frames, _ = P._time_segmentation(sig, fs, nperseg, noverlap)
    L = nperseg ÷ 2
    freq_axis = collect((1:L) .* (fs / nperseg))
    H = P._h_weighting(2L, fs)
    gzi47 = P._gzi.((1:47) ./ 2)

    println("=== $name (fc=$(fc) Hz, fmod=$(fm) Hz, overlap=$overlap, ",
            "$(size(frames, 2)) frames) ===")
    for s in 1:size(frames, 2)
        spec = P._blackman_spectrum(@view frames[:, s])
        R, h0s, maxabs_hBP, nzeros_hBP, mdepth, ki, R_spec =
            probe_frame(spec, freq_axis, fs, gzi47, H)
        # drift guard: instrumented replica must match the package pipeline
        R == ref.roughness_over_time[s] || error(
            "probe drift: frame $s R = $R vs roughness_dw " *
            "$(ref.roughness_over_time[s]) — re-sync probe_frame with " *
            "_roughness_dw_frame")
        for c in 1:47
            if maxabs_hBP[c] > 0 && maxabs_hBP[c] < DUST_RATIO * h0s[c]
                @printf("frame %d channel %2d: DUST  h0 = %.6e  maxabs(hBP) = %.3e (%.0f orders below h0)  exact zeros = %d/%d\n",
                        s, c, h0s[c], maxabs_hBP[c],
                        log10(h0s[c] / maxabs_hBP[c]), nzeros_hBP[c], 2L)
                for k in (c - 2, c)          # ki entries correlating channel c
                    if 1 <= k <= 45
                        @printf("    ki[%2d] (ch %2d x ch %2d) = %+.4e   R_spec[%2d] = %.4e\n",
                                k, k, k + 2, ki[k], k, R_spec[k])
                    end
                end
            end
        end
    end
    println()
end
println("Interpretation: DUST channels have no exact zeros, so the ki exact-",
        "zero guard passes them; their ki (and thus R_spec) is a correlation ",
        "against FFT rounding noise and differs across FFT implementations. ",
        "See test/test_crosscheck_mosqito_roughness.jl for the adjudication.")
