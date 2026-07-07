# Daniel & Weber roughness, transcribed 2026-07-07 from MoSQITo commit
# d990c33f94f1 (roughness_dw.py, _roughness_dw_main_calc.py,
# sound_level_meter/comp_spectrum.py, utils/time_segmentation.py).
# Reference: R. Daniel, D. Weber, "Psychoacoustical roughness:
# implementation of an optimized model", Acustica 83, 1997.

"""
    RoughnessResult

Result of [`roughness_dw`](@ref). Fields: `roughness` (overall [asper] —
arithmetic mean over frames; MoSQITo returns only per-frame values),
`roughness_over_time` (per 200 ms frame), `specific_roughness`
(47 half-Bark channels × frames [asper/Bark]), `bark_axis`
(0.5:0.5:23.5), `time_axis` (frame centers [s]).
"""
struct RoughnessResult
    roughness::Float64
    roughness_over_time::Vector{Float64}
    specific_roughness::Matrix{Float64}
    bark_axis::Vector{Float64}
    time_axis::Vector{Float64}
end

# MoSQITo time_segmentation.py quirk, preserved: the `noverlap` argument is
# the HOP SIZE (block l covers samples l*hop+1 : l*hop+nperseg), and
# noverlap == 0 is remapped to hop = nperseg (no overlap). Frame times are
# the mean of the block's sample times, i.e. frame centers.
function _time_segmentation(sig::Vector{Float64}, fs::Real, nperseg::Int, noverlap::Int)
    hop = noverlap == 0 ? nperseg : noverlap
    nseg = length(sig) < nperseg ? 0 : (length(sig) - nperseg) ÷ hop + 1
    frames = Matrix{Float64}(undef, nperseg, nseg)
    times = Vector{Float64}(undef, nseg)
    for l in 0:(nseg - 1)
        frames[:, l + 1] = @view sig[l * hop + 1 : l * hop + nperseg]
        times[l + 1] = (l * hop + (nperseg - 1) / 2) / fs
    end
    return frames, times
end

# One-sided complex spectrum per MoSQITo comp_spectrum.py (window="blackman",
# db=False): Blackman window normalized by its sum, fft, keep bins 1..n÷2
# (python 0-based 0:n//2), times the empirical 1.42 amplitude factor.
# numpy.blackman(M): 0.42 - 0.5cos(2πk/(M-1)) + 0.08cos(4πk/(M-1)), k=0..M-1.
function _blackman_spectrum(frame::AbstractVector{Float64})
    n = length(frame)
    w = [0.42 - 0.5 * cos(2π * k / (n - 1)) + 0.08 * cos(4π * k / (n - 1)) for k in 0:(n - 1)]
    w ./= sum(w)
    return fft(frame .* w)[1:n ÷ 2] .* 1.42
end

# Per-frame Daniel & Weber calculation (_roughness_dw_main_calc.py).
# `spec` is the one-sided complex spectrum (length L), `freq_axis` its L
# frequencies, `gzi47` the 47 channel weights, `H` the 47 × 2L weighting
# matrix. Returns (R, R_spec[47]).
function _roughness_dw_frame(spec::Vector{ComplexF64}, freq_axis::AbstractVector{Float64},
                             fs::Real, gzi47::Vector{Float64}, H::Matrix{Float64})
    L = length(spec)
    spec2 = vcat(spec, reverse(spec))            # python concatenate((spec, spec[len(spec)::-1]))
    n = 2L
    bark_axis = freq2bark.(freq_axis)

    # a0 ear weighting on the lower half; upper half of spec2 becomes 0 —
    # upstream behaviour (a0 is zero there), preserved.
    a0 = zeros(n)
    a0[1:L] = _db2amp.(_a0_ear_db.(bark_axis); ref = 1.0)
    spec2 = a0 .* spec2

    mod_ = abs.(spec2[1:L])
    spec_dB = _amp2db.(mod_; ref = 2e-5)
    threshold = _ltq_roughness.(bark_axis)
    audible = findall(spec_dB .> threshold)
    n_aud = length(audible)

    # Upstream (MoSQITo _roughness_dw_main_calc.py) crashes with an
    # IndexError on the same inputs: audible components at z ≥ 24 Bark
    # (≥ 15.5 kHz at high level) index past the 47-channel arrays. The
    # Daniel & Weber channel structure ends at 23.5 Bark, so we fail
    # loudly instead of inventing behavior upstream never validated.
    if !isempty(audible) && bark_axis[audible[end]] >= 24
        throw(DomainError(freq_axis[audible[end]],
            "audible component above 24 Bark (~15.5 kHz) — outside the Daniel & Weber channel range; low-pass the signal below 15.5 kHz first"))
    end

    # ---- stage 1: excitation patterns (Terhardt slopes) ----
    s1 = -27.0
    s2 = [min(-24.0 - 230.0 / freq_axis[audible[k]] + 0.2 * spec_dB[audible[k]], 0.0)
          for k in 1:n_aud]
    zi = (1:47) ./ 2
    zb = bark2freq.(zi) .* n ./ fs               # channel centers in bin units
    binaxis = collect(1.0:L)                     # python nZ = arange(1, n//2+1)
    minexcit = [_lininterp(z, binaxis, threshold) for z in zb]

    # python ch_low/ch_high are 0-based channel ids; keep the VALUES and
    # add 1 only when indexing Julia arrays.
    ch_low = [floor(Int, 2 * bark_axis[audible[i]]) - 1 for i in 1:n_aud]
    ch_high = [ceil(Int, 2 * bark_axis[audible[i]]) - 1 for i in 1:n_aud]

    slopes = zeros(n_aud, 47)
    for k in 1:n_aud
        lev = spec_dB[audible[k]]
        b = bark_axis[audible[k]]
        for j in 0:ch_low[k]                     # python range(0, ch_low+1)
            sl = s1 * (b - (j + 1) * 0.5) + lev
            sl > minexcit[j + 1] && (slopes[k, j + 1] = _db2amp(sl; ref = 2e-5))
        end
        for j in ch_high[k]:46                   # python range(ch_high, 47)
            sl = s2[k] * ((j + 1) * 0.5 - b) + lev
            sl > minexcit[j + 1] && (slopes[k, j + 1] = _db2amp(sl; ref = 2e-5))
        end
    end

    # ---- stage 2: envelopes and modulation depths ----
    hBP = zeros(47, n)
    mdepth = zeros(47)
    for i in 0:46                                # python channel loop, 0-based
        exc = zeros(ComplexF64, n)
        for j in 1:n_aud
            ind = audible[j]
            ampl = if ch_low[j] == i || ch_high[j] == i
                1.0
            elseif ch_high[j] > i
                slopes[j, i + 2] / mod_[ind]     # python slopes[j, i+1]
            else
                slopes[j, i] / mod_[ind]         # python slopes[j, i-1]
            end
            exc[ind] = ampl * spec2[ind]
        end
        temporal = abs.(n .* real.(ifft(exc)))
        h0 = mean(temporal)
        env_spec = fft(temporal .- h0) .* @view(H[i + 1, :])
        hBP[i + 1, :] = 2 .* real.(ifft(env_spec))
        if h0 > 0
            mdepth[i + 1] = min(sqrt(mean(abs2, @view(hBP[i + 1, :]))) / h0, 1.0)
        end
    end

    # ---- stage 3: cross-correlation and summation ----
    # The all-nonzero guard is the ONLY thing preventing R = NaN for silent
    # channels: cor() of constant vectors is 0/0 = NaN and mdepth = 0 does
    # not rescue 0*NaN. Transcribed from python line 173.
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

    return 0.25 * sum(R_spec), R_spec
end

"""
    roughness_dw(signal, fs; overlap=0.5, pa_per_unit=1.0) -> RoughnessResult

Psychoacoustic roughness [asper] of a time signal, per Daniel & Weber
(Acustica 83, 1997) **as implemented by MoSQITo**, which differs from the
canonical MATLAB implementation in the placement of the gzi weighting
(outside the square). Signal is multiplied by `pa_per_unit` to obtain
pascals. `overlap` is MoSQITo's frame-overlap coefficient (its `noverlap`
parameter is a hop size upstream; `overlap=0.5` gives 50 % overlapping
200 ms frames, `overlap=0` gives none). For values other than 0 and 0.5,
note the inversion: the effective overlap fraction between consecutive
frames is `1 − overlap`. Anchor: 1 kHz, 60 dB tone, 100 %
amplitude-modulated at 70 Hz → 1 asper.

!!! warning
    The Daniel & Weber channel structure covers 0.5–23.5 Bark. Signals
    with audible content above 24 Bark (roughly ≥ 15.5 kHz at high SPL)
    throw a `DomainError`; upstream MoSQITo crashes (`IndexError`) on the
    same inputs.
"""
function roughness_dw(signal::AbstractVector{<:Real}, fs::Real;
                      overlap::Real = 0.5, pa_per_unit::Real = 1.0)
    fs >= 44100 || throw(ArgumentError(
        "fs = $fs Hz is below the validated range; resample to ≥ 44.1 kHz first"))
    0 <= overlap < 1 || throw(ArgumentError("overlap must be in [0, 1), got $overlap"))
    pa_per_unit > 0 || throw(ArgumentError("pa_per_unit must be positive, got $pa_per_unit"))
    nperseg = floor(Int, 0.2 * fs)
    length(signal) >= nperseg || throw(ArgumentError(
        "signal too short: need at least one 200 ms frame ($nperseg samples), got $(length(signal))"))

    sig = Float64.(signal) .* pa_per_unit
    noverlap = floor(Int, overlap * nperseg)
    frames, time_axis = _time_segmentation(sig, fs, nperseg, noverlap)
    L = nperseg ÷ 2
    freq_axis = collect((1:L) .* (fs / nperseg))
    H = _h_weighting(2L, fs)
    gzi47 = _gzi.((1:47) ./ 2)

    nseg = size(frames, 2)
    R = Vector{Float64}(undef, nseg)
    R_spec = Matrix{Float64}(undef, 47, nseg)
    for s in 1:nseg
        spec = _blackman_spectrum(@view frames[:, s])
        R[s], R_spec[:, s] = _roughness_dw_frame(spec, freq_axis, fs, gzi47, H)
    end
    return RoughnessResult(mean(R), R, R_spec, collect(0.5:0.5:23.5), time_axis)
end
