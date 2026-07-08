# Stimulus synthesis for fluctuation-strength tests and generators.
#
# Conventions pinned against the SQAT reference (github.com/ggrecow/SQAT
# @ 00b449e40599f1c1ef4abe0596094552213d57eb), run under Octave 11.3.0 with
# the octave-signal 1.4.7 package, 2026-07-07. `FluctuationStrength_Osses2016`
# itself does not synthesize its validation stimuli (the validation scripts
# under validation/FluctuationStrength_Osses2016/*/run_validation_FS_*.m load
# pre-rendered .wav files from a Zenodo dataset, doi:10.5281/zenodo.7933206,
# which was not downloaded here — 390 MB, out of scope). The synthesis
# convention below was instead reverse-engineered from the *shipped* file
# `sound_files/reference_signals/RefSignal_FluctuationStrength_Osses2016.wav`
# (60 dB SPL, 1 kHz tone, 100% modulated at 4 Hz — the model's own worked
# example, documented to yield ~1 vacil), by direct sample inspection:
#
#   - rms(wav) == 0.02 Pa to 1e-10 (== 2e-5*10^(60/20)): the target SPL
#     calibrates the RMS of the FULL post-modulation signal, not the
#     pre-modulation carrier.
#   - wav[1:10] == reverse(wav[end-9:end]) (even symmetry about t=0) and
#     wav[1] is within 1e-7 of the signal's global max: both the envelope
#     and the carrier are COSINES (not sines) — a "raised-cosine" AM,
#     (1 + m*cos(2*pi*fmod*t)) * cos(2*pi*fc*t), so both factors peak
#     together at t=0. A sine convention would put x[1] near zero, which
#     the data rules out.
#   - No onset/offset ramp is present in the reference file: the global
#     max occurs at sample 1, and local RMS is periodic with period
#     1/fmod across the full 4s duration.
#
# This gives the amplitude-modulated tone's calibration constant: for
# x(t) = (1 + m*cos(2*pi*fmod*t)) * cos(2*pi*fc*t), E[x^2] = 1/2 + m^2/4
# (averaging sin^2/cos^2 to 1/2 and using E[cos^2(fmod)] = 1/2), so scaling
# by k = p / sqrt(1/2 + m^2/4) (p = 2e-5 * 10^(SPL/20)) gives rms(k*x) == p.
# Verified against the running oracle: FluctuationStrength_Osses2016 on this
# convention's anchor (1 kHz, 4 Hz, m=1, 60 dB, 4 s @ 44100 Hz) yields
# FSmean = 1.005363 vacil (expected ~1, tolerance 0.9-1.1 per the model's
# own documentation).
#
# FM and BBN conventions are not independently verifiable this way (no FM or
# noise reference file ships with SQAT), so they are inferred by analogy:
# cosine carrier/phase-integral for FM (continuing the "everything is a
# cosine, t=0 is an extremum" convention), and cosine modulator with
# post-modulation RMS calibration for BBN (mirroring the AM-tone finding
# that the target SPL binds the final, not the pre-modulation, signal).
#
# BBN uses a fixed plain-xorshift64 stream so the raw noise samples are
# reproduced exactly forever, independent of Julia's stdlib RNG. (The
# subsequent FFT band-limiting round-trip in fs_am_noise is deterministic
# for a given FFTW build, but may differ at the ULP level across FFTW
# versions/architectures; compare fixtures with a small relative
# tolerance, not bit equality.)
using FFTW: fft, ifft
using Statistics: mean

function _fs_xorshift64(seed::UInt64, n::Int)
    out = Vector{Float64}(undef, n)
    s = seed
    for i in 1:n
        s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17
        out[i] = 2.0 * (s >> 11) / (1 << 53) - 1.0   # uniform in [-1, 1)
    end
    return out
end

function fs_am_tone(fc, fmod, spl, dur, fs; mdepth = 1.0)
    t = (0:round(Int, dur * fs)-1) ./ fs
    x = (1 .+ mdepth .* cos.(2π * fmod .* t)) .* cos.(2π * fc .* t)
    p = 2e-5 * 10.0^(spl / 20)
    k = p / sqrt(0.5 + mdepth^2 / 4)   # calibrates rms(x) to p (post-modulation)
    return x .* k
end

function fs_fm_tone(fc, fmod, fdev, spl, dur, fs)
    t = (0:round(Int, dur * fs)-1) ./ fs
    x = cos.(2π * fc .* t .+ (fdev / fmod) .* sin.(2π * fmod .* t))
    p = 2e-5 * 10.0^(spl / 20)
    return x .* (p * sqrt(2))   # constant-envelope tone: rms(cos) = 1/sqrt(2)
end

function fs_am_noise(bw_hz, fmod, spl, dur, fs; seed = UInt64(0x9e3779b97f4a7c15))
    n = round(Int, dur * fs)
    noise = _fs_xorshift64(seed, n)
    # band-limit via FFT brick wall to [20, bw_hz] (both sides deterministic)
    X = fft(noise)
    freqs = (0:n-1) .* (fs / n)
    keep = (freqs .>= 20.0) .& (freqs .<= bw_hz)
    keep .|= (freqs .>= fs - bw_hz) .& (freqs .<= fs - 20.0)
    X[.!keep] .= 0
    noise = real.(ifft(X))
    t = (0:n-1) ./ fs
    x = (1 .+ cos.(2π * fmod .* t)) .* noise
    p = 2e-5 * 10.0^(spl / 20)
    return x .* (p / sqrt(mean(abs2, x)))         # level-normalized post-modulation
end

# Dispatcher shared by the fixture generator and by tests: maps a fixture
# case name to its exact synthesis call, so both sides see identical
# samples. Returns (signal, fs, method).
function synthesize_case(name::AbstractString)
    if name == "anchor_44k"
        return fs_am_tone(1000, 4, 60, 4.0, 44100), 44100.0, 1
    elseif name == "anchor_48k"
        return fs_am_tone(1000, 4, 60, 4.0, 48000), 48000.0, 1
    elseif name == "am_8hz_70db"
        return fs_am_tone(1000, 8, 70, 4.0, 44100), 44100.0, 1
    elseif name == "fm_tone_4hz"
        return fs_fm_tone(1500, 4, 700, 70, 4.0, 44100), 44100.0, 1
    elseif name == "am_bbn_4hz"
        return fs_am_noise(16000, 4, 60, 4.0, 44100), 44100.0, 1
    elseif name == "stationary_anchor"
        return fs_am_tone(1000, 4, 60, 4.0, 44100), 44100.0, 0
    elseif name == "short_fallback"
        sig, fs, _ = synthesize_case("anchor_44k")
        n = round(Int, 1.5 * fs)
        return sig[1:n], fs, 1
    elseif name == "tone_25bark"
        t = (0:round(Int, 4.0 * 44100)-1) ./ 44100
        p = 2e-5 * 10.0^(70 / 20)
        return (p * sqrt(2)) .* cos.(2π * 15800 .* t), 44100.0, 1
    else
        error("unknown fluctuation-strength fixture case: $name")
    end
end
