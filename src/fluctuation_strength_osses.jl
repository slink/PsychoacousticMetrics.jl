# Fluctuation strength per Osses, García & Kohlrausch (2016),
# doi:10.1121/2.0000410, as maintained in SQAT's
# FluctuationStrength_Osses2016 (@ commit 00b449e).
# Implemented from the published description (Osses 2018 PhD thesis,
# App. B) with behaviors and parameter values pinned against the SQAT
# reference — see fluctuation_strength_tables.jl for the licensing note.

# MATLAB buffer(x, N, round(0.9N), 'nodelay') semantics (pinned against
# Octave, see .superpowers/sdd/buffer_pins.json): frames start at sample 1
# and advance by hop = N - round(0.9N); trailing partial frames are
# ZERO-PADDED, not dropped. The frame-count closed form is
# floor((L-V)/hop) + 1 — NOT ceil((L-V)/hop), which disagrees whenever
# L-V is an exact multiple of hop (verified against all 4 pinned (L,N,V)
# triples, 3 of 4 hit this case; see fs-oracle-pins.md Step 3.1). This
# also means a signal with L == N (the "stationary" case, where hop ==
# L-V exactly) produces 2 frames, not 1 — a documented pinned quirk, not
# a bug. Returns the frame matrix and 1-based start indices.
#
# Inputs shorter than V are OUTSIDE the reference's reachable domain (the
# reference forces N = length(insig) for signals under 2 s, so its buffer
# call always sees L >= N > V). There the oracle has no useful behavior
# to pin — verified against octave-signal 1.4.7's buffer on the pinned
# checkout: V-hop <= L < V yields an EMPTY N×0 matrix, and L < V-hop
# throws an internal "subscripts must be ... integers" error. We clamp to
# one zero-padded frame instead (a deliberate graceful extension, not an
# oracle behavior; it matches MATLAB's documented buffer padding).
function _fs_frames(sig::AbstractVector{Float64}, N::Int)
    V = round(Int, 0.9 * N)
    hop = N - V
    L = length(sig)
    # pinned formula (buffer_pins.json); max(1, ·) only fires for L < V,
    # i.e. outside the pinned domain — see the comment block above.
    nf = max(1, fld(L - V, hop) + 1)
    frames = zeros(N, nf)
    starts = Vector{Int}(undef, nf)
    for k in 1:nf
        i0 = (k - 1) * hop
        starts[k] = i0 + 1
        n = clamp(L - i0, 0, N)
        n > 0 && (frames[1:n, k] = @view sig[i0+1:i0+n])
    end
    return frames, starts
end

# Flat window with raised-cosine ramps, quirks of cos_ramp.m @ SQAT
# 00b449e preserved: 1-based attack (w[1] > 0, w[a] = 1); release spans
# a+1 samples and ends at EXACTLY 0; attack >= N falls back to
# round((N-1)/2).
function _cos_ramp_window(N::Int, fs::Real, ramp_ms::Real)
    a = round(Int, fs * ramp_ms / 1000)
    a >= N && (a = round(Int, (N - 1) / 2))
    w = ones(N)
    for i in 1:a
        w[i] *= 0.5 * (1 - cos(π * i / a))
    end
    for i in (N - a):N
        w[i] *= 0.5 * (1 - cos(π * (i - N) / a))
    end
    return w
end

# FFT overlap-free convolution equivalent to MATLAB filter(B,1,·) plus the
# reference's K/2 zero-pad-and-trim linear-phase delay compensation
# (il_PeripheralHearingSystem_t @ SQAT 00b449e).
function _apply_a0_fir(x::AbstractVector{Float64}, B::Vector{Float64})
    K = length(B) - 1
    L = length(x) + K ÷ 2
    n = nextpow(2, L + K)
    X = fft(vcat(x, zeros(n - length(x))))
    Bf = fft(vcat(B, zeros(n - length(B))))
    y = real.(ifft(X .* Bf))[1:L]
    return y[K÷2+1:end]
end

# Terhardt excitation patterns, FS variant (thesis App. B.2.2; behavioral
# pins from TerhardtExcitationPatterns.m @ SQAT 00b449e): dBFS = 94 with a
# +3 dB spectrum adjustment; audible bins 20 Hz–20 kHz; slopes S1 =
# -27 dB/Bark and S2 = -24 - 230/f + 0.2·L (only where negative);
# channel amplitudes from the floor/ceil half-Bark channels of each
# component; contributions gated by the MinBf threshold vector (a 49-entry
# band-edge/center list indexed by channel number — upstream quirk,
# preserved); e_i = 2N·real(ifft(·)); finally all 47 channels are rescaled
# so the level of their sum matches the level of the input frame.
# A frame with NO audible components (all bins at/below the hearing
# threshold, e.g. silence) short-circuits to an all-zero 47×N pattern
# BEFORE the level recalibration — preserved quirk: the recalibration
# would otherwise divide by an all-zero sum (0/0), and the pinned
# behavior for inaudible input is exact zeros, not NaN (see
# fs-oracle-pins.md Step 3.3: exact-zero output, silently, no error).
function _terhardt_excitation_fs(frame::AbstractVector{Float64}, fs::Real)
    N = length(frame)
    dBFS = 94.0
    df = fs / N
    N0 = round(Int, 20 / df) + 1
    Ntop = round(Int, 20e3 / df) + 1
    qb = N0:Ntop
    freqs = (qb .- 1) .* df
    barks = freq2bark.(freqs)                # Task-2 equivalence to Get_Bark

    spec = (10.0^((dBFS + 3) / 20) / N) .* fft(frame)
    Lg = abs.(@view spec[qb])
    LdB = 20 .* log10.(max.(Lg, floatmin()))
    minexc = _fs_threshold_db.(barks)
    audible = findall(LdB .> minexc)
    isempty(audible) && return zeros(47, N)

    # MinBf: threshold sampled at the sorted Zwicker band edge/center bins,
    # then indexed by CHANNEL number (quirk; il_calculate_MinBf @ 00b449e).
    minbf = _fs_minbf(N0, df, minexc)

    nL = length(audible)
    s2 = zeros(nL)
    for (j, k) in enumerate(audible)
        steep = -24.0 - 230.0 / freqs[k] + 0.2 * LdB[k]
        steep < 0 && (s2[j] = steep)
    end
    zlow = floor.(Int, 2 .* barks[audible])
    zhigh = ceil.(Int, 2 .* barks[audible])

    slopes = zeros(nL, 47)
    for j in 1:nL
        lev = LdB[audible[j]]
        b = barks[audible[j]]
        for ch in 1:min(zlow[j], 47)
            s = -27.0 * (b - ch / 2) + lev
            ch <= length(minbf) && s > minbf[ch] && (slopes[j, ch] = 10.0^(s / 20))
        end
        for ch in max(zhigh[j], 1):47
            s = s2[j] * (ch / 2 - b) + lev
            ch <= length(minbf) && s > minbf[ch] && (slopes[j, ch] = 10.0^(s / 20))
        end
    end

    ei = zeros(47, N)
    etmp = zeros(ComplexF64, N)
    for ch in 1:47
        fill!(etmp, 0.0 + 0.0im)
        for (j, k) in enumerate(audible)
            amp = if zlow[j] == ch || zhigh[j] == ch
                1.0
            elseif zhigh[j] > ch
                (ch < 47 ? slopes[j, ch+1] : 0.0) / Lg[k]
            else
                (ch > 1 ? slopes[j, ch-1] : 0.0) / Lg[k]
            end
            etmp[k+N0-1] = amp * spec[k+N0-1]
        end
        ei[ch, :] = 2 * N .* real.(ifft(etmp))
    end

    # per-frame level recalibration (no roughness counterpart)
    outsum = vec(sum(ei; dims = 1))
    ref_db = 20 * log10(sqrt(mean(abs2, frame))) # + dBFS on both sides cancels
    out_db = 20 * log10(sqrt(mean(abs2, outsum)))
    return ei .* 10.0^((ref_db - out_db) / 20)
end

# Sorted Zwicker band edge/center bin indices for the MinBf gate
# (il_calculate_MinBf @ SQAT 00b449e). `minexc` is indexed on the qb grid.
function _fs_minbf(N0::Int, df::Float64, minexc::Vector{Float64})
    centers = Float64[100, 200, 300, 400, 510, 630, 770, 920, 1080, 1270,
        1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700,
        9500, 12000, 15500]                     # Bark_raw rows 2:25, col 2
    edges = Float64[50, 150, 250, 350, 450, 570, 700, 840, 1000, 1170,
        1370, 1600, 1850, 2150, 2500, 2900, 3400, 4000, 4800, 5800, 7000,
        8500, 10500, 13500, 20000]              # Bark_raw rows 1:25, col 3
    zb = sort(vcat(round.(Int, centers ./ df) .- (N0 - 1) .+ 1,
        round.(Int, edges ./ df) .- (N0 - 1) .+ 1))
    return minexc[zb]
end
