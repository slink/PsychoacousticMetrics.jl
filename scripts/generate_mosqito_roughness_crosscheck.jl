# One-time generator for test/data/mosqito_roughness_crosscheck.jl.
# Synthesizes the cross-check signals with the package's own test generator,
# hands them to MoSQITo (via uv), vendors the per-frame results.
# Requires uv and network on first run.
# Usage: julia scripts/generate_mosqito_roughness_crosscheck.jl
using Statistics: std

include(joinpath(@__DIR__, "..", "test", "support", "am_generator.jl"))

# (name, fs_synth, fc, fmod, duration, spl, fs_analyze, overlap)
# "mislabeled_fs" replicates MoSQITo's literal pytest call (synthesized at
# 48 kHz, analyzed as 44.1 kHz, overlap 0) — fidelity check of the fs path.
# "fc2000_fm100" probes a mid-cluster point of the fig-3 fc=2000 Hz
# discrepancy vs Zwicker & Fastl: if MoSQITo agrees with us here, the
# fig-3 deviation is upstream's model behavior, not our transcription.
cases = [
    ("anchor_1k_70",   48000, 1000.0, 70.0,  1.0, 60, 48000.0, 0.5),
    ("mislabeled_fs",  48000, 1000.0, 70.0,  1.0, 60, 44100.0, 0.0),
    ("fc250_fm40",     48000,  250.0, 40.0,  1.5, 60, 48000.0, 0.0),
    ("fc2000_fm100",   48000, 2000.0, 100.0, 1.5, 60, 48000.0, 0.0),
]

jsonvec(v) = string("[", join(string.(v), ","), "]")
tmp = tempname() * ".json"
open(tmp, "w") do io
    entries = String[]
    for (name, fss, fc, fm, dur, spl, fsa, ov) in cases
        t = range(0, dur; length = Int(round(dur * fss)))
        sig = am_sine(sin.(2π * fm .* t), fss, fc, spl)
        push!(entries, """{"name": "$name", "fs": $fsa, "overlap": $ov, "signal": $(jsonvec(sig))}""")
    end
    print(io, "[", join(entries, ","), "]")
end
out = joinpath(@__DIR__, "..", "test", "data", "mosqito_roughness_crosscheck.jl")
run(`uv run --with mosqito --with numpy --with matplotlib python $(joinpath(@__DIR__, "crosscheck_roughness.py")) $tmp $out`)
