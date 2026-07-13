# One-time generator for test/data/sqat_pa_crosscheck.jl.
#
# Two fixture families, both via the pinned SQAT reference
# (github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb,
# CC BY-NC) through Octave (scripts/run_sqat_pa.m):
#
#   1. SQAT_PA_FORMULA_CASES: a non-negative (N, S, R, FS) grid (840
#      cases; the reference goes complex on negative inputs — see
#      .superpowers/sdd/pa-oracle-pins.md — so the grid excludes them)
#      through `PsychoacousticAnnoyance_Widmann1992_from_percentile`.
#   2. SQAT_PA_SIGNAL_CASES: steady/AM tones synthesized with this
#      package's own stimulus generator (test/support/fs_stimuli.jl's
#      `fs_am_tone`, reused as-is: `mdepth=0` collapses it to a
#      correctly-calibrated pure tone since the modulation term
#      `1 + mdepth*cos(...)` becomes the constant 1 — verified by direct
#      RMS check, see pins file) through the signal-level
#      `PsychoacousticAnnoyance_Widmann1992`, vendoring PA together with
#      SQAT's own N5/S5/R5/FS5 percentile values so any wrapper deviation
#      (Task 3) is attributable per-component before any tolerance is set.
#
# Requires: Octave with the octave-signal package installed, and the
# pinned checkout at /tmp/sqat-pinned (see
# .superpowers/sdd/pa-oracle-pins.md Step 1).
#
# Usage: julia --project=. scripts/generate_sqat_pa_crosscheck.jl
include(joinpath(@__DIR__, "..", "test", "support", "fs_stimuli.jl"))

const SQAT_DIR = get(ENV, "SQAT_DIR", "/tmp/sqat-pinned")
isdir(SQAT_DIR) || error("pinned SQAT checkout not found at $SQAT_DIR — clone+checkout per Step 1 first")

driver = joinpath(@__DIR__, "run_sqat_pa.m")
workdir = mktempdir()

octave_version = strip(read(`octave --eval "disp(OCTAVE_VERSION)"`, String))
signal_version = strip(read(`octave --eval "pkg load signal; disp(ver('signal').Version)"`, String))

jsonvec(v) = "[" * join(v, ", ") * "]"

function parse_json_number_array(text, key)
    # Minimal hand-rolled extraction (matches generate_sqat_fs_crosscheck.jl):
    # the driver's JSON is a flat object of scalars/number arrays, so a
    # simple bracket scan suffices without a JSON dependency.
    m = match(Regex("\"$key\"\\s*:\\s*(\\[[^\\]]*\\]|[-0-9.eE]+)"), text)
    m === nothing && error("key $key not found")
    raw = m.captures[1]
    if startswith(raw, "[")
        return parse.(Float64, split(strip(raw, ['[', ']']), ","))
    else
        return parse(Float64, raw)
    end
end

# --- 1. Formula grid -------------------------------------------------------
#
# All non-negative; N=0 exercises the reference's Inf/NaN zeroing; S values
# straddle the ws threshold (1.75 strict `>`, see pins file) including one
# value 1e-7 above it; R/FS span zero (no roughness/fluctuation influence)
# up to values above the AM-tone anchor magnitude.
const N_VALS = [0.0, 1e-6, 0.5, 1.0, 4.0, 16.0, 64.0]
const S_VALS = [0.5, 1.0, 1.75, 1.7500001, 2.5, 4.0]
const R_VALS = [0.0, 0.05, 0.3, 1.0]
const FS_VALS = [0.0, 0.05, 0.3, 1.0, 2.0]

formula_grid = NamedTuple{(:N, :S, :R, :FS), NTuple{4, Float64}}[]
for N in N_VALS, S in S_VALS, R in R_VALS, FS in FS_VALS
    push!(formula_grid, (N = N, S = S, R = R, FS = FS))
end
@assert length(formula_grid) == 840

formula_manifest_path = joinpath(workdir, "formula_manifest.json")
write(formula_manifest_path, "[" * join(
    ["""{"N": $(c.N), "S": $(c.S), "R": $(c.R), "FS": $(c.FS)}""" for c in formula_grid], ",") * "]")

formula_output_path = joinpath(workdir, "formula_output.json")
run(`octave --no-gui $driver formula $formula_manifest_path $formula_output_path $SQAT_DIR`)

formula_text = read(formula_output_path, String)
formula_N = parse_json_number_array(formula_text, "N")
formula_S = parse_json_number_array(formula_text, "S")
formula_R = parse_json_number_array(formula_text, "R")
formula_FS = parse_json_number_array(formula_text, "FS")
formula_pa = parse_json_number_array(formula_text, "pa")
@assert length(formula_pa) == 840
# Echo-back sanity check: Octave received exactly the grid we sent (order
# and values), so vendoring the Julia-side literals below is safe.
for (i, c) in enumerate(formula_grid)
    @assert formula_N[i] == c.N && formula_S[i] == c.S && formula_R[i] == c.R && formula_FS[i] == c.FS
end

# --- 2. Signal cases --------------------------------------------------------
#
# name => (fc, fmod, spl_db, dur_s, fs, mdepth)
const SIGNAL_CASE_DEFS = [
    ("steady_1k_40db", 1000, 0, 40, 5.0, 44100),      # published anchor: 1 kHz, 40 dB -> PA = 1 au
    ("steady_1k_60db", 1000, 0, 60, 5.0, 44100),
    ("am_4hz_60db", 1000, 4, 60, 5.0, 44100),         # FS anchor signal (mdepth=1) — informational, modulated
    ("steady_1k_40db_48k", 1000, 0, 40, 5.0, 48000),  # anchor tone at the other supported fs
]
signal_mdepth(name) = name == "am_4hz_60db" ? 1.0 : 0.0

signal_manifest_entries = String[]
for (name, fc, fmod, spl, dur, fs) in SIGNAL_CASE_DEFS
    sig = fs_am_tone(fc, fmod, spl, dur, fs; mdepth = signal_mdepth(name))
    bin_path = joinpath(workdir, "$name.bin")
    open(bin_path, "w") do io
        write(io, sig)
    end
    push!(signal_manifest_entries, """{"name": "$name", "fs": $(Float64(fs)), "signal_path": "$bin_path"}""")
end
signal_manifest_path = joinpath(workdir, "signal_manifest.json")
write(signal_manifest_path, "[" * join(signal_manifest_entries, ",") * "]")

signal_output_dir = joinpath(workdir, "signal_out")
run(`octave --no-gui $driver signal $signal_manifest_path $signal_output_dir $SQAT_DIR`)

# --- 3. Write fixtures -------------------------------------------------------

lines = String[]
push!(lines, "# AUTOGENERATED by scripts/generate_sqat_pa_crosscheck.jl — do not edit.")
push!(lines, "#")
push!(lines, "# SQAT_PA_FORMULA_CASES: (N [sone], S [acum], R [asper], FS [vacil], pa [au])")
push!(lines, "# computed by SQAT's PsychoacousticAnnoyance_Widmann1992_from_percentile (CC")
push!(lines, "# BY-NC, github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb)")
push!(lines, "# under Octave $octave_version (octave-signal $signal_version). Pure scalar")
push!(lines, "# arithmetic — any disagreement with `psychoacoustic_annoyance_widmann` is a")
push!(lines, "# transcription bug, not a tolerance question. Grid is non-negative only (the")
push!(lines, "# reference goes complex on negative inputs, see pa-oracle-pins.md).")
push!(lines, "const SQAT_PA_FORMULA_CASES = [")
for i in 1:length(formula_grid)
    c = formula_grid[i]
    push!(lines, "    (N = $(c.N), S = $(c.S), R = $(c.R), FS = $(c.FS), pa = $(formula_pa[i])),")
end
push!(lines, "]")
push!(lines, "")
push!(lines, "# SQAT_PA_SIGNAL_CASES: (name, fs, pa, N5 [sone], S5 [acum], R5 [asper],")
push!(lines, "# FS5 [vacil]) computed by SQAT's signal-level PsychoacousticAnnoyance_Widmann1992")
push!(lines, "# (same pin), LoudnessField=0 (free field, matches this package's")
push!(lines, "# ZwickerLoudness field_type=:free default), time_skip=0 (steady-state signals,")
push!(lines, "# no onset transient to discard), on signals synthesized by this package's own")
push!(lines, "# fs_am_tone (test/support/fs_stimuli.jl) with mdepth=0 for the steady tones")
push!(lines, "# (collapses the AM envelope to the constant 1, verified by direct RMS check —")
push!(lines, "# see pa-oracle-pins.md) and mdepth=1 for am_4hz_60db (the FS anchor signal,")
push!(lines, "# vendored for information; not a steady tone). pa = SQAT's OUT.ScalarPA; N5/S5/")
push!(lines, "# R5/FS5 = OUT.L.N5 / OUT.S.S5 / OUT.R.R5 / OUT.FS.FS5 — vendored so any wrapper")
push!(lines, "# deviation (Task 3) is attributable per-component before any tolerance is set.")
push!(lines, "const SQAT_PA_SIGNAL_CASES = [")
for (name, fc, fmod, spl, dur, fs) in SIGNAL_CASE_DEFS
    text = read(joinpath(signal_output_dir, "$name.json"), String)
    pa = parse_json_number_array(text, "pa")
    N5 = parse_json_number_array(text, "N5")
    S5 = parse_json_number_array(text, "S5")
    R5 = parse_json_number_array(text, "R5")
    FS5 = parse_json_number_array(text, "FS5")
    push!(lines, "    (name = \"$name\", fs = $(Float64(fs)), pa = $pa, N5 = $N5, S5 = $S5, R5 = $R5, FS5 = $FS5),")
end
push!(lines, "]")

out_path = joinpath(@__DIR__, "..", "test", "data", "sqat_pa_crosscheck.jl")
write(out_path, join(lines, "\n") * "\n")
println("wrote $out_path: $(length(formula_grid)) formula cases, $(length(SIGNAL_CASE_DEFS)) signal cases")
for (name, fc, fmod, spl, dur, fs) in SIGNAL_CASE_DEFS
    text = read(joinpath(signal_output_dir, "$name.json"), String)
    println("  $name: pa = $(parse_json_number_array(text, "pa")), N5 = $(parse_json_number_array(text, "N5")), S5 = $(parse_json_number_array(text, "S5")), R5 = $(parse_json_number_array(text, "R5")), FS5 = $(parse_json_number_array(text, "FS5"))")
end
