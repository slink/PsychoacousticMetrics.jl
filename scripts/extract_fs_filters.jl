# One-time filter-coefficient extraction for the fluctuation-strength
# oracle rig. Reads two kinds of FACTS out of the pinned SQAT reference
# (github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb,
# CC BY-NC — numeric values only are extracted here, never source code):
#
#   1. The `Hweight` band-pass second-order-section filters used to
#      estimate modulation depth (private/Hweight-{fs}-Hz-{HP,LP}.mat),
#      read via scipy.io.loadmat — these are plain MATLAB v5 .mat data
#      files, not code.
#   2. The `a0` (outer/middle-ear transmission) FIR filter taps, produced
#      by calling the reference's own `calculate_a0(fs, 4096,
#      'fluctuationstrength_osses2016')` under Octave and dumping the
#      returned coefficient vector — a numeric artifact, not a
#      transcription of the function's logic.
#
# Requires: the pinned checkout at /tmp/sqat-pinned (see
# .superpowers/sdd/fs-oracle-pins.md Step 1), `uv` (for scipy), and Octave
# with the octave-signal package installed.
#
# Usage: julia scripts/extract_fs_filters.jl
using Dates: now

const SQAT_DIR = get(ENV, "SQAT_DIR", "/tmp/sqat-pinned")
const PRIVATE_DIR = joinpath(SQAT_DIR, "psychoacoustic_metrics", "FluctuationStrength_Osses2016", "private")
const SCRATCH_DIR = joinpath(@__DIR__, "..", ".superpowers", "sdd")

isdir(SQAT_DIR) || error("pinned SQAT checkout not found at $SQAT_DIR — clone+checkout per Step 1 first")
mkpath(SCRATCH_DIR)

octave_version = strip(read(`octave --eval "disp(OCTAVE_VERSION)"`, String))
signal_version = strip(read(`octave --eval "pkg load signal; disp(ver('signal').Version)"`, String))
provenance = "SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb " *
    "(github.com/ggrecow/SQAT, CC BY-NC — values only, no code copied), " *
    "Octave $octave_version, octave-signal $signal_version, $(now())"
header = """
# Extracted $(now()) from SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb
# (github.com/ggrecow/SQAT, CC BY-NC — values only, no code copied).
# Octave $(octave_version), octave-signal $(signal_version).
"""

# --- 1. Hweight SOS matrices, via scipy (values only) -----------------------
py = """
import json
import scipy.io as sio

out = {"_provenance": "$provenance"}
for fs in (44100, 48000):
    hp = sio.loadmat("$(PRIVATE_DIR)/Hweight-%d-Hz-HP.mat" % fs)["Hweight_HP"]
    lp = sio.loadmat("$(PRIVATE_DIR)/Hweight-%d-Hz-LP.mat" % fs)["Hweight_LP"]
    assert hp.shape == (1, 6), "HP SOS must be 1 section (2nd order): got %r" % (hp.shape,)
    assert lp.shape == (2, 6), "LP SOS must be 2 sections (4th order): got %r" % (lp.shape,)
    assert (hp[:, 3] == 1).all() and (lp[:, 3] == 1).all(), "a0 (column 4) must be 1 in every SOS row"
    out["Hweight_HP_%d" % fs] = hp.tolist()
    out["Hweight_LP_%d" % fs] = lp.tolist()

with open("$(joinpath(SCRATCH_DIR, "hweight_sos.json"))", "w") as f:
    json.dump(out, f, indent=2)

def julia_lit(name, rows):
    body = "; ".join(" ".join("%.17g" % v for v in row) for row in rows)
    print("const %s = [%s]" % (name, body))

for key, rows in out.items():
    if not key.startswith("_"):
        julia_lit(key, rows)
print("SANITY_OK: Hweight SOS shapes and a0==1 verified for fs in {44100, 48000}")
"""
hweight_out = read(`uv run --with scipy python3 -c $py`, String)
occursin("SANITY_OK", hweight_out) || error("Hweight extraction sanity checks failed:\n$hweight_out")
println(hweight_out)

# --- 2. a0 FIR taps, via Octave's own calculate_a0 (values only) -----------
mscript = """
pkg load signal;
addpath('$(joinpath(SQAT_DIR, "utilities"))');
fid_json = fopen('$(joinpath(SCRATCH_DIR, "a0_fir_b.json"))', 'w');
fprintf(fid_json, '{\\n  \"_provenance\": \"$provenance\",\\n');
fs_list = [44100 48000];
for i = 1:numel(fs_list)
  fs = fs_list(i);
  B = calculate_a0(fs, 4096, 'fluctuationstrength_osses2016');
  assert(length(B) == 4097, 'expected 4097 taps, got %d', length(B));
  symerr = max(abs(B - fliplr(B)));
  assert(symerr < 1e-9, 'B not symmetric to 1e-9: max err %.3e', symerr);
  printf('const A0_FIR_B_%d = [', fs);
  printf('%.17g', B(1));
  printf(', %.17g', B(2:end));
  printf(']\\n');
  fprintf(fid_json, '  \"%d\": [', fs);
  fprintf(fid_json, '%.17g', B(1));
  fprintf(fid_json, ', %.17g', B(2:end));
  if i < numel(fs_list)
    fprintf(fid_json, '],\\n');
  else
    fprintf(fid_json, ']\\n');
  end
  printf('SANITY_OK: fs=%d length(B)=%d symmetry_err=%.3e\\n', fs, length(B), symerr);
end
fprintf(fid_json, '}\\n');
fclose(fid_json);
"""
a0_out = read(`octave --no-gui --eval $mscript`, String)
occursin("SANITY_OK", a0_out) || error("a0 FIR extraction sanity checks failed:\n$a0_out")
println(a0_out)

println(header)
println("Wrote $(joinpath(SCRATCH_DIR, "hweight_sos.json")) and $(joinpath(SCRATCH_DIR, "a0_fir_b.json"))")
