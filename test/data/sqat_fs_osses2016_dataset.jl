# SQAT's FluctuationStrength_Osses2016 (CC BY-NC, github.com/ggrecow/SQAT
# @ 00b449e40599f1c1ef4abe0596094552213d57eb) under Octave 11.3.0
# (octave-signal 1.4.7), method 1 (time-varying, 2 s frames, 90 % overlap),
# run on the REAL Osses et al. (2016) validation stimuli -- the 18 files
# SQAT's own validation scripts read (`AM-tone-fc-1000_fmod-*_mdept-100-SPL-
# 70-dB.wav`, `FM-tone-fc-1500_fmod-*_deltaf-700-SPL-70-dB.wav`,
# `randomnoise-Fc-8010_BW-15980_Fmod-*_Mdept-100_SPL-60.wav`), obtained from
# github.com/aosses-tue/fluctuation-strength-TUe @ c56701c
# (auxdata/osses2016a/Stimuli/; the same files are in doi:10.5281/zenodo.7933206,
# CC BY 4.0). Loaded as in SQAT's validation scripts: samples x 10^((100-94)/20).
# The audio itself is never vendored; these are the oracle's mean fluctuation
# strength [vacil] per file, computed 2026-09-02 via scripts/run_sqat_fs.m.
# Provenance and the reverse-engineered stimulus laws: docs/oracle-pins.md §2.13.
# (kind, fmod [Hz]) => fs_mean
const SQAT_FS_OSSES2016_DATASET = Dict(
    ("am", 1.0) => 0.38178856975167985,
    ("am", 2.0) => 1.1062228530515006,
    ("am", 4.0) => 1.3262007080178402,
    ("am", 8.0) => 1.3129458400020242,
    ("am", 16.0) => 0.33694882539338394,
    ("am", 32.0) => 0.016357179038920276,
    ("fm", 1.0) => 0.880307578393832,
    ("fm", 2.0) => 1.488294062208034,
    ("fm", 4.0) => 2.2068063048644273,
    ("fm", 8.0) => 2.5988149464179013,
    ("fm", 16.0) => 1.2123364800349792,
    ("fm", 32.0) => 0.05445593660564533,
    ("bbn", 1.0) => 0.45650522975417923,
    ("bbn", 2.0) => 1.5573905807738735,
    ("bbn", 4.0) => 1.9035027786891245,
    ("bbn", 8.0) => 1.8192177237189446,
    ("bbn", 16.0) => 0.3722354924377307,
    ("bbn", 32.0) => 0.0348629036991254,
)
