# Input spectra (28-band one-third-octave SPL [dB]) for the MoSQITo
# weighting cross-check. Shared by the test and by
# scripts/generate_mosqito_crosscheck.jl.
const CROSSCHECK_CASES = [
    ("readme_example", [60.0, 62, 65, 68, 70, 72, 74, 75, 73, 71,
                        69, 67, 65, 63, 61, 59, 57, 55, 53, 50,
                        47, 44, 41, 38, 35, 32, 29, 26]),
    ("rising_ramp",  collect(range(30.0, 75.0; length=28))),
    ("falling_ramp", collect(range(75.0, 30.0; length=28))),
]
