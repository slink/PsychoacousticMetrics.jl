# Formula cross-check vs SQAT's PsychoacousticAnnoyance_Widmann1992_from_percentile
# (CC BY-NC, github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb),
# oracle only — no code transcribed (see .superpowers/sdd/pa-oracle-pins.md).
# Pure scalar arithmetic: any disagreement beyond the measured tolerance is a
# transcription bug, not a tolerance question.

isdefined(@__MODULE__, :SQAT_PA_FORMULA_CASES) ||
    include(joinpath(@__DIR__, "data", "sqat_pa_crosscheck.jl"))

@testset "psychoacoustic_annoyance_widmann formula" begin
    @testset "cross-check vs SQAT (840-case grid)" begin
        # Measured max relative deviation across all 840 cases, this
        # machine: 0.0 exactly (identical bit patterns) — see
        # .superpowers/sdd/pa-task-2-report.md. rtol kept at 1e-12 (pure
        # scalar arithmetic; exceeding 1e-9 would indicate a transcription
        # bug, not a tolerance question).
        RTOL = 1e-12
        for c in SQAT_PA_FORMULA_CASES
            pa = psychoacoustic_annoyance_widmann(c.N, c.S, c.R, c.FS)
            @test isapprox(pa, c.pa; rtol = RTOL, atol = 1e-12)
        end
    end

    @testset "N = 0 => PA = 0 exactly" begin
        # Pinned (pa-oracle-pins.md Step 4): all 120 grid rows with N == 0.0
        # give pa == 0.0 exactly via the reference's Inf/NaN zeroing,
        # including the R = FS = 0 sub-case where wfr = Inf*0 = NaN before
        # zeroing.
        zero_rows = filter(c -> c.N == 0.0, SQAT_PA_FORMULA_CASES)
        @test length(zero_rows) == 120
        for c in zero_rows
            @test psychoacoustic_annoyance_widmann(c.N, c.S, c.R, c.FS) == 0.0
        end
        @test psychoacoustic_annoyance_widmann(0.0, 2.5, 0.3, 0.3) == 0.0
        @test psychoacoustic_annoyance_widmann(0.0, 2.5, 0.0, 0.0) == 0.0
        @test psychoacoustic_annoyance_widmann(0.0, 0.5, 0.0, 0.0) == 0.0
    end

    @testset "threshold edge: strict S > 1.75" begin
        # Pinned (pa-oracle-pins.md Step 4): S = 1.75 takes the ws = 0
        # branch (else); only S > 1.75 (e.g. 1.7500001) engages ws != 0.
        pa_175 = psychoacoustic_annoyance_widmann(4.0, 1.75, 0.3, 0.3)
        pa_below = psychoacoustic_annoyance_widmann(4.0, 1.0, 0.3, 0.3)
        pa_above = psychoacoustic_annoyance_widmann(4.0, 1.7500001, 0.3, 0.3)
        @test pa_175 == pa_below
        @test pa_175 != pa_above
        @test isapprox(pa_175, 5.502497448336122; rtol = 1e-12)
        @test isapprox(pa_above, 5.502497448336126; rtol = 1e-12)
    end

    @testset "tiny N stays finite" begin
        # Pinned (pa-oracle-pins.md Step 4): no epsilon guard needed beyond
        # the documented Inf/NaN zeroing.
        pa = psychoacoustic_annoyance_widmann(1e-12, 2.5, 0.3, 0.3)
        @test isfinite(pa)
        @test pa > 0
        pa_zero_rfs = psychoacoustic_annoyance_widmann(1e-12, 2.5, 0.0, 0.0)
        @test isfinite(pa_zero_rfs)
    end

    @testset "input validation" begin
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(-1.0, 1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, -1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, -0.1, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, 0.0, -0.1)

        @test_throws ArgumentError psychoacoustic_annoyance_widmann(NaN, 1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(Inf, 1.0, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, NaN, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, Inf, 0.0, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, NaN, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, Inf, 0.0)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, 0.0, NaN)
        @test_throws ArgumentError psychoacoustic_annoyance_widmann(1.0, 1.0, 0.0, Inf)
    end
end
