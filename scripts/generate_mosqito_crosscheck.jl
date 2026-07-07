# One-time generator for test/data/mosqito_crosscheck.jl. Requires `uv`
# (https://docs.astral.sh/uv/) and network access on first run.
#
# Usage: julia scripts/generate_mosqito_crosscheck.jl

using Pkg
Pkg.activate(; temp=true)
Pkg.add(Pkg.PackageSpec(name="ZwickerLoudness", version="0.2"))

using ZwickerLoudness

include(joinpath(@__DIR__, "..", "test", "data", "crosscheck_cases.jl"))

jsonvec(v) = string("[", join(string.(v), ","), "]")

tmp = tempname() * ".json"
open(tmp, "w") do io
    entries = String[]
    for (name, spl) in CROSSCHECK_CASES
        r = zwicker_loudness(spl)
        push!(entries, """{"name": "$name", "specific_loudness": $(jsonvec(r.specific_loudness))}""")
    end
    print(io, "[", join(entries, ","), "]")
end

out = joinpath(@__DIR__, "..", "test", "data", "mosqito_crosscheck.jl")
run(`uv run --with mosqito --with numpy --with matplotlib python $(joinpath(@__DIR__, "crosscheck_mosqito.py")) $tmp $out`)
