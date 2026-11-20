using Documenter, Oscar

include("make_work.jl")

formats = "--pdf" in ARGS ? [:pdf] : [:html]
@invokelatest BuildDoc.doit(Oscar; warnonly=false, local_build=false, doctest=false, formats)
