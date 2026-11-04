using Documenter, Oscar

include("make_work.jl")

@invokelatest BuildDoc.doit(Oscar; warnonly=false, local_build=false, doctest=false)
