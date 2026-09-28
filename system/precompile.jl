import Pkg
Pkg.precompile()
ENV["OSCAR_TEST_SUBSET"] = "short"
include(joinpath(pkgdir(Oscar), "test", "runtests.jl"))
Hecke.system("precompile.jl")
