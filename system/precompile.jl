import Pkg
Pkg.add("Documenter")
Pkg.add("PrettyTables")
Pkg.add("Aqua")

Pkg.precompile()

include(joinpath(pkgdir(Oscar), "test", "runtests.jl"))
Hecke.system("precompile.jl")
