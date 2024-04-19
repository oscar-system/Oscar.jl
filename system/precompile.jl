import Pkg
Pkg.add("Documenter")
Pkg.add("PrettyTables")
Pkg.add("Aqua")
Pkg.add("DeepDiffs")

Pkg.precompile()

include(joinpath(pkgdir(Oscar), "test", "runtests.jl"))
Hecke.system("precompile.jl")
