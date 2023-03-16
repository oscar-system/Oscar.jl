import Pkg
Pkg.add("Documenter")
Pkg.add("PrettyTables")
Pkg.add("Printf")

include(joinpath(pkgdir(Oscar), "test", "runtests.jl"))
Hecke.system("precompile.jl")

