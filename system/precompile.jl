import Pkg
Pkg.add("Documenter")
Pkg.add("PrettyTables")
Pkg.add("Printf")
Pkg.add("Aqua")

Pkg.precompile()

println("Path of Oscar is $(pkgdir(Oscar))")

exit(1)

include(joinpath(pkgdir(Oscar), "test", "runtests.jl"))
Hecke.system("precompile.jl")

