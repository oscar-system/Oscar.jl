using Pkg
Pkg.add("PackageCompiler")

using PackageCompiler

write("/tmp/CompileOscar.jl", """
using Oscar
using Pkg, Test
include(joinpath(Oscar.pkgdir, "test", "runtests.jl"))
Hecke.test_module("runtests", false)
""")

PackageCompiler.create_sysimage([:Oscar, :Hecke, :Nemo, :AbstractAlgebra, :Singular, :Polymake, :GAP], sysimage_path="/tmp/Oscar.so", precompile_execution_file="/tmp/CompileOscar.jl")        

println("(re)start julia as")
println("\tjulia -J /tmp/Oscar.so")



