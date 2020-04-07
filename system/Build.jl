using Pkg
Pkg.add("PackageCompiler")

using PackageCompiler

write("/tmp/CompileOscar.jl", """
using Oscar
using Pkg, Test
include(joinpath(Oscar.pkgdir, "test", "runtests.jl"))
Hecke.test_module("runtests", false)
""")

sysimage="/tmp/Oscar.$(Libdl.dlext)"
PackageCompiler.create_sysimage([:Oscar], sysimage_path=sysimage, precompile_execution_file="/tmp/CompileOscar.jl")

println("(re)start julia as")
println("\tjulia -J $(sysimage)")
