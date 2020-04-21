using Pkg
Pkg.add("PackageCompiler")
Pkg.add("Libdl")

using PackageCompiler, Libdl

write("/tmp/CompileOscar.jl", """
using Oscar
using Pkg, Test
Oscar.system("precompile.jl")
""")

sysimage="/tmp/Oscar.$(Libdl.dlext)"
#PackageCompiler.create_sysimage([:Oscar], sysimage_path=sysimage, precompile_execution_file="/tmp/CompileOscar.jl")
PackageCompiler.create_sysimage([:Oscar], sysimage_path=sysimage, precompile_statements_file="/tmp/Oscar")

println("(re)start julia as")
println("\tjulia -J $(sysimage)")
