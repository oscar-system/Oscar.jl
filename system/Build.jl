using Pkg
Pkg.add("PackageCompiler")
Pkg.add("Libdl")

using PackageCompiler, Libdl

tmp = mktempdir(cleanup = false)
CO = joinpath(tmp, "CompileOscar.jl")

write(CO, """
using Oscar
using Pkg, Test
Oscar.system("precompile.jl")
""")

sysimage=joinpath(tmp, "Oscar.$(Libdl.dlext)")
PackageCompiler.create_sysimage([:Oscar], sysimage_path=sysimage, precompile_execution_file=CO)

println("(re)start julia as")
println("\tjulia -J $(sysimage)")
