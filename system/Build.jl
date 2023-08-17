if VERSION < v"1.9.0-DEV"
  error("Julia >= 1.9 required")
end

oscarpath = pkgdir(Oscar)

println("=====================")
println(oscarpath)
println("=====================")

using Pkg
Pkg.activate(temp=true)
Pkg.add(path="$(oscarpath)")
using Oscar
Pkg.add("PackageCompiler")
Pkg.add("Libdl")

using PackageCompiler, Libdl

tmp = mktempdir(cleanup = false)
CO = joinpath(tmp, "CompileOscar.jl")


write(CO, """
using Pkg
Pkg.add(path="$(oscarpath)")
Pkg.precompile()
using Oscar
using Test
println("oscarpath is $(oscarpath) , while Oscar path is $(pkgdir(Oscar))")
Oscar.system("precompile.jl")
""")

sysimage=joinpath(tmp, "Oscar.$(Libdl.dlext)")
PackageCompiler.create_sysimage([:Oscar], sysimage_path=sysimage, precompile_execution_file=CO)

println("(re)start julia as")
println("\tjulia -J $(sysimage)")
