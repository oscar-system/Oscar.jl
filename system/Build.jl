if VERSION < v"1.9.0-DEV"
  error("Julia >= 1.9 required")
end

oscarpath = pkgdir(Oscar)

using Pkg
Pkg.activate(temp=true)
Pkg.develop(path="$(oscarpath)")
using Oscar
Pkg.add("PackageCompiler")
Pkg.add("Libdl")

using PackageCompiler, Libdl

tmp = mktempdir(cleanup = false)
CO = joinpath(tmp, "CompileOscar.jl")


write(CO, """
using Pkg
Pkg.develop(path="$(oscarpath)")
using Oscar
using Test
Oscar.system("precompile.jl")
""")

sysimage=joinpath(tmp, "Oscar.$(Libdl.dlext)")
if !("JULIA_CPU_TARGET" in keys(ENV)) || (ENV["JULIA_CPU_TARGET"] == "")
  println("Building system image for generic target. Use JULIA_CPU_TARGET to change.")
  target = PackageCompiler.default_app_cpu_target()
else
  target = ENV["JULIA_CPU_TARGET"]
end
PackageCompiler.create_sysimage([:Oscar], sysimage_path=sysimage, precompile_execution_file=CO; cpu_target=target)

println("(re)start julia as")
println("\tjulia -J $(sysimage)")
