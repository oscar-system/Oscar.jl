using Pkg
if VERSION >= v"1.4"
  d = Pkg.dependencies()
  f = filter(x->occursin("PackageCompiler", x.name), collect(values(d)))
else
  d = Pkg.installed()
  f = filter(x->occursin("PackageCompiler", x), keys(d))
end

if length(f) == 0
  Pkg.add("PackageCompiler")
end

using PackageCompiler

f = open("/tmp/CompileOscar.jl", "w")
println(f, "using Oscar")
println(f, "using Pkg, Test")
println(f, "include(joinpath(Oscar.pkgdir, \"test\", \"runtests.jl\"))")
println(f, "Hecke.test_module(\"runtests\", false)")
close(f)

PackageCompiler.create_sysimage([:Oscar, :Hecke, :Nemo, :AbstractAlgebra, :Singular, :Polymake, :GAP], sysimage_path="/tmp/Oscar.so", precompile_execution_file="/tmp/CompileOscar.jl")        

println("(re)start julia as")
println("\tjulia -J /tmp/Oscar.so")



