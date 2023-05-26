@doc raw"""
    test_module(file::AbstractString; new::Bool = true)

Run the Oscar tests in the file `test/<file>.jl` where `file` may be a path.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.

For experimental modules, use [`test_experimental_module`](@ref) instead.
"""
function test_module(file::AbstractString; new::Bool=true)
  julia_exe = Base.julia_cmd()
  project_path = Base.active_project()
  rel_test_file = normpath("test", "$file.jl")
  test_file = joinpath(oscardir, rel_test_file)

  if new
    cmd = "using Test; using Oscar; Hecke.assertions(true); include(\"$test_file\");"
    @info("spawning ", `$julia_exe --project=$project_path -e \"$cmd\"`)
    run(`$julia_exe --project=$project_path -e $cmd`)
  else
    @req isdefined(Base.Main, :Test) "You need to do \"using Test\""
    @info("Running tests for $rel_test_file in same session")
    Base.include(Base.Main, test_file)
  end
end

@doc raw"""
    test_experimental_module(project::AbstractString; file::AbstractString="runtests", new::Bool = true)

Run the Oscar tests in the file `experimental/<project>/test/<file>.jl`
where `file` may be a path.
The default is to run the entire test suite of the module `project`.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.
"""
function test_experimental_module(
  project::AbstractString; file::AbstractString="runtests", new::Bool=true
)
  test_file = "../experimental/$project/test/$file"
  test_module(test_file; new)
end

