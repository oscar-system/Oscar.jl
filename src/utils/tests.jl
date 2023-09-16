function _timed_include(str::String, mod::Module=Main; use_ctime::Bool=VERSION >= v"1.9.0")
  if use_ctime
    compile_elapsedtimes = Base.cumulative_compile_time_ns()
  end
  stats = @timed Base.include(identity, mod, str)
  fullpath = abspath(joinpath(Base.source_dir(), str))
  # skip files which just include other files and ignore
  # files outside of the oscar folder
  if startswith(fullpath, Oscar.oscardir)
    path = relpath(fullpath, Oscar.oscardir)
    if use_ctime
      compile_elapsedtimes = Base.cumulative_compile_time_ns() .- compile_elapsedtimes
      compile_elapsedtimes = compile_elapsedtimes ./ 10^9
    end
    rtime=NaN
    if use_ctime
      comptime = first(compile_elapsedtimes)
      rcomptime = last(compile_elapsedtimes)
      println("-> Testing $path took: runtime $(round(stats.time-comptime; digits=3)) seconds + compilation $(round(comptime-rcomptime; digits=3)) seconds + recompilation $(round(rcomptime; digits=3)) seconds, $(Base.format_bytes(stats.bytes))")
      return (path=>(time=stats.time-comptime, ctime=comptime-rcomptime, rctime=rcomptime, alloc=stats.bytes/2^20))
    else
      println("-> Testing $path took: $(round(stats.time; digits=3)) seconds, $(Base.format_bytes(stats.bytes))")
      return (path=>(time=stats.time, alloc=stats.bytes/2^20))
    end
  else
    return ()
  end
end

function _gather_tests(path::AbstractString; ignore=[])
  if !isabspath(path)
    path = joinpath(Oscar.oscardir, path)
  end
  # default ignore dirs, these are compared as suffix with endswith
  # i.e. make sure they are unique
  ignorepaths = [
                  # these two files seem obsolete
                  "Modules/GradedModules.jl",
                  "Modules/FreeModules-graded.jl",
                  # FIXME: temporarily disable AlgClosureFp tests until we resolve
                  # issue https://github.com/oscar-system/Oscar.jl/issues/2691
                  "Rings/AlgClosure.jl",
                ]
  append!(ignorepaths, ignore)

  any(s->endswith(path, s) ,ignorepaths) && return String[]
  isfile(path) && return [path]

  # if there is a runtests.jl we ignore everything else in that folder
  # except for the main Oscar test dir
  isfile(joinpath(path, "runtests.jl")) &&
    path != joinpath(Oscar.oscardir, "test") &&
    return [joinpath(path, "runtests.jl")]

  tests = String[]
  for entry in readdir(path; join=true)
    any(s->endswith(entry, s), ignorepaths) && continue
    if isdir(entry)
      append!(tests, _gather_tests(entry; ignore=ignore))
    elseif isfile(entry) && endswith(entry, ".jl") &&
        !endswith(entry, "setup_tests.jl")
      push!(tests, entry)
    end
  end
  return tests
end




@doc raw"""
    test_module(path::AbstractString; new::Bool = true, timed::Bool=false)

Run the Oscar tests in `path`:
- if `path` is relative then it will be set to `<oscardir>/test/<path>`
- if `path` is a directory, run all test files in that directory and below
- if `path` or `path.jl` is a file, run this file

If a directory contains a `runtests.jl` file only this file will be executed, otherwise
all files will be included independently.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.

With the optional parameter `timed` the function will return a dict mapping file
names to a named tuple with compilation times and allocations.
This only works for `new=false`.

For experimental modules, use [`test_experimental_module`](@ref) instead.
"""
function test_module(path::AbstractString; new::Bool=true, timed::Bool=false)
  julia_exe = Base.julia_cmd()
  project_path = Base.active_project()
  if !isabspath(path)
    rel_test_path = normpath("test", "$path")
    path = joinpath(oscardir, rel_test_path)
  end
  if new
    cmd = "using Test; using Oscar; Hecke.assertions(true); Oscar.test_module(\"$entry\"; new=false);"
    @info("spawning ", `$julia_exe --project=$project_path -e \"$cmd\"`)
    run(`$julia_exe --project=$project_path -e $cmd`)
  else
    if isfile(path)
      testlist = [path]
    elseif isfile("$path.jl")
      testlist = ["$path.jl"]
    elseif isdir(path)
      testlist = _gather_tests(path)
    else
      @req false "no such file or directory: $path[.jl]"
    end

    @req isdefined(Base.Main, :Test) "You need to do \"using Test\""

    use_ctime = timed && VERSION >= v"1.9.0-DEV"
    if use_ctime
      Base.cumulative_compile_timing(true)
    end
    stats = Dict{String,NamedTuple}()
    for entry in testlist
      dir = dirname(entry)
      if isfile(joinpath(dir,"setup_tests.jl"))
        Base.include(identity, Main, joinpath(dir,"setup_tests.jl"))
      end
      if timed
        push!(stats, _timed_include(entry; use_ctime=use_ctime))
      else
        Base.include(identity, Main, entry)
      end
    end

    if timed
      use_ctime && Base.cumulative_compile_timing(false)
      return stats
    else
      return nothing
    end
  end
end

@doc raw"""
    test_experimental_module(project::AbstractString; file::AbstractString="",
      new::Bool=true, timed::Bool=false)

Run the Oscar tests in `experimental/<project>/test/<path>`:
- if `path` is empty then all tests in that module are run, either via `runtests.jl` or directly.
- if `path` or `path.jl` is a file in that directory only this file is run.

The default is to run the entire test suite of the module `project`.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.

With the optional parameter `timed` the function will return a dict mapping file
names to a named tuple with compilation times and allocations.

"""
function test_experimental_module(
  project::AbstractString; file::AbstractString="", new::Bool=true, timed::Bool=false
)
  test_file = "../experimental/$project/test/$file"
  test_module(test_file; new, timed=timed)
end

