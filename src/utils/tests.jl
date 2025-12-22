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
  # default ignore patterns
  ignorepatterns = Regex[
                     # this can only run on the main process and not on distributed workers
                     # so it is included directly in runtests
                     r"Serialization/IPC(\.jl)?$",
                     r"MultigradedImplicitization(\.jl)?$",              
                     # ignore book example code (except for main file)
                     r"(^|/)book/.*/.*\.jl$",
                   ]
  for i in ignore
    if i isa Regex
      push!(ignorepatterns, i)
    elseif i isa AbstractString
      if endswith(i, ".jl")
        push!(ignorepatterns, Regex("$i\$"))
      else
        push!(ignorepatterns, Regex("$i(\\.jl)?\$"))
      end
    else
      throw(ArgumentError("invalid ignore pattern $i"))
    end
  end

  if any(p->contains(path, p), ignorepatterns)
    @info "ignore: $(relpath(path, Oscar.oscardir))"
    return String[]
  end

  if !isabspath(path)
    path = joinpath(Oscar.oscardir, path)
  end

  isfile(path) && return [path]
  isfile("$path.jl") && return ["$path.jl"]

  # if there is a runtests.jl we ignore everything else in that folder
  # except for the main Oscar test dir
  isfile(joinpath(path, "runtests.jl")) &&
    path != joinpath(Oscar.oscardir, "test") &&
    return [joinpath(path, "runtests.jl")]

  tests = String[]
  for entry in readdir(path; join=true)
    if any(s->contains(relpath(entry, Oscar.oscardir), s), ignorepatterns)
      @info "ignore: $(relpath(entry, Oscar.oscardir))"
      continue
    end
    endswith(entry, "setup_tests.jl") && continue
    # this is only for the main test/runtests.jl
    endswith(entry, "runtests.jl") && continue
    if isdir(entry)
      append!(tests, _gather_tests(entry; ignore=ignore))
    elseif isfile(entry) && endswith(entry, ".jl")
      push!(tests, entry)
    end
  end
  return tests
end

@doc raw"""
    Oscar.@_AuxDocTest "Name of the test set", (fix = bool),
    raw\"\"\"
    docstring content here
    \"\"\"

This macro is used to define a doctest from within a testfile.
The first argument is the name of the test set introduced for the doctest.
`bool` may be any expression that evaluates to a boolean value (e.g. `true` or `false`),
and is used to determine whether the doctest should be fixed automatically or not.

This macro is dependent on the correct usage of commas and parentheses. Please refer to the
example above.
"""
macro _AuxDocTest(data::Expr)
  @assert data.head == :tuple
  @assert length(data.args) == 3
  testset = data.args[1]
  @assert data.args[2].head == :(=)
  @assert data.args[2].args[1] == :fix
  fix = data.args[2].args[2]
  docstring = data.args[3]
  @assert docstring.head == :macrocall
  @assert docstring.args[1] == Symbol("@raw_str")
  module_name = Symbol("AuxDocTest_", lstrip(string(gensym()), '#'))
  logging_flag_var = Symbol("logging_flag_", gensym())
  result = Expr(
    :toplevel,
    :(import Documenter),
    :(
      module $(esc(module_name))
      @doc $(docstring) function dummy_placeholder end
      end # module
    ),
    # temporarily disable GC logging to avoid glitches in the doctests
    if isdefined(GC, :logging_enabled)
      esc(
      quote
        $(logging_flag_var) = GC.logging_enabled()
        GC.enable_logging(false)
      end,
    )
    else
      esc(:(GC.enable_logging(false)))
    end,
    esc(
      quote
        Documenter.doctest(
          nothing,
          [$(module_name)];
          fix=$(fix),
          testset=$(testset),
          doctestfilters=[
            Oscar.doctestfilters()...,
            r"(?:^.*Warning: .* is deprecated, use .* instead.\n.*\n.*Core.*\n)?"m,  # removes deprecation warnings
          ],
        )
      end,
    ),
    if isdefined(GC, :logging_enabled)
      esc(
      quote
        GC.enable_logging($(logging_flag_var))
      end,
    )
    else
      esc(:(GC.enable_logging(true)))
    end,
  )
  Meta.replace_sourceloc!(__source__, result)
  return result
end

@doc raw"""
    test_module(path::AbstractString; new::Bool = true, timed::Bool=false, ignore=[])

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

The parameter `ignore` can be used to pass a list of `String` or `Regex` patterns.
Test files or folders matching these will be skipped. Strings will be compared as
suffixes.
This only works for `new=false`.

For experimental modules, use [`test_experimental_module`](@ref) instead.
"""
function test_module(path::AbstractString; new::Bool=true, timed::Bool=false, tempproject::Bool=true, ignore=[])
  with_unicode(false) do
    julia_exe = Base.julia_cmd()
    project_path = Base.active_project()
    if !isabspath(path)
      if !startswith(path, "test")
        path = joinpath("test", path)
      end
      rel_test_path = normpath(path)
      path = joinpath(oscardir, rel_test_path)
    end
    if new
      @req isempty(ignore) && !timed "The `timed` and `ignore` options only work for `new=false`."
      cmd = """
            using Test;
            using Oscar;
            Hecke.assertions(true);
            Oscar.test_module("$path"; new=false);
            """
      @info("spawning ", `$julia_exe --project=$project_path -e \"$cmd\"`)
      run(`$julia_exe --project=$project_path -e $cmd`)
    else
      testlist = _gather_tests(path; ignore=ignore)
      @req !isempty(testlist) "no such file or directory: $path[.jl]"

      use_ctime = timed && VERSION >= v"1.9.0-DEV"
      if use_ctime
        Base.cumulative_compile_timing(true)
      end
      stats = Dict{String,NamedTuple}()

      if tempproject
        # we preserve the old load path
        oldloadpath = copy(LOAD_PATH)
        # make a copy of the test environment to make sure any existing manifest doesn't interfere
        tmpproj = joinpath(mktempdir(), "Project.toml")
        cp(joinpath(Oscar.oscardir, "test", "Project.toml"), tmpproj)
        # activate the temporary project
        Base.set_active_project(tmpproj)
        # and make sure the current project is still available to allow e.g. `using Oscar`
        pushfirst!(LOAD_PATH, dirname(project_path))
        Pkg.resolve()
      else
        @req isdefined(Base.Main, :Test) "You need to do \"using Test\""
      end

      try
        for entry in testlist
          dir = dirname(entry)
          if isfile(joinpath(dir, "setup_tests.jl"))
            Base.include(identity, Main, joinpath(dir, "setup_tests.jl"))
          end
          if timed
            push!(stats, _timed_include(entry; use_ctime=use_ctime))
          else
            Base.include(identity, Main, entry)
          end
        end
      finally
        # restore load path and project
        if tempproject
          copy!(LOAD_PATH, oldloadpath)
          Base.set_active_project(project_path)
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
end

@doc raw"""
    test_experimental_module(project::AbstractString; file::AbstractString="",
      new::Bool=true, timed::Bool=false, ignore=[])

Run the Oscar tests of the experimental module `project`:
- if `file` is empty, run the entire test suite in `experimental/<project>/test`
- if `experimental/<project>/test/<file>` is a directory which contains a `runtests.jl`,
  run this file, otherwise run all test files in that directory and below
- if `experimental/<project>/test/<file>` or `experimental/<project>/test/<file>.jl`
  is a file, run this file

The default is to run the entire test suite of the module `project`.

The optional parameter `new` takes the values `false` and `true` (default). If
`true`, then the tests are run in a new session, otherwise the currently active
session is used.

With the optional parameter `timed` the function will return a dict mapping file
names to a named tuple with compilation times and allocations.
This only works for `new=false`.

The parameter `ignore` can be used to pass a list of `String` or `Regex` patterns.
Test files or folders matching these will be skipped. Strings will be compared as
suffixes.
This only works for `new=false`.
"""
function test_experimental_module(
  project::AbstractString; file::AbstractString="", new::Bool=true, timed::Bool=false, ignore=[]
)
  test_file = "../experimental/$project/test/$file"
  test_module(test_file; new, timed=timed, ignore=ignore)
end

