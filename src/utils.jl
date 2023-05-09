###############################################################################
###############################################################################
##
##  versioninfo
##
###############################################################################
###############################################################################
# When a specific branch is loaded via `]add Package#branch` julia will only
# create a checkout and keep a bare git repo in a separate directory.
# In a bare repo HEAD will not point to the correct commit so we use the git
# tree-hash that Pkg.jl provides and manually map this to a corresponding
# commit.
function _lookup_commit_from_cache(url::AbstractString, tree::AbstractString)
   if Sys.which("git") != nothing
      try
         path = Pkg.Types.add_repo_cache_path(url)
         if isdir(path)
            commit = readchomp(`sh -c "git -C $path log --oneline --all --pretty='tree %T;%H' | grep \"^tree $tree\" | cut -d\; -f2 | head -n1"`)
            return readchomp(`git -C $path show -s --format=", %h -- %ci" $commit`)
         end
      catch
      end
   end
   return ""
end

function _lookup_git_branch(dir::AbstractString; commit=false)
   info = ""
   if Sys.which("git") != nothing &&
         isdir(joinpath(dir,".git"))
      try
         ref = readchomp(`git -C $dir rev-parse --abbrev-ref HEAD`)
         info = " - #$(ref)"
         if commit
            c = readchomp(`git -C $dir show -s --format="%h -- %ci" HEAD`)
            info = "$info, $c"
         end
      catch
      end
   end
   return info
end

function _deps_git_info(dep::Pkg.API.PackageInfo; commit=false)
   if dep.is_tracking_repo
      info = commit ? _lookup_commit_from_cache(dep.git_source, dep.tree_hash) : ""
      return " - #$(dep.git_revision)$info"
   elseif dep.is_tracking_path
      return _lookup_git_branch(dep.source; commit=commit)
   end
   return ""
end

function _print_dependency_versions(io::IO, deps::AbstractArray{<:AbstractString}; padding="    ", suffix="", branch=false, commit=false)
   width = maximum(length.(deps))+length(suffix)+2
   deps = filter(d->d.name in deps, collect(values(Pkg.dependencies())))
   deps = sort!(deps; by=x->x.name)
   for dep in deps
      print(io, "$(padding)$(rpad(dep.name*suffix, width, ' ')) v$(dep.version)")
      println(io, branch ? _deps_git_info(dep; commit=commit) : "")
   end
end

@doc raw"""
    Oscar.versioninfo(io::IO=stdout; branch=false, jll=false, julia=false, commit=false, full=false)

Print the versions of all Oscar-related dependencies.

# Arguments
- `branch::Bool=false`: include git branch name in the output
- `commit::Bool=false`: include git commit hash and date where applicable
- `jll::Bool=false`   : include binary packages (jll) in the output
- `julia::Bool=false` : include julia `versioninfo` output
- `full::Bool=false`  : include all of the above
"""
function versioninfo(io::IO=stdout; branch=false, jll=false, julia=false, commit=false, full=false)
   if full
      branch = jll = julia = commit = true
   end
   print(io, "OSCAR version $(VERSION_NUMBER)")
   println(io, branch ? _lookup_git_branch(dirname(@__DIR__); commit=commit) : "")
   println(io, "  combining:")
   _print_dependency_versions(io, cornerstones; suffix=".jl", branch=branch, commit=commit)
   if jll
      println(io, "  building on:")
      _print_dependency_versions(io, jll_deps; branch=branch, commit=commit)
      println(io, "See `]st -m` for a full list of dependencies.")
   end
   if julia
      println(io, "")
      Main.InteractiveUtils.versioninfo(io)
      println(io, Base.TAGGED_RELEASE_BANNER)
   end
end

###############################################################################
###############################################################################
##
##  Documentation helpers
##
###############################################################################
###############################################################################

# use tempdir by default to ensure a clean manifest (and avoid modifying the project)
function doc_init(;path=mktempdir())
  global docsproject = path
  if !isfile(joinpath(docsproject,"Project.toml"))
    cp(joinpath(oscardir, "docs", "Project.toml"), joinpath(docsproject,"Project.toml"))
  end
  Pkg.activate(docsproject) do
    # we dev all packages with the paths from where they are currently loaded
    for dir in [aadir, nemodir, heckedir, oscardir]
      Pkg.develop(path=dir)
    end
    Pkg.instantiate()
    Base.include(Main, joinpath(oscardir, "docs", "make_work.jl"))
  end
end

#function doc_update_deps()
#  Pkg.activate(Pkg.update, joinpath(oscardir, "docs"))
#end

function open_doc()
    filename = normpath(Oscar.oscardir, "docs", "build", "index.html")
    @static if Sys.isapple()
        run(`open $(filename)`; wait = false)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`; wait = false)
    elseif Sys.iswindows()
        cmd = get(ENV, "COMSPEC", "cmd.exe")
        run(`$(cmd) /c start $(filename)`; wait = false)
    else
        @warn("Opening files the default application is not supported on this OS.",
              KERNEL = Sys.KERNEL)
    end
end


@doc raw"""
    build_doc(; doctest=false, strict=false, open_browser=true)

Build the manual of `Oscar.jl` locally and open the front page in a
browser.

The optional parameter `doctest` can take three values:
  - `false`: Do not run the doctests (default).
  - `true`: Run the doctests and report errors.
  - `:fix`: Run the doctests and replace the output in the manual with
    the output produced by Oscar. Please use this option carefully.

In GitHub Actions the Julia version used for building the manual is 1.8 and
doctests are run with >= 1.7. Using a different Julia version may produce
errors in some parts of Oscar, so please be careful, especially when setting
`doctest=:fix`.

The optional parameter `strict` is passed on to `makedocs` of `Documenter.jl`
and if set to `true` then according to the manual of `Documenter.jl` "a
doctesting error will always make makedocs throw an error in this mode".

To prevent the opening of the browser at the end, set the optional parameter
`open_browser` to `false`.

When working on the manual the `Revise` package can significantly sped
up running `build_doc`. First, install `Revise` in the following way:
```
using Pkg ; Pkg.add("Revise")
```
Second, restart Julia and load `Revise` before Oscar:
```
using Revise, Oscar;
```
The first run of `build_doc` will take the usual few minutes, subsequently runs
will be significantly faster.
"""
function build_doc(; doctest=false, strict=false, open_browser=true)
  versioncheck = (VERSION.major == 1) && (VERSION.minor >= 7)
  versionwarn = 
"The Julia reference version for the doctests is 1.7 or later, but you are using
$(VERSION). Running the doctests will produce errors that you do not expect."
  if doctest != false && !versioncheck
    @warn versionwarn
  end
  if !isdefined(Main, :BuildDoc)
    doc_init()
  end
  Pkg.activate(docsproject) do
    Base.invokelatest(Main.BuildDoc.doit, Oscar; strict=strict, local_build=true, doctest=doctest)
  end
  if open_browser
    open_doc()
  end
  if doctest != false && !versioncheck
    @warn versionwarn
  end
end

###############################################################################
###############################################################################
##
##  Testsuite helpers
##
###############################################################################
###############################################################################
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
  rel_test_file = normpath("test", "$file.jl")
  test_file = joinpath(oscardir, rel_test_file)

  if new
    cmd = "using Test; using Oscar; Hecke.assertions(true); include(\"$test_file\");"
    @info("spawning ", `$julia_exe -e \"$cmd\"`)
    run(`$julia_exe -e $cmd`)
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
