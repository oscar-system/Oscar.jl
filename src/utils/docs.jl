###############################################################################
###############################################################################
##
##  Documentation helpers
##
###############################################################################
###############################################################################

################################################################################
#
#  DocTestSetup
#
################################################################################

# Oscar needs some complicated setup to get the printing right. This provides a
# helper function to set this up consistently.
function doctestsetup()
  link_experimental_docs()
  return :(using Oscar; Oscar.AbstractAlgebra.set_current_module(@__MODULE__))
end

# The settings every documentation page needs. Passed to Documenter as the
# default for its `@meta` blocks, so that no page has to repeat them.
function docmeta()
  return Dict(
    :CurrentModule => Oscar,
    :CollapsedDocStrings => true,
    :DocTestSetup => doctestsetup(),
  )
end

function doctestfilters()
  # this returns a list of doctest filters that should be passed to all doctest invocations
  return [
    # due to JuliaLang/julia#57898, oscar-system/Oscar.jl#4872:
    r"^0-element SubObjectIterator{RayVector{QQFieldElem}}|^RayVector{QQFieldElem}\[\]"m => "RayVector{QQFieldElem}[]",
    r"^0-element SubObjectIterator{PointVector{QQFieldElem}}|^PointVector{QQFieldElem}\[\]"m => "PointVector{QQFieldElem}[]",
    r"^ 0-element SubObjectIterator{RayVector{QQFieldElem}}|^ \[\]"m => " []",
    r"^ 0-element SubObjectIterator{PointVector{QQFieldElem}}|^ \[\]"m => " []",
  ]
end

# https://github.com/JuliaLang/julia/pull/57509 changed hashing for many base types.
# We thus filter doctestoutputs that show Sets, Dicts, etc. for julia versions with
# this change to still have the doctests passing.
function doctestfilter_hash_changes_in_1_13()
  if VERSION >= v"1.13.0-DEV.570"
    # The filter in place here checks that the number of lines in the output matches the input.
    [r".*"]
  else
    return []
  end
end

################################################################################
#
#  The documentation environment
#
################################################################################
#
# Building the manual needs Documenter, which OSCAR does not depend on. `docs/`
# is an environment providing it, stacked onto whatever project the user works
# in, so that OSCAR itself keeps being loaded from there.
#
# It deliberately contains no OSCAR package. That keeps its manifest independent
# of OSCAR's dependency graph, so it neither has to be re-resolved whenever that
# graph moves nor can it disagree with the OSCAR that is already loaded.

# Instantiating is cheap once done, but not free; only do it once per session.
const _docs_env_instantiated = Ref(false)

@doc raw"""
    Oscar.docs_env()

Make `Documenter` and the `OscarDocs` build driver loadable by appending the
`docs/` environment to `LOAD_PATH`, instantiating it if necessary, and return
its path.

This is called for you by [`Oscar.@build_doc`](@ref) and friends.
"""
function docs_env()
  env = joinpath(oscardir, "docs")

  if !_docs_env_instantiated[]
    Pkg.activate(Pkg.instantiate, env)
    _docs_env_instantiated[] = true
  end

  env in LOAD_PATH || push!(LOAD_PATH, env)
  return env
end

################################################################################
#
#  Entry points
#
################################################################################
#
# Each entry point has to set up the environment, load `OscarDocs` and then call
# into it. The macros do that as a `:toplevel` block evaluated in the caller's
# module, so that each step gets its own world age. The functions of the same
# names do the same work but have to resort to `invokelatest` for it.

function _split_macro_args(args)
  kws = Any[]
  pos = Any[]
  for arg in args
    if arg isa Expr && (arg.head === :(=) || arg.head === :kw)
      push!(kws, Expr(:kw, arg.args[1], arg.args[2]))
    else
      push!(pos, arg)
    end
  end
  return pos, kws
end

function _docs_toplevel(__source__, callee::Symbol, pos, kws)
  # `Oscar` goes in as the module object itself, so that the caller does not
  # need to have it in scope
  call = Expr(:call, Expr(:., :OscarDocs, QuoteNode(callee)),
              Expr(:parameters, kws...), Oscar, pos...)
  result = Expr(:toplevel,
                :($(Oscar).docs_env()),
                :(import OscarDocs),
                esc(call))
  Meta.replace_sourceloc!(__source__, result)
  return result
end

@doc raw"""
    Oscar.@build_doc
    Oscar.@build_doc doctest=false warnonly=true
    Oscar.@build_doc open_browser=true start_server=false

Build the manual of `Oscar.jl` locally and open the front page in a browser.

The optional parameter `doctest` can take three values:
  - `false`: Do not run the doctests (default).
  - `true`: Run the doctests and report errors.
  - `:fix`: Run the doctests and replace the output in the manual with
    the output produced by Oscar. Please use this option carefully.

The optional parameter `warnonly` is passed on to `makedocs` of `Documenter.jl`
and if set to `false` then according to the manual of `Documenter.jl` "a
doctesting error will always make makedocs throw an error in this mode".
Alternatively, one can pass a list of symbols to `warnonly` to suppress
errors for the given error types.

To prevent the opening of the browser at the end, set the optional parameter
`open_browser` to `false`.

Alternatively, one can use the optional parameter `start_server` to start a web
server in the build directory which is accessible via `127.0.0.1:8000`. If both
`start_server` and `open_browser` are `true`, the browser will show this page.
The server keeps running in the background until the `julia` session is
terminated, so the proposed usage for this option is to run
`Oscar.@build_doc start_server=true` for the first build and
`Oscar.@build_doc open_browser=false` for following builds and only refresh the
browser tab.

When working on the manual the `Revise` package can significantly speed
up rebuilding. First, install `Revise` in the following way:
```julia
using Pkg ; Pkg.add("Revise")
```
Second, restart Julia and load `Revise` before Oscar:
```julia
using Revise, Oscar;
```
The first build will take the usual few minutes, subsequent ones will be
significantly faster.

!!! note
    This must be called at the top level, it does not work inside a function.
"""
macro build_doc(args...)
  pos, kws = _split_macro_args(args)
  isempty(pos) || throw(ArgumentError("@build_doc only takes keyword arguments"))
  return _docs_toplevel(__source__, :build, pos, kws)
end

@doc raw"""
    Oscar.@doctest
    Oscar.@doctest f
    Oscar.@doctest path

Run the doctests of OSCAR. Without an argument, all of them are run, both those
in docstrings and those in the manual sources; this does not build the manual.

Given a function `f`, only the doctests of that function are run. Given a
string `path`, only the doctests in files whose full pathname contains `path`
are run.

# Examples
```julia
julia> Oscar.@doctest symmetric_group

julia> Oscar.@doctest "/Rings/"
```

!!! note
    This must be called at the top level, it does not work inside a function.
"""
macro doctest(args...)
  pos, kws = _split_macro_args(args)
  length(pos) <= 1 || throw(ArgumentError("@doctest takes at most one argument"))
  return _docs_toplevel(__source__, :doctest, pos, kws)
end

@doc raw"""
    Oscar.@doctest_fix
    Oscar.@doctest_fix f
    Oscar.@doctest_fix path

Like [`Oscar.@doctest`](@ref), but replace the expected output of every doctest
that is run by the output OSCAR actually produces.

# Examples
The following call fixes all doctests for the function `symmetric_group`:
```julia
julia> Oscar.@doctest_fix symmetric_group
```
The following call fixes all doctests in files that live in a directory called
`Rings` (or a subdirectory thereof), so e.g. everything in `src/Rings/`:
```julia
julia> Oscar.@doctest_fix "/Rings/"
```

!!! danger
    Please use this carefully:
    - Make sure to only commit the changes to the doctests originating from
      your changes to the code.
    - The doctests also serve as actual tests, so make absolutely sure that the
      output is still mathematically correct.

!!! note
    This must be called at the top level, it does not work inside a function.
"""
macro doctest_fix(args...)
  pos, kws = _split_macro_args(args)
  length(pos) <= 1 || throw(ArgumentError("@doctest_fix takes at most one argument"))
  filter!(kw -> kw.args[1] !== :fix, kws)
  push!(kws, Expr(:kw, :fix, true))
  return _docs_toplevel(__source__, :doctest, pos, kws)
end

# The functions are kept alongside the macros because they are what everyone
# has in their fingers, and because downstream code such as OscarDevTools calls
# them. `OscarDocs` may only become loadable during the call, hence
# `invokelatest`.
const _oscardocs_uuid = UUID("39b4aa4e-cfe7-45c1-b576-a6789dd0b603")

function _oscardocs()
  docs_env()
  return Base.require(Base.PkgId(_oscardocs_uuid, "OscarDocs"))
end

@doc raw"""
    build_doc(; doctest=false, warnonly=true, open_browser=true, start_server=false)

Build the OSCAR manual. Equivalent to [`Oscar.@build_doc`](@ref), which should
be preferred at the REPL.
"""
build_doc(; kwargs...) = Base.invokelatest(_oscardocs().build, Oscar; kwargs...)

@doc raw"""
    doctest(; fix::Bool = false)
    doctest(f::Function; fix::Bool = false, set_meta::Bool = false)
    doctest(path::String; fix::Bool = false, set_meta::Bool = false)

Run the doctests of OSCAR, of a single function, or of every file whose path
contains `path`. Equivalent to [`Oscar.@doctest`](@ref), which should be
preferred at the REPL.
"""
doctest(; kwargs...) = Base.invokelatest(_oscardocs().doctest, Oscar; kwargs...)
doctest(what::Union{Function,String}; kwargs...) =
  Base.invokelatest(_oscardocs().doctest, Oscar, what; kwargs...)

@doc raw"""
    doctest_fix(f::Function; set_meta::Bool = false)
    doctest_fix(path::String; set_meta::Bool = false)

Run doctests as [`Oscar.doctest`](@ref) does, replacing their expected output
with what OSCAR produces. Please use carefully.
"""
doctest_fix(what::Union{Function,String}; kwargs...) = doctest(what; fix=true, kwargs...)

################################################################################
#
#  Experimental documentation
#
################################################################################

# Create symbolic links from any documentation directory in `experimental` into
# `docs/src`
function link_experimental_docs()
  # Remove symbolic links from earlier runs
  expdocdir = joinpath(oscardir, "docs", "src", "Experimental")
  for x in readdir(expdocdir; join=true)
    !islink(x) && continue
    pkg = splitpath(x)[end]
    if !(pkg in exppkgs)
      # We don't know this link, let's remove it
      rm(x)
    end
  end

  for pkg in exppkgs
    # Set symlink inside docs/src/experimental
    symlink_link = joinpath(oscardir, "docs", "src", "Experimental", pkg)
    symlink_target = joinpath(oscardir, "experimental", pkg, "docs", "src")

    if !ispath(symlink_target)
      continue
    end

    if !ispath(symlink_link)
      symlink(symlink_target, symlink_link)
    elseif !islink(symlink_link) || readlink(symlink_link) != symlink_target
      error("""$symlink_link already exists, but is not a symlink to $symlink_target
      Please investigate the contents of $symlink_link,
      optionally move them somewhere else and delete the directory once you are done.""")
    end
  end

  return nothing
end
