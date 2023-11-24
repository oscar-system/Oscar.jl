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
doctestsetup() = :(using Oscar; Oscar.AbstractAlgebra.set_current_module(@__MODULE__))

# use tempdir by default to ensure a clean manifest (and avoid modifying the project)
function doc_init(;path=mktempdir())
  global docsproject = path
  if !isfile(joinpath(docsproject,"Project.toml"))
    cp(joinpath(oscardir, "docs", "Project.toml"), joinpath(docsproject,"Project.toml"))
  end
  Pkg.activate(docsproject) do
    # we dev all "our" packages with the paths from where they are currently
    # loaded
    for pkg in [AbstractAlgebra, Nemo, Hecke, Singular, GAP, Polymake]
      Pkg.develop(path=Base.pkgdir(pkg))
    end
    Pkg.develop(path=oscardir)
    Pkg.instantiate()
    Base.include(Main, joinpath(oscardir, "docs", "make_work.jl"))
  end
end

function get_document(set_meta::Bool)
  if !isdefined(Main, :Documenter)
    error("you need to do `using Documenter` first")
  end

  if isdefined(Main.Documenter, :Document)
    Document = Main.Documenter.Document
  else
    Document = Main.Documenter.Documents.Document
  end
  doc = Document(root = joinpath(oscardir, "docs"), doctest = :fix)

  if Main.Documenter.DocMeta.getdocmeta(Oscar, :DocTestSetup) === nothing || set_meta
    Main.Documenter.DocMeta.setdocmeta!(Oscar, :DocTestSetup, Oscar.doctestsetup(); recursive=true)
  end

  if isdefined(Main.Documenter, :DocTests)
    doctest = Main.Documenter.DocTests.doctest
  else
    doctest = Main.Documenter._doctest
  end

  return doc, doctest
end

"""
    doctest_fix(f::Function; set_meta::Bool = false)

Fixes all doctests for the given function `f`.
"""
function doctest_fix(f::Function; set_meta::Bool = false)
  S = Symbol(f)
  doc, doctest = get_document(set_meta)

  #essentially inspired by Documenter/src/DocTests.jl
  pm = parentmodule(f)
  bm = Base.Docs.meta(pm)
  md = bm[Base.Docs.Binding(pm, S)]
  for s in md.order
    doctest(md.docs[s], Oscar, doc)
  end
end

"""
    doctest_fix(n::String; set_meta::Bool = false)

Fixes all doctests for the file `n`, ie. all files in Oscar where
`n` occurs in the full pathname of.
"""
function doctest_fix(n::String; set_meta::Bool=false)
  doc, doctest = get_document(set_meta)

  walkmodules(Oscar) do m
    #essentially inspired by Documenter/src/DocTests.jl
    bm = Base.Docs.meta(m)
    for (k, md) in bm
      for s in md.order
        if occursin(n, md.docs[s].data[:path])
          doctest(md.docs[s], Oscar, doc)
        end
      end
    end
  end
end

# copied from JuliaTesting/Aqua.jl
function walkmodules(f, x::Module)
  f(x)
  for n in names(x; all=true)
    # `isdefined` and `getproperty` can trigger deprecation warnings
    if Base.isbindingresolved(x, n) && !Base.isdeprecated(x, n)
      isdefined(x, n) || continue
      y = getproperty(x, n)
      if y isa Module && y !== x && parentmodule(y) === x
        walkmodules(f, y)
      end
    end
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

function start_doc_preview_server(;open_browser::Bool = true, port::Int = 8000)
  build_dir = normpath(Oscar.oscardir, "docs", "build")
  cmd = """
        using Pkg;
        Pkg.activate(temp = true);
        Pkg.add("LiveServer");
        using LiveServer;
        LiveServer.serve(dir = "$build_dir", launch_browser = $open_browser, port = $port);
        """
  live_server_process = run(`$(Base.julia_cmd()) -e $cmd`, wait = false)
  atexit(_ -> kill(live_server_process))
  @info "Starting server with PID $(getpid(live_server_process)) listening on 127.0.0.1:$port"
  return nothing
end

@doc raw"""
    build_doc(; doctest=false, warnonly=true, open_browser=true, start_server=false)

Build the manual of `Oscar.jl` locally and open the front page in a
browser.

The optional parameter `doctest` can take three values:
  - `false`: Do not run the doctests (default).
  - `true`: Run the doctests and report errors.
  - `:fix`: Run the doctests and replace the output in the manual with
    the output produced by Oscar. Please use this option carefully.

In GitHub Actions the Julia version used for building the manual is 1.9 and
doctests are run with >= 1.7. Using a different Julia version may produce
errors in some parts of Oscar, so please be careful, especially when setting
`doctest=:fix`.

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
`build_doc(start_server = true)` for the first build and
`build_doc(open_browser = false)` for following builds and only refresh the
browser tab.

When working on the manual the `Revise` package can significantly speed
up running `build_doc`. First, install `Revise` in the following way:
```
using Pkg ; Pkg.add("Revise")
```
Second, restart Julia and load `Revise` before Oscar:
```
using Revise, Oscar;
```
The first run of `build_doc` will take the usual few minutes, subsequent runs
will be significantly faster.
"""
function build_doc(; doctest::Union{Symbol, Bool} = false, warnonly = true, open_browser::Bool = true, start_server::Bool = false)
  versioncheck = (VERSION.major == 1) && (VERSION.minor >= 7)
  versionwarn = """The Julia reference version for the doctests is 1.7 or later, but you are using
                $(VERSION). Running the doctests will produce errors that you do not expect."""
  if doctest != false && !versioncheck
    @warn versionwarn
  end
  if !isdefined(Main, :BuildDoc)
    doc_init()
  end
  Pkg.activate(docsproject) do
    Base.invokelatest(Main.BuildDoc.doit, Oscar; warnonly=warnonly, local_build=true, doctest=doctest)
  end
  if start_server
    start_doc_preview_server(open_browser = open_browser)
  elseif open_browser
    open_doc()
  end
  if doctest != false && !versioncheck
    @warn versionwarn
  end
end
