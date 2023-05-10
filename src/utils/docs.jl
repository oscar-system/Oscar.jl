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

