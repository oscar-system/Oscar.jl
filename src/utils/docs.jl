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
    #ugly: needs to be in sync with the docs/make_docs.jl file
    Main.Documenter.DocMeta.setdocmeta!(Oscar, :DocTestSetup, :(using Oscar; Oscar.AbstractAlgebra.set_current_module(@__MODULE__)); recursive = true)
  end

  return doc
end

"""
    doctest_fix(f::Function; set_meta::Bool = false)

Fixes all doctests for the given function `f`.
"""
function doctest_fix(f::Function; set_meta::Bool = false)
  S = Symbol(f)
  doc = get_document(set_meta)

  #essentially inspired by Documenter/src/DocTests.jl
  bm = Main.Documenter.DocSystem.getmeta(Oscar)
  md = bm[Base.Docs.Binding(Oscar, S)]
  for s in md.order
    Main.Documenter.DocTests.doctest(md.docs[s], Oscar, doc)
  end
end

"""
    doctest_fix(n::String; set_meta::Bool = false)

Fixes all doctests for the file `n`, ie. all files in Oscar where
  `n` occurs in the full pathname of.
"""
function doctest_fix(n::String; set_meta::Bool = false)
  doc = get_document(set_meta)

  #essentially inspired by Documenter/src/DocTests.jl
  bm = Main.Documenter.DocSystem.getmeta(Oscar)
  for (k, md) = bm
    for s in md.order
      if occursin(n, md.docs[s].data[:path])
        Main.Documenter.DocTests.doctest(md.docs[s], Oscar, doc)
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


@doc raw"""
    build_doc(; doctest=false, strict=false, open_browser=true)

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
