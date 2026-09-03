#
# Driver for building the OSCAR manual and for running its doctests.
#
# This package deliberately does *not* depend on Oscar: Oscar is always passed
# in as a `Module` argument. That keeps the whole OSCAR dependency graph out of
# the environment providing this package, so that environment need not be
# re-resolved whenever OSCAR's dependencies change. See `Oscar.docs_env`.
#
module OscarDocs

using Documenter, DocumenterCitations

include("citation_style.jl")
include("documenter_helpers.jl")
include("doctest.jl")
include("preview.jl")

# Overwrite printing to make the header not full of redundant nonsense
# Turns
#   Hecke.Order - Method
# into
#   Order - Method

# To remove the '-'
# Documenter.Utilities.print_signature(io::IO, signature)        = print(io, signature)

# To remove the "Method", "Type", "Module" use the following
# Documenter.Utilities.doccat(b::Base.Docs.Binding, ::Type)  = ""
# doccat(::Type)     = ""
# doccat(::Module)   = ""

# Remove the module prefix
Base.print(io::IO, b::Base.Docs.Binding) = print(io, b.var)

# When we read a `doc.main` from an experimental package, we need to equip all
# its entries with a prefix to fit with our docs. The doc.main of an
# experimental package will contain paths relative to
# `experimental/PACKAGE_NAME/docs/src`. When generating the docs a symlink is
# set in `docs/src/Experimental/PACKAGE_NAME` pointing to
# `experimental/PACKAGE_NAME/docs/src`. Hence the paths in `doc.main` need to
# get the prefix `Experimental/PACKAGE_NAME`.
#
# Example:
# cat experimental/FTheoryTools/docs/doc.main 
# [
#    "F-Theory Tools" => [
#       "introduction.md",
#       "weierstrass.md",
#       "tate.md",
#    ],
# ]
# after `add_prefix_to_experimental_docs` becomes
# [
#    "F-Theory Tools" => [
#       "Experimental/FTheoryTools/introduction.md",
#       "Experimental/FTheoryTools/weierstrass.md",
#       "Experimental/FTheoryTools/tate.md",
#    ],
# ]
#
# Since the entries of a `doc.main` vary in type, we have split this up into
# three functions.
add_prefix_to_experimental_docs(Oscar::Module, docs::String, prefix::String) = joinpath(prefix, docs)
add_prefix_to_experimental_docs(Oscar::Module, docs::Pair{String,String}, prefix::String) = Pair{String,String}(docs.first, add_prefix_to_experimental_docs(Oscar, docs.second, prefix))
add_prefix_to_experimental_docs(Oscar::Module, docs::Pair{String, Vector{T}}, prefix::String) where T = Pair{String, Vector{T}}(docs.first, add_prefix_to_experimental_docs(Oscar, docs.second, prefix))
add_prefix_to_experimental_docs(Oscar::Module, docs::Vector{T}, prefix::String) where T = T[add_prefix_to_experimental_docs(Oscar, entry, prefix) for entry in docs]


function setup_experimental_package(Oscar::Module, package_name::String)
  oscardir = Base.pkgdir(Oscar)
  doc_main_path = joinpath(oscardir, "experimental", package_name, "docs", "doc.main")
  if !isfile(doc_main_path)
    return []
  end

  # Assumes that a symbolic link from `experimental/package_name/docs/src`
  # to `docs/src/Experimental/package_name` has been created (or there is no
  # documentation for this package)
  if !ispath(joinpath(oscardir, "docs", "src", "Experimental", package_name))
    return []
  end

  # Read doc.main of package
  exp_s = read(doc_main_path, String)
  exp_doc = try
    eval(Meta.parse(exp_s))
  catch
    println("error while parsing $doc_main_path:")
    rethrow()
  end

  # Prepend path
  prefix = "Experimental/" * package_name * "/"
  result = add_prefix_to_experimental_docs(Oscar, exp_doc, prefix)
  return result
end

# The revision a package's sources correspond to, used for the "source" links in
# the manual. `Pkg.dependencies()` would be the obvious tool, but it lists the
# dependencies of the active project only, and OSCAR itself is not among those
# when the manual is built from OSCAR's own project.
function package_revision(pkg::Module)
  dir = Base.pkgdir(pkg)
  # `.git` is a file, not a directory, in worktrees and submodules
  if dir !== nothing && ispath(joinpath(dir, ".git"))
    try
      return readchomp(`git -C $dir rev-parse HEAD`)
    catch
      # not a usable checkout, fall through to the version number
    end
  end
  return "v$(pkgversion(pkg))"
end

@doc """
    build(Oscar::Module; doctest=false, warnonly=true, local_build=true,
                         open_browser=true, start_server=false)

Build the OSCAR manual. See `Oscar.@build_doc` for a description of the
keyword arguments.
"""
function build(
  Oscar::Module;
  doctest::Union{Bool,Symbol}=false,
  warnonly=true,
  local_build::Bool=true,
  open_browser::Bool=true,
  start_server::Bool=false,
)
  doctest === false || setup_doctest_stats(Oscar)
  with_docs_terminal(Oscar) do
    _build(Oscar; doctest, warnonly, local_build)
  end

  if start_server
    start_doc_preview_server(Oscar; open_browser)
  elseif open_browser
    open_doc(Oscar)
  end

  return nothing
end

# Always reached through `build`, which owns the defaults.
function _build(
  Oscar::Module;
  warnonly,
  local_build::Bool,
  doctest::Union{Bool,Symbol},
)
  oscardir = Base.pkgdir(Oscar)

  # include the list of pages, performing substitutions
  s = read(joinpath(oscardir, "docs", "doc.main"), String)
  doc = eval(Meta.parse(s))

  # Link experimental docs to `docs/src` and collect the documentation pages
  # Experimental documentation order:
  # 1. intro.md (how to add new projects)
  # 2. ExperimentalTemplate
  # 3. all other experimental packages (alphabetical)
  Oscar.link_experimental_docs()
  collected = Any["Experimental/intro.md"]
  append!(collected, setup_experimental_package(Oscar, "ExperimentalTemplate"))
  new_collected = Any[]
  for pkg in Oscar.exppkgs
    pkg == "ExperimentalTemplate" && continue
    pkgdocs = setup_experimental_package(Oscar, pkg)
    if length(pkgdocs) > 0
      push!(new_collected, pkgdocs)
    end
  end
  sort!(new_collected; by = pkgdocs -> lowercase(replace(first(pkgdocs[1]), r"\W" => "")))
  for pkgdocs in new_collected
    append!(collected, pkgdocs)
  end
  pos = findfirst(d -> d isa Pair && startswith(d[1], "Experimental"), doc)
  append!(doc[pos].second, collected)

  # Load the bibliography
  bib = CitationBibliography(
    joinpath(oscardir, "docs", "oscar_references.bib"); style=oscar_style
  )

  # Copy documentation from Hecke, Nemo, AbstractAlgebra
  other_packages = [Oscar.Hecke, Oscar.Nemo, Oscar.AbstractAlgebra]
  for pkg in other_packages
    srcbase = normpath(Base.pkgdir(pkg), "docs", "src")
    dstbase = normpath(oscardir, "docs", "src", string(nameof(pkg)))

    # clear the destination directory first
    rm(dstbase; recursive=true, force=true)

    for (root, dirs, files) in walkdir(srcbase)
      for dir in dirs
        d = normpath(dstbase, relpath(root, srcbase), dir)
        mkpath(d)
      end
      for file in files
        # HACK: delete Hecke's bibliography, to avoid warnings of the
        # form "Warning: 'Eis95' is not unique" which actually turn into
        # errors down the road
        if file == "references.md"
          continue
        end
        src = normpath(root, file)
        dst = normpath(dstbase, relpath(root, srcbase), file)
        if endswith(file, ".md")
          symlink(src, dst)
        else
          cp(src, dst; force=true)
        end
        chmod(dst, 0o644)
      end
    end
  end

  aarev = package_revision(Oscar.AbstractAlgebra)
  nemorev = package_revision(Oscar.Nemo)
  heckerev = package_revision(Oscar.Hecke)
  singularrev = package_revision(Oscar.Singular)
  oscarrev = package_revision(Oscar)

  cd(joinpath(oscardir, "docs")) do
    DocMeta.setdocmeta!(Oscar, :DocTestSetup, Oscar.doctestsetup(); recursive=true, warn=false)
    DocMeta.setdocmeta!(Oscar.Hecke, :DocTestSetup, :(using Hecke); recursive=true, warn=false)
    DocMeta.setdocmeta!(Oscar.AbstractAlgebra, :DocTestSetup, :(using AbstractAlgebra); recursive=true, warn=false)
    DocMeta.setdocmeta!(Oscar.Nemo, :DocTestSetup, :(using Nemo); recursive=true, warn=false)

    if doctest !== false
      Documenter.doctest(Oscar; fix = doctest === :fix, doctestfilters=Oscar.doctestfilters(), meta=Oscar.docmeta())
    end

    makedocs(;
      format=Documenter.HTML(;
        prettyurls=!local_build,
        collapselevel=1,
        size_threshold=409600,
        size_threshold_warn=204800,
        size_threshold_ignore=["manualindex.md"],
        canonical="https://docs.oscar-system.org/stable/",
      ),
      sitename="Oscar.jl",
      modules=[Oscar, Oscar.Hecke, Oscar.Nemo, Oscar.AbstractAlgebra, Oscar.Singular],
      clean=true,
      doctest=false,
      meta=Oscar.docmeta(),
      warnonly=warnonly,
      treat_markdown_warnings_as_error=!warnonly,
      checkdocs=:none,
      pages=doc,
      remotes=Dict(
        Base.pkgdir(Oscar.AbstractAlgebra) => (Remotes.GitHub("Nemocas", "AbstractAlgebra.jl"), aarev),
        Base.pkgdir(Oscar.Nemo) => (Remotes.GitHub("Nemocas", "Nemo.jl"), nemorev),
        Base.pkgdir(Oscar.Hecke) => (Remotes.GitHub("thofma", "Hecke.jl"), heckerev),
        Base.pkgdir(Oscar.Singular) => (Remotes.GitHub("oscar-system", "Singular.jl"), singularrev),
        Base.pkgdir(Oscar) => (Remotes.GitHub("oscar-system", "Oscar.jl"), oscarrev),
      ),
      plugins=[bib],
    )
  end

  # remove the copied documentation again
  for pkg in other_packages
    dstbase = normpath(oscardir, "docs", "src", string(nameof(pkg)))
    rm(dstbase; recursive=true, force=true)
  end

  patch_search_index(Oscar, doc, local_build)
end

# Documenter indexes every page it finds below `docs/src`, including the ones we
# only symlinked in. Restrict the index to the pages actually listed in
# `doc.main`.
function patch_search_index(Oscar::Module, doc, local_build::Bool)
  JSON = Oscar.JSON
  docspath = normpath(Base.pkgdir(Oscar), "docs")
  @info "Patching search index."
  # extract valid json from search_index.js
  run(pipeline(`sed -n '2p;3q' $(joinpath(docspath, "build", "search_index.js"))`, stdout=(joinpath(docspath, "build", "search_index.json")))) # imperfect file, but JSON parses it

  # extract paths from doc
  filelist=String[]
  docmain = doc
  while !isempty(docmain)
    n = pop!(docmain)
    if n isa Pair
      push!(docmain, last(n))
    elseif n isa String
      push!(filelist, n)
    elseif n isa Array{String}
      append!(filelist,n)
    elseif n isa Array
      append!(docmain,n)
    else
      error("err: $(typeof(n))")
    end
  end
  suffix = local_build ? ".html" : "/"
  filelist = replace.(filelist, r"\.md$"=>suffix)

  # read these files
  iosearchindex = open(joinpath(docspath, "build", "search_index.json"), "r")
  searchindex = JSON.parse(iosearchindex)
  close(iosearchindex)

  newsearchindex = []

  for item in searchindex
    if split(item["location"], "#")[1] in filelist
      push!(newsearchindex, item)
    end
  end


  # combine this to valid javascript again, and overwrite input
  ionewsearchindex = open(joinpath(docspath, "build", "search_index.js"), "w")
  write(ionewsearchindex, """var documenterSearchIndex = {"docs":\n""")
  JSON.print(ionewsearchindex, newsearchindex)
  write(ionewsearchindex, "\n}")
  close(ionewsearchindex)

  # clean up
  rm(joinpath(docspath, "build", "search_index.json"))
end

end # module OscarDocs
