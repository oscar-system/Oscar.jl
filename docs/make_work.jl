#
# This file is included by docs/make.jl and by a helper function
# in src/utils/docs.jl
#
module BuildDoc

using Documenter, DocumenterCitations, JSON

include("documenter_helpers.jl")
include("citation_style.jl")

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

function doit(
  Oscar::Module;
  warnonly=false,
  local_build::Bool=false,
  doctest::Union{Bool,Symbol}=true,
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
  for pkg in sort(Oscar.exppkgs)
    pkg == "ExperimentalTemplate" && continue
    pkgdocs = setup_experimental_package(Oscar, pkg)
    if length(pkgdocs) > 0
      append!(collected, pkgdocs)
    end
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

  function get_rev(uuid::Base.UUID)
    deps = Documenter.Pkg.dependencies()
    @assert haskey(deps, uuid)
    if !isnothing(deps[uuid].git_revision)
      return deps[uuid].git_revision
    else
      return "v$(deps[uuid].version)"
    end
  end
  aarev = get_rev(Base.PkgId(Oscar.AbstractAlgebra).uuid)
  nemorev = get_rev(Base.PkgId(Oscar.Nemo).uuid)
  heckerev = get_rev(Base.PkgId(Oscar.Hecke).uuid)
  singularrev = get_rev(Base.PkgId(Oscar.Singular).uuid)

  cd(joinpath(oscardir, "docs")) do
    DocMeta.setdocmeta!(Oscar, :DocTestSetup, Oscar.doctestsetup(); recursive=true)
    DocMeta.setdocmeta!(Oscar.Hecke, :DocTestSetup, :(using Hecke); recursive=true)
    DocMeta.setdocmeta!(Oscar.AbstractAlgebra, :DocTestSetup, :(using AbstractAlgebra); recursive=true)
    DocMeta.setdocmeta!(Oscar.Nemo, :DocTestSetup, :(using Nemo); recursive=true)

    if doctest !== false
      Documenter.doctest(Oscar; fix = doctest === :fix, doctestfilters=Oscar.doctestfilters())
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
      warnonly=warnonly,
      treat_markdown_warnings_as_error=!warnonly,
      checkdocs=:none,
      pages=doc,
      remotes=Dict(
        Base.pkgdir(Oscar.AbstractAlgebra) => (Remotes.GitHub("Nemocas", "AbstractAlgebra.jl"), aarev),
        Base.pkgdir(Oscar.Nemo) => (Remotes.GitHub("Nemocas", "Nemo.jl"), nemorev),
        Base.pkgdir(Oscar.Hecke) => (Remotes.GitHub("thofma", "Hecke.jl"), heckerev),
        Base.pkgdir(Oscar.Singular) => (Remotes.GitHub("oscar-system", "Singular.jl"), singularrev),
      ),
      plugins=[bib],
    )
  end

  # remove the copied documentation again
  for pkg in other_packages
    dstbase = normpath(oscardir, "docs", "src", string(nameof(pkg)))
    rm(dstbase; recursive=true, force=true)
  end
  
  # postprocessing, for the search index
  docspath = normpath(joinpath(oscardir, "docs"))
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

end # module BuildDoc

using .BuildDoc
