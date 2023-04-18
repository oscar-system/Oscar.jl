#
# This file is included by docs/make.jl and by a helper function
# in src/Oscar.jl
#
module BuildDoc

using Documenter, DocumenterCitations


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

# We monkey-patch Base.walkdir to use true as default value for follow_symlinks
# (normally false is the default), in order to "trick" the Documenter code into
# following those symlinks.
# See also:
# https://github.com/JuliaDocs/Documenter.jl/pull/552
# https://github.com/JuliaLang/julia/blob/master/doc/make.jl#L19
Base.walkdir(str::String) = Base.walkdir(str; follow_symlinks=true)


# When we read a `doc.main` from an experimental package, we need to equip all
# its entries with a prefix to fit with our docs. The doc.main of an
# experimental package will contain paths relative to
# `experimental/PACKAGE_NAME/docs/src`. When generating the docs a symlink is
# set in `docs/src/Experimental/PACKAGE_NAME` pointing to
# `experimental/PACKAGE_NAME/docs/src`. Hence the paths in `doc.main` need to
# get the prefix `Experimental/PACKAGE_NAME`.
#
# Example:
# 1. cat experimental/PlaneCurve/docs/doc.main:
# [
#    "plane_curves.md",
# ]
# after `add_prefix_to_experimental_docs` becomes
# [
#    "Experimental/PlaneCurve/plane_curves.md",
# ]
#
# 2. cat experimental/FTheoryTools/docs/doc.main 
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
add_prefix_to_experimental_docs(Oscar::Module, docs::Pair{String, Vector{T}}, prefix::String) where T = Pair{String, Vector{T}}(docs[1], add_prefix_to_experimental_docs(Oscar, docs[2], prefix))
add_prefix_to_experimental_docs(Oscar::Module, docs::Vector{T}, prefix::String) where T = T[add_prefix_to_experimental_docs(Oscar, entry, prefix) for entry in docs]


function setup_experimental_package(Oscar::Module, package_name::String)
  doc_main_path = joinpath(Oscar.oscardir, "experimental", package_name, "docs/doc.main")
  if !isfile(doc_main_path)
    return[]
  end

  # Set symlink inside docs/src/experimental
  symlink_link = joinpath(Oscar.oscardir, "docs/src/Experimental", package_name)
  symlink_target = joinpath(Oscar.oscardir, "experimental", package_name, "docs", "src")

  if !ispath(symlink_target)
    return []
  end

  if !ispath(symlink_link)
    symlink(symlink_target, symlink_link)
  elseif !islink(symlink_link) || readlink(symlink_link) != symlink_target
      error("$symlink_link already exists, but is not a symlink to $symlink_target
Please investigate the contents of $symlink_link,
optionally move them somewhere else and delete the directory once you are done.")
  end

  # Read doc.main of package
  exp_s = read(doc_main_path, String)
  exp_doc = eval(Meta.parse(exp_s))

  # Prepend path
  prefix = "Experimental/" * package_name * "/"
  result = add_prefix_to_experimental_docs(Oscar, exp_doc, prefix)
  return result
end

function doit(Oscar::Module; strict::Bool = true, local_build::Bool = false, doctest::Union{Bool,Symbol} = true)

  # Remove symbolic links from earlier runs
  expdocdir = joinpath(Oscar.oscardir, "docs", "src", "Experimental")
  for x in readdir(expdocdir; join=true)
    islink(x) && rm(x)
  end

  # include the list of pages, performing substitutions
  s = read(joinpath(Oscar.oscardir, "docs", "doc.main"), String)
  doc = eval(Meta.parse(s))
  collected = Any["Experimental/intro.md"]
  for pkg in Oscar.exppkgs
    pkgdocs = setup_experimental_package(Oscar, pkg)
    if length(pkgdocs) > 0
      append!(collected, pkgdocs)
    end
  end
  push!(doc, ("Experimental" => collected))

  # Load the bibliography
  bib = CitationBibliography(joinpath(Oscar.oscardir, "docs", "oscar_references.bib"), sorting = :nyt)

  # Copy documentation from Hecke, Nemo, AnstratAlgebra
  other_packages = [
    (Oscar.Hecke,Oscar.heckedir),
    (Oscar.Nemo, Oscar.nemodir),
    (Oscar.AbstractAlgebra, Oscar.aadir),
    ]
  for (pkg, pkgdir) in other_packages
      srcbase = normpath(pkgdir, "docs", "src")
      dstbase = normpath(Oscar.oscardir, "docs", "src", string(nameof(pkg)))

      # clear the destination directory first
      rm(dstbase, recursive=true, force=true)

      for (root, dirs, files) in walkdir(srcbase)
          for dir in dirs
              d = normpath(joinpath(dstbase, relpath(root, srcbase), dir))
              mkpath(d)
          end
          for file in files
              # HACK: delete Hecke's bibliography, to avoid warnings of the
              # form "Warning: 'Eis95' is not unique" which actually turn into
              # errors down the road
              if file == "references.md"
                continue
              end
              src = normpath(joinpath(root, file))
              dst = normpath(joinpath(dstbase, relpath(root, srcbase), file))
              cp(src, dst; force = true)
              chmod(dst, 0o644)
          end
      end
  end

  cd(joinpath(Oscar.oscardir, "docs")) do

    DocMeta.setdocmeta!(Oscar, :DocTestSetup, :(using Oscar); recursive = true)
    DocMeta.setdocmeta!(Oscar.Hecke, :DocTestSetup, :(using Hecke); recursive = true)
    DocMeta.setdocmeta!(Oscar.AbstractAlgebra, :DocTestSetup, :(using AbstractAlgebra); recursive = true)
    DocMeta.setdocmeta!(Oscar.Nemo, :DocTestSetup, :(using Nemo); recursive = true)


    makedocs(bib,
           format   = Documenter.HTML(prettyurls = !local_build, collapselevel = 1),
           sitename = "Oscar.jl",
           modules = [Oscar, Oscar.Hecke, Oscar.Nemo, Oscar.AbstractAlgebra, Oscar.Singular],
           clean = true,
           doctest = doctest,
           strict = strict,
           checkdocs = :none,
           pages    = doc)
  end

  # remove the copied documentation again
  for (pkg, pkgdir) in other_packages
      dstbase = normpath(Oscar.oscardir, "docs", "src", string(nameof(pkg)))
      rm(dstbase, recursive=true, force=true)
  end

end

end # module BuildDoc

using .BuildDoc
