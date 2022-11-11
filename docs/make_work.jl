#
# This file is included by docs/make.jl and by a helper function
# in src/Oscar.jl
#
module BuildDoc

using Documenter, DocumenterMarkdown, DocumenterCitations


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

function doit(Oscar::Module; strict::Bool = true, local_build::Bool = false, doctest::Union{Bool,Symbol} = true)

  # include the list of pages, performing substitutions
  s = read(joinpath(Oscar.oscardir, "docs", "doc.main"), String)
  doc = eval(Meta.parse(s))

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
    DocMeta.setdocmeta!(Oscar.Graphs, :DocTestSetup, :(using Oscar; using Oscar.Graphs); recursive = true)
    DocMeta.setdocmeta!(Oscar.AbstractAlgebra, :DocTestSetup, :(using AbstractAlgebra); recursive = true)
    DocMeta.setdocmeta!(Oscar.Nemo, :DocTestSetup, :(using Nemo); recursive = true)


    makedocs(bib,
           format   = Documenter.HTML(prettyurls = !local_build, collapselevel = 1),
  #         format   = Documenter.HTML(),
  #         format   = Markdown(),
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
