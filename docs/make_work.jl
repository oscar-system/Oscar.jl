using Documenter, Oscar, DocumenterMarkdown, DocumenterCitations

# Copy documentation from Hecke, Nemo, AnstratAlgebra
for pkg in [Oscar.Hecke, Oscar.Nemo, Oscar.AbstractAlgebra]
    build = normpath(pkgdir(Oscar), "docs", "src", string(nameof(pkg)))
    source = normpath(pkgdir(pkg), "docs", "src")
    for (root, dirs, files) in walkdir(source)
        for dir in dirs
            d = normpath(joinpath(build, relpath(root, source), dir))
            mkpath(d)
        end
        for file in files
            src = normpath(joinpath(root, file))
            dst = normpath(joinpath(build, relpath(root, source), file))
            cp(src, dst; force = true)
        end
    end
end


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


const hecke = "Hecke"
const nemo = "Nemo"
const aa = "AbstractAlgebra"

bib = CitationBibliography(joinpath(@__DIR__, "oscar_references.bib"), sorting = :nyt)

function doit(strict::Bool = true, local_build::Bool = false)

  s = read(joinpath(Oscar.oscardir, "docs", "doc.main"), String)
  doc = eval(Meta.parse(s))

  cd(joinpath(Oscar.oscardir, "docs")) do

  DocMeta.setdocmeta!(Oscar, :DocTestSetup, :(using Oscar); recursive = true)
  DocMeta.setdocmeta!(Hecke, :DocTestSetup, :(using Hecke); recursive = true)

  makedocs(bib,
         format   = Documenter.HTML(prettyurls = !local_build, collapselevel = 1),
#         format   = Documenter.HTML(),
#         format   = Markdown(),
         sitename = "Oscar.jl",
         modules = [Oscar, Hecke, Nemo, AbstractAlgebra, Singular],
         clean = true,
         doctest = true,
         strict = strict,
         checkdocs = :none,
         pages    = doc)
end

end
