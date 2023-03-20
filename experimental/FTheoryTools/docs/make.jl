push!(LOAD_PATH,"../src/")
using Documenter, DocumenterCitations, FTheoryTools, DocumenterMarkdown, Oscar

# ensure FTheoryTools is loaded for the doc tests
DocMeta.setdocmeta!(FTheoryTools, :DocTestSetup, :(using FTheoryTools, Random); recursive = true)

# read the bibliography
bib = CitationBibliography(joinpath(FTheoryTools.ftheorytoolsdir, "docs", "references.bib"), sorting = :nyt)

# build documentation
makedocs(
    bib,
    format = Documenter.HTML(),
    sitename = "FTheoryTools",
    modules = [FTheoryTools],
    clean = true,
    doctest = false,
    pages = ["index.md",
             "weierstrass.md",
             "tate.md",
             "developer.md",
             "References" => "references.md",
            ]
    )

# deploy the documentation
deploydocs(repo = "github.com/Julia-meets-String-Theory/FTheoryTools.jl.git")
