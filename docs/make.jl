push!(LOAD_PATH,"../src/")
using Documenter, FTheoryTools, DocumenterMarkdown, Oscar

# ensure FTheoryTools is loaded for the doc tests
DocMeta.setdocmeta!(FTheoryTools, :DocTestSetup, :(using FTheoryTools, Random; using FTheoryTools.Oscar); recursive = true)

makedocs(
    format = Documenter.HTML(),
    sitename = "FTheoryTools",
    modules = [FTheoryTools],
    clean = true,
    doctest = false,
    pages = ["index.md",
            "FTheoryTools/Introduction.md",
            "Tools" => ["FTheoryTools/Weierstrass.md"],
            ]
    )

deploydocs(repo = "github.com/HereAround/FTheoryTools.jl.git")
