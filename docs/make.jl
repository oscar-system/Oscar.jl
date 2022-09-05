push!(LOAD_PATH,"../src/")
using Documenter, FTheoryTools, DocumenterMarkdown

# ensure FTheoryTools is loaded for the doc tests
DocMeta.setdocmeta!(FTheoryTools, :DocTestSetup, :(using FTheoryTools, Random); recursive = true)

makedocs(
    format = Documenter.HTML(),
    sitename = "FTheoryTools",
    modules = [FTheoryTools],
    clean = true,
    doctest = false,
    pages = ["index.md",
            "FTheoryTools/Introduction.md",
            "Tools" => ["FTheoryTools/Functions.md"],
            ]
    )

deploydocs(repo = "github.com/HereAround/FTheoryTools.jl.git")
