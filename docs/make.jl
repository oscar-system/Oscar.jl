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
            "Tools" => ["FTheoryTools/WeierstrassGeneralBase.md"],
            "Tools" => ["FTheoryTools/Tate.md"],
            "Tools" => ["FTheoryTools/TateGeneralBase.md"],
            ]
    )

deploydocs(repo = "github.com/Julia-meets-String-Theory/FTheoryTools.jl.git")
