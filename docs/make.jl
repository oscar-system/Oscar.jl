push!(LOAD_PATH,"../src/")
using Documenter, FTheoryTools, DocumenterMarkdown

# ensure FTheoryTools is loaded for the doc tests
DocMeta.setdocmeta!(FTheoryTools, :DocTestSetup, :(using FTheoryTools, Random); recursive = true)

makedocs(sitename="FTheoryTools")

deploydocs(repo = "github.com/HereAround/FTheoryTools.jl.git")
