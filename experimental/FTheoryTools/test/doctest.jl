using Pkg
using Documenter
using FTheoryTools

DocMeta.setdocmeta!(FTheoryTools, :DocTestSetup, :(using FTheoryTools, Random); recursive = true)

doctest(FTheoryTools)
