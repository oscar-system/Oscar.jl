using Pkg
Pkg.add("Oscar")
using Oscar
using Documenter
using FTheoryTools

DocMeta.setdocmeta!(FTheoryTools, :DocTestSetup, :(using FTheoryTools, Random; using Oscar); recursive = true)

doctest(FTheoryTools)
