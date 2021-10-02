using Documenter, JToric

DocMeta.setdocmeta!(JToric, :DocTestSetup, :(using JToric, Random; using JToric.Oscar); recursive = true)

doctest(JToric)
