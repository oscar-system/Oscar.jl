using Documenter, JToric

DocMeta.setdocmeta!(JToric, :DocTestSetup, :(using JToric, Random); recursive = true)

doctest(JToric)
