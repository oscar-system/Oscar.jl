```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id linear_algebra)

The linear algebra part of OSCAR provides functionality for handling matrices,
linear solving, modules, lattices and related structures.

The documentation in this section combines functionality provided directly by
OSCAR together with functionality originating from several of its underlying
packages, including AbstractAlgebra.jl, Nemo.jl and Hecke.jl.

The main topics covered in this section are:
- [Matrix functionality](@ref matrix_functionality_chapter)
- [Matrix spaces](@ref "Matrix Spaces")
- [Linear solving](@ref solving_chapter)
- [Eigenvalues and -spaces](@ref "Eigenvalues and -spaces")
- [Generic matrix algebras](@ref "Generic matrix algebras")
- [Sparse linear algebra](@ref)
- [Modules](@ref "Finitely presented modules")

For developers interested in the implementation details of matrices, see
[Matrix implementation](@ref "Matrix implementation").

General textbooks offering details on theory and algorithms include:
- [Lan71](@cite)
- [Kir16](@cite)

## Tutorials

We encourage you to take a look at the tutorials on linear algebra in
OSCAR, which can be found [here](https://www.oscar-system.org/tutorials/LinearAlgebra/).


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Claus Fieker](https://math.rptu.de/en/wgs/agag/people/head/fieker),
* [Tommy Hofmann](https://www.thofma.com/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
