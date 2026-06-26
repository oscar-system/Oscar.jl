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

General textbooks offering details on theory and algorithms include:
- [Lan71](@cite)
- [Kir16](@cite)


## Matrix implementations

OSCAR does not rely on Julia's built-in matrix types for two independent
reasons:

* In empty matrices (with zero rows or columns), only the type of the matrix
  entries is known. For the algebraic types used throughout OSCAR, this
  information is generally insufficient to construct new elements, so
  operations such as `zero(T)` cannot always be implemented correctly.

* Julia's linear algebra is designed around types that embed into the real or
  complex numbers. For example, `det(ones(Int, (1,1))) == 1.0`, so the fact
  that the result is exactly the integer `1` is lost. More generally, many of
  the rings used in OSCAR cannot be embedded into the real or complex numbers
  at all.

Instead, OSCAR makes use of several matrix implementations. Most exact dense
matrix types, together with many associated operations, are provided by Nemo.jl,
which in turn builds on the generic matrix interfaces of AbstractAlgebra.jl. In
addition, OSCAR defines matrix types for specialized mathematical structures,
and wraps matrix types provided
by external libraries.

Consequently, OSCAR uses different matrix implementations depending on the
underlying mathematical structures and performance requirements. While these
implementations share many common operations, they do not all support exactly
the same functionality. Where appropriate, conversion between OSCAR matrix
types and Julia's native `Matrix` type is supported.


## Topics covered

The main topics covered in this section are:
- [Matrix functionality](@ref matrix_functionality_chapter)
- [Matrix spaces](@ref "Matrix Spaces")
- [Linear solving](@ref solving_chapter)
- [Eigenvalues and -spaces](@ref "Eigenvalues and -spaces")
- [Generic matrix algebras](@ref "Generic matrix algebras")
- [Sparse linear algebra](@ref)
- [Modules](@ref "Finitely presented modules")


## For developers

For developers interested in the implementation details of matrices, see
[Matrix implementation](@ref "Matrix implementation").


## Tutorials

We encourage you to take a look at the tutorials on linear algebra in
OSCAR, which can be found [here](https://www.oscar-system.org/tutorials/LinearAlgebra/).


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Claus Fieker](https://math.rptu.de/en/wgs/agag/people/head/fieker),
* [Tommy Hofmann](https://www.thofma.com/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
