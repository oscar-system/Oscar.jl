```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id modules_linear_algebra)

OSCAR provides two implementations of finitely presented modules, with
different mathematical scopes and algorithmic foundations.

* **Linear algebra modules** are discussed in this section. They originate
  from AbstractAlgebra.jl, work over Euclidean domains and fields, and
  employ algorithms from linear algebra, such as the Smith Normal Form.

* **Commutative algebra modules** over multivariate polynomial rings,
  quotients of such rings and their localizations are discussed in the
  [Commutative Algebra](@ref modules_multivariate) chapter. They employ
  Gröbner basis methods and are designed for multivariate polynomial rings,
  as well as quotients and localizations thereof.


## Scope and functionality of linear algebra modules

The implementation described in this section is limited to finitely presented
modules over fields and Euclidean domains.

Free modules and vector spaces are available over fields and Euclidean domains,
respectively. Submodule, quotient module and direct sum constructions can then
be applied recursively to these.

Invariant decompositions can be computed using the Smith Normal Form. The
system also provides module homomorphisms and isomorphisms, building on top of
the map interface.
