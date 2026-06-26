```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id modules_multivariate)

OSCAR provides two implementations of finitely presented modules, with
different mathematical scopes and algorithmic approaches.

* **Commutative algebra modules** are discussed in this section. They work
  over multivariate polynomial rings, quotients of such rings and their
  localizations, and employ Gröbner basis methods.

* **Linear algebra modules** are documented in the
  [Linear Algebra modules](@ref modules_linear_algebra) section. They
  originate from AbstractAlgebra.jl, work over Euclidean domains and fields,
  and employ algorithms from linear algebra, such as the Smith Normal Form.


## Scope and functionality of commutative algebra modules

The implementation described in this section is limited to finitely presented
modules over the following rings:

* multivariate polynomial rings (OSCAR type `MPolyRing`),
* quotients of multivariate polynomial rings (OSCAR type `MPolyQuoRing`), and
* localizations of the above rings (OSCAR types `MPolyLocRing`, `MPolyQuoLocRing`).

The OSCAR type `FreeMod` provides a sparse implementation of free modules.
More generally, modules are represented as *subquotients*, that is, as
submodules of quotients of free modules. This representation naturally
includes both submodules and quotients of free modules.

!!! note
    Most functions in this section rely on Gröbner (standard) basis techniques and
    should therefore only be applied to modules over the rings listed above.

!!! note
    For simplicity of the presentation, many examples use modules over
    multivariate polynomial rings. The same functionality is generally available
    for modules over the other supported ring types.
