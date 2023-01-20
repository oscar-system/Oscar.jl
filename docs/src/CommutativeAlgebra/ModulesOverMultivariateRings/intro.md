```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```

# Introduction

Our focus in this section is on finitely presented modules over rings from the following list:
- multivariate polynomial rings (OSCAR type `MPolyRing`),
- quotients of multivariate polynomial rings  (OSCAR type `MPolyQuo`), and
- localizations of the above rings (OSCAR types `MPolyLocalizedRing`, `MPolyQuoLocalizedRing`).
Hence, if not mentioned otherwise, the word *module* refers to a finitely presented module over a
ring of one of the above types. Offering a sparse way of implementing free modules, the
OSCAR type `FreeMod` provides the basis for implementing all modules discussed here. More concretely,
the general way of implementing a module is to represent it as a *subquotient*, that is, as a
submodule of a quotient of a free module. Note that subquotients form the smallest class of
modules which naturally  includes both submodules and quotients of free modules.

!!! note
    Most functions in this section rely on Gröbner (standard) bases techniques. Thus, the functions
    should not be applied to modules over rings other than those from the list above. See the Linear
	Algebra chapter for module types which are designed to handle modules over Euclidean
	domains.

!!! note
    For simplicity of the presentation in what follows, functions are often only illustrated by examples  
    with focus on modules over multivariate polynomial rings, but work similarly for modules over
	a ring of any of the above types.
