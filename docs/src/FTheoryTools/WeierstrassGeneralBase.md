```@meta
CurrentModule = FTheoryTools
```

```@contents
Pages = ["Weierstrass.md"]
```

# Global Weierstrass models

This method constructs a global Weierstrass model over a base space that is not
fully specified. Rather, it assums that a base space exists such that
the Weierstrass sections ``f`` and ``g`` are well-defined, so that the
Weierstrass model in question is well-defined.

For many practical applications, one wishes to assume a further specialize
the Weierstrass sections ``f`` and ``g``. This has the advantage that one can
engineer singularity loci or even the singularity type over a specific locus.
To some extend, this is the backbone of many F-theory constructions.

Consequently, the construction of such models accepts a polynomial ring whose
variables are the sections used in the desired factorization of the Weierstrass
sections ``f`` and ``g``.

In theory, one can consider the indeterminates of this auxiliary ring as a
(possilby redundant) set of local coordinate of the base space. For the computer
implementation, this ring will therefore serve as the coordinate ring of an auxiliary
toric base space, namely an affine space with those coordinates.

For such geometries, we support the following functionality.

## Constructors

We support the following constructors:
```@docs
GlobalWeierstrassModel(f::fmpq_mpoly, g::fmpq_mpoly, auxiliary_base_ring::MPolyRing)
```

## Attributes

```@docs
weierstrass_section_f(w::WeierstrassModelOverGeneralBaseSpace)
weierstrass_section_g(w::WeierstrassModelOverGeneralBaseSpace)
weierstrass_polynomial(w::WeierstrassModelOverGeneralBaseSpace)
auxiliary_base_space(w::WeierstrassModelOverGeneralBaseSpace)
auxiliary_ambient_space(w::WeierstrassModelOverGeneralBaseSpace)
```
