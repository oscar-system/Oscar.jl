```@meta
CurrentModule = FTheoryTools
```

```@contents
Pages = ["Weierstrass.md"]
```

# Global Weierstrass models

## ... over concrete bases

A global Weierstrass model describes a particular form of an elliptic fibration.
We focus on elliptic fibrations over base 3-folds ``B3``. Consider
the weighted projective space ``\mathbb{P}^{2,3,1}`` with coordinates
``x, y, z``. In addition, consider
- ``f \in H^0( B_3, \overline{K}_{B_3}^{\otimes 4} )``,
- ``g \in H^0( B_3, \overline{K}_{B_3}^{\otimes 6} )``,
Then form a ``\mathbb{P}^{2,3,1}``-bundle over ``B3`` such that
- ``x`` transforms as a section of ``2 \overline{K}_{B_3}``,
- ``y`` transforms as a section of ``3 \overline{K}_{B_3}``,
- ``z`` transforms as a section of ``0 \overline{K}_{B_3} = \mathcal{O}_{B_3}``.
In this 5-fold ambient space, a global Weierstrass model is the hypersurface defined
by the vanishing of the Weierstrass polynomial ``P_W = x^3 - y^2 + f x z^4 + g z^6``.

To construct a Weierstrass model as a hypersurface in an ambient space,
we first need to construct the ambient space in question. For a toric base,
one way to achieve this is to first focus on the Cox ring of the toric
ambient space. This ring must be graded such that the Weierstrass polynomial
is homogeneous and cuts out a Calabi-Yau hypersurface. Given this grading,
one can perform a triangulation. This triangulation will typically
take a long time to complete and yield a large number of candidate ambient
spaces. Typically, one wishes to focus on those spaces which contain the
toric base space in a manifest way. But even this criterion usually
allows for many ambient spaces to be used.

To circumvent this obstacle, all our constructors operate in the opposite direction.
That is, they begins by extracting the rays and maximal cones of the chosen
toric base space. Subsequently, those rays and cones are extended
to form one of the many toric ambient spaces. Note that the so-constructed
ambient space need not be smooth.

We support the following constructors:
```@docs
GlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)
GlobalWeierstrassModelOverProjectiveSpace()
GlobalWeierstrassModel(f::MPolyElem_dec{fmpq, fmpq_mpoly}, g::MPolyElem_dec{fmpq, fmpq_mpoly}, base::Oscar.AbstractNormalToricVariety)
```

## ... over not fully specified bases

A global Weierstrass model can also be constructed over a base space that
is not fully specified. Rather, it assumes that a base space exists such that
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
```@docs
GlobalWeierstrassModel(f::fmpq_mpoly, g::fmpq_mpoly, auxiliary_base_ring::MPolyRing)
```


## Attributes

```@docs
weierstrass_section_f(w::GlobalWeierstrassModel)
weierstrass_section_g(w::GlobalWeierstrassModel)
weierstrass_polynomial(w::GlobalWeierstrassModel)
cy_hypersurface(w::GlobalWeierstrassModel)
Oscar.:discriminant(w::GlobalWeierstrassModel)
```
In case the Weierstrass model is constructed over a not fully specified base,
it is nontheless possibly to construct an auxiliary base space as well as an
auxiliary ambient space. The (auxiliary) base and ambient space can
be accessed with the following functions:
```@docs
toric_base_space(w::GlobalWeierstrassModel)
toric_ambient_space(t::GlobalWeierstrassModel)
```
To tell if they are auxiliary or not, one can use `base_fully_specified(w)`,
which returns `true` in case the Tate model is defined over a concrete base and
`false` otherwise.


## Properties

```@docs
base_fully_specified(t::GlobalWeierstrassModel)
```
