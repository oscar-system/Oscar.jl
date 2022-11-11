```@meta
CurrentModule = FTheoryTools
```

```@contents
Pages = ["Tate.md"]
```

# Global Tate models without specified base space

This method constructs a global Tate model over a base space that is not
fully specified. Rather, it assums that a base space exists such that
the Tate sections ``a_i`` are well-defined so that the Tate model in
question is well-defined.

For many practical applications, one wishes to assume a further factorization
of the Tate sections ``a_i``. This has the advantage that one can engineer
singularity loci or even the singularity type over a specific locus. This is
the backbone of many F-theory constructions.

To this end, this method accepts a polynomial ring whose variables are the sections
used in the desired factorization of the Tate sections ``a_i``. For example, if we
desired a factorization:
- ``a_1 = a_{10} w^0``,
- ``a_2 = a_{21} w^1``,
- ``a_3 = a_{32} w^2``,
- ``a_4 = a_{43} w^3``,
- ``a_6 = a_{65} w^5``,
then the polynomial ring in question is the ring with indeterminates
``a_{10}``, ``a_{21}``, ``a_{32}``, ``a_{43}``, ``a_{65}`` and ``w``.

In theory, one can consider these indeterminates as local coordinate of the base space.
For the computer implementation, this ring will therefore serve as the coordinate
ring of an auxiliary toric base space, namely an affine space with those coordinates.

For such geometries, we support the following functionality.


## Constructors

```@docs
GlobalTateModel(ais::Vector{fmpq_mpoly}, auxiliary_base_ring::MPolyRing)
```


## Attributes

```@docs
tate_section_a1(t::TateModelOverGeneralBaseSpace)
tate_section_a2(t::TateModelOverGeneralBaseSpace)
tate_section_a3(t::TateModelOverGeneralBaseSpace)
tate_section_a4(t::TateModelOverGeneralBaseSpace)
tate_section_a6(t::TateModelOverGeneralBaseSpace)
tate_polynomial(t::TateModelOverGeneralBaseSpace)
auxiliary_base_space(t::TateModelOverGeneralBaseSpace)
auxiliary_ambient_space(t::TateModelOverGeneralBaseSpace)
```
