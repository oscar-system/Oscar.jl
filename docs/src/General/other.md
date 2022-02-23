# Notes for users of other computer algebra systems

## General differences

- TODO: ``2^{100}`` is evaluated to 0; write ``fmpz(2)^100`` to get
  a long

## Notes for GAP users

- TODO
- TODO: also talk about how to use it from Oscar?

## Notes for Singular users

- TODO
- TODO: also talk about how to use it from Oscar?

## Notes for Polymake users

- OSCAR (and Julia) is `1`-based, meaning that it counts from `1`, rather than
  from `0` like polymake. For most properties we have taken care of the
  translation but be aware that it might pop up at some point and generate
  confusion.
- Polyhedra and polyhedral complexes in OSCAR are represented inhomogeneously,
  i.e. without the leading `1` for vertices or `0` for rays. Hence constructors
  take points, rays, and lineality generators separately.
- TODO: also talk about how to use it from Oscar?

## Notes for Magma users

- TODO

## Notes for SageMath users

- TODO
