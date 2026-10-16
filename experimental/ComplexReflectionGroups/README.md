# Complex Reflection Groups

## Aim

The aim is to implement [complex reflection groups](https://en.wikipedia.org/wiki/Complex_reflection_group), their types as complex linear groups, type-level invariants, and several exact matrix realizations.

By [Ulrich Thiel](https://ulthiel.com/math), 2023.

## Example

```julia
julia> T = complex_reflection_group_type([4, (10, 5, 3)]; fixed_space_dimension=2)
Complex reflection group type G(10,5,3) x G4 with 2-dimensional fixed space

julia> (order(T), rank(T), fixed_space_dimension(T), degree(T))
(28800, 5, 2, 7)

julia> W = complex_reflection_group(T)
Matrix group of degree 7
  over cyclotomic field of order 30

julia> complex_reflection_group_type(W) == T
true

julia> length(complex_reflections(W))
41
```

## Status

- [x] `ComplexReflectionGroupType` for full ambient conjugacy types. A type records the dimension of the fixed space together with the normalized unordered collection of its nontrivial irreducible Shephard--Todd types.

- [x] Type-level invariants including `order`, `degree`, `rank`, `fixed_space_dimension`, `is_imprimitive`, `number_of_reflections`, `number_of_hyperplanes`, `number_of_reflection_classes`, `degrees`, and `codegrees`.

- [x] Explicit models via `complex_reflection_group` (Magma, Lehrer--Taylor, and CHEVIE), including trivial ambient summands and exact block-diagonal products over compatible coefficient fields.

- [x] Type-level and concrete reflecting-hyperplane orbits, conjugacy classes of reflections, and a hyperplane-first `reflection_library`.

- [ ] Conjugacy classes of arbitrary elements and character tables (main problem will be compatibility with CHEVIE; it may be easiest to [compute](https://webusers.imj-prg.fr/~jean.michel/gap3/htm/chap087.htm) the labeling).

- [ ] Explicit models of the irreducible representations (same problem as above + problem of how to store the data)

- [x] Cotangent doubling of complex reflection groups, with the source, block decomposition, and alternating form retained as realization data.

- [ ] General symplectic reflection groups and their reflection data.

- [ ] With the same explanation as above: Drinfeld–Hecke algebras (which includes symplectic reflection algebras and rational Cherednik algebras)

- [ ] Recognition of the type of an unmarked matrix group (deliberately separate from the type-aware constructors).
