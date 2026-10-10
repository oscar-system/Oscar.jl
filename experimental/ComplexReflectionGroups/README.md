# Complex Reflection Groups

## Aim

The aim is to implement [complex reflection groups](https://en.wikipedia.org/wiki/Complex_reflection_group), their type-level invariants, and several exact matrix realizations.

By [Ulrich Thiel](https://ulthiel.com/math), 2023.

## Example

```julia
julia> W = ComplexReflectionGroupType([4, (10, 5, 3)])
Complex reflection group type G(10,5,3) x G4

julia> order(W)
28800
```

## Status

- [x] `ComplexReflectionGroupType` for conjugacy types of essential reflection groups, represented in normalized Shephard–Todd notation.

- [x] Type-level invariants including `order`, `rank`, `is_imprimitive`, `number_of_reflections`, `number_of_hyperplanes`, `number_of_reflection_classes`, `degrees`, and `codegrees`.

- [x] Explicit models via `complex_reflection_group` (Magma, Lehrer–Taylor, and CHEVIE), including exact block-diagonal products over compatible coefficient fields.

- [x] Type-level and concrete reflecting-hyperplane orbits, conjugacy classes of reflections, and a hyperplane-first `reflection_library`.

- [ ] Conjugacy classes of arbitrary elements and character tables (main problem will be compatibility with CHEVIE; it may be easiest to [compute](https://webusers.imj-prg.fr/~jean.michel/gap3/htm/chap087.htm) the labeling).

- [ ] Explicit models of the irreducible representations (same problem as above + problem of how to store the data)

- [x] Initial symplectic-reflection-group support via cotangent doubling, with the source, block decomposition, and alternating form retained as realization data.

- [ ] With the same explanation as above: Drinfeld–Hecke algebras (which includes symplectic reflection algebras and rational Cherednik algebras)

- [ ] Recognition of the type of an unmarked matrix group (deliberately separate from the type-aware constructors).
