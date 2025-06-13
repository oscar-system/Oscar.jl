# Complex Reflection Groups

## Aim 
The aim is to implement [complex reflection groups](https://en.wikipedia.org/wiki/Complex_reflection_group) along with data that is known or can be computed.

By [Ulrich Thiel](https://ulthiel.com/math), 2023.

## Example

```
julia> W = ComplexReflectionGroupType([ (4,) , (10,5,3)] )
Complex reflection group type G4 x G(10,5,3)

julia> order(W)
28800
```

## Status

- [X] Structure ```ComplexReflectionGroupType``` for type (equivalence up to conjugacy) of complex reflection group following Shephard–Todd notation.

- [X] Basic data for complex reflection group types: ```order```, ```is_imprimitive```, ```num_reflections```, ```is_well_generated```, ```is_pseudo_real```, ```is_real```, ```is_rational```, ```coxeter_number```, ```is_spetsial```, ```num_reflections```, ```num_hyperplanes```, ```num_reflection_classes``` 

- [X] Explicit models via ```complex_reflection_group``` (Magma, Lehrer–Taylor, and CHEVIE)

- [ ] "Reflection library" (list of reflections, conjugacy classes of reflections, hyperplane orbits etc.) 

- [ ] Conjugacy classes and character tables (main problem will be to have compatibility with CHEVIE; it may be easiest/best to [compute](https://webusers.imj-prg.fr/~jean.michel/gap3/htm/chap087.htm) the labeling)

- [ ] Explicit models of the irreducible representations (same problem as above + problem of how to store the data)

- [ ] Because it somehow fits in here as well: symplectic reflection groups

- [ ] With the same explanation as above: Drinfeld–Hecke algebras (which includes symplectic reflection algebras and rational Cherednik algebras)

 - [ ] Identification of the type of a complex reflection group (challenging)
 