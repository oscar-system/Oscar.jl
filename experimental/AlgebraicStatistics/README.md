# AlgebraicStatistics

## Aims

This experimental package provides tools for *algebraic statistics*,
facilitating the use of commutative algebra to solve problems in
statistics.

The interface is loosely modeled after the [`GraphicalModels`](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/GraphicalModels/html/index.html)
package of [`Macaulay2`](https://macaulay2.com) but it fully embraces
features of Julia's type system. We intend to first get the functionality
of this package on par with its `Macaulay2` predecessor before adding
new functionality.

## Status

- [X] Conditional independence models
  - [X] Gaussian ideals
  - [X] Discrete ideals
- [ ] Directed graphical models (Bayesian networks)
  - [ ] Markov properties
  - [ ] Trek separation
  - [X] Gaussian parametrization
  - [ ] Discrete parametrization
- [X] Undirected graphical models (Markov networks)
  - [ ] Markov properties
  - [X] Gaussian parametrization
  - [X] Discrete parametrization
- [ ] Phylogenetics
  - [X] Cavender--Farris--Neyman model
  - [X] Jukes--Cantor model
  - [X] Kimura2 / Kimura3
  - [X] Fourier parameters
- [ ] Maximum likelihood
