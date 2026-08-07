# Origami

## Aims

The aim of this project is to make the functions from our GAP package `Origami` available in OSCAR.
In the long run, the goal is also to re-implement more and more functions natively in Julia.

## Acknowledgments

The `possible_lengths_and_heights` function, used to generate all origamis in a given stratum,
is a direct translation of `CylinderDiagram.widths_and_heights_iterator` from [Vincent Delecroix's
`surface-dynamics` package](https://flatsurf.github.io/surface-dynamics/).

## Status

Major todos are the following:

* Extend tests and documentation.
* Make the Veech group and other function return values proper OSCAR objects. Use the newly implemented `experimental/ModularGroup` package for this.
* Extend and overwork the `origamis` function so that it lists H(1, 1, 1, 1) correctly and also eventually supports other strata.
* Implement database functionality.
