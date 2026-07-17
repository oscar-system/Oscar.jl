# Origami

## Aims

The aim of this project is to make the functions from our GAP package `Origami` available in OSCAR.
In the long run, the goal is also to re-implement more and more functions natively in Julia.

## Status

Major todos are the following:

* Extend tests and documentation.
* Make the Veech group and other function return values proper OSCAR objects. Use the newly implemented `experimental/ModularGroup` package for this.
* Extend and overwork the `origamis` function so that it lists H(1, 1, 1, 1) correctly and also eventually supports other strata.
* Implement database functionality.
