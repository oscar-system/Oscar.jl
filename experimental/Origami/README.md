# Origami

## Aims

The aim of this project is to make the functions from our GAP package `Origami` available in OSCAR.
In the long run, the goal is also to re-implement more and more functions natively in Julia.

## Status

Major todos are the following:
* As soon as the GAP packages `Origami` and `ModularGroup` have new clean releases, include them as artifacts.
* Extend tests and documentation.
* Make the Veech group and other function return values proper OSCAR objects. It might be necessary to create an `experimental/ModularGroup` package for this.

