# TotalImage

This experimental OSCAR module computes the actual image of a rational map as
an alternating tree of closed sets. It is a Julia/OSCAR translation of the
[Macaulay2 TotalImage package](https://github.com/coreysharris/TotalImage),
version 0.2, by Corey Harris, Mateusz Michałek, and Emre Sertöz.

The implementation uses OSCAR ideals and ring homomorphisms throughout. Its
main entry points are `partial_image`, `total_image`, `is_contained`,
and `is_closed_image`.

Load the module with `using Oscar.TotalImage`.

The API is experimental and may change before promotion into OSCAR's main
source tree.
