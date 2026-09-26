# Subalgebra (SAGBI) Bases

contributor/author: Kaleb Ruscitti (kaleb (dot) ruscitti (at) uwaterloo.ca)

## Aims

This package aims to expose a Julia interface for working with SAGBI 
(Subalgebra Analogue to Grobner Bases for Ideals) bases. The
['SubalgebraBases'](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/SubalgebraBases/html/index.html)
package in Macaulay2 is an example of a fully-featured package for working 
with SAGBI bases.

The goal is implement three main features:
1. Checking if a generating set of a subalgebra is a SAGBI basis.
2. Subducting a polynomial by a (possibly partial) SAGBI basis.
3. Completing a generating set for a subalgebra to a SAGBI basis. 

Though there are other SAGBI related features that could be useful, I think these three features form a useful baseline functionality.

## Status

At time of writing all three main features are implemented. The main computational bottlenecks call to 4ti2.

More detailed testing and documentation needs to be written, and the code hasn't been reviewed in detail by anyone other than the author.