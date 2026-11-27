# Orbits on Subspaces

Let ``G`` be a group of ``d \times d`` matrices over a prime field ``F``.
The functions described here compute representatives, lengths and stabilizers
of the orbits of ``G`` on the ``k``-dimensional subspaces of ``F^d``.

Each subspace is identified with its basis matrix in reduced row echelon
form, and that matrix with an integer between ``1`` and the number of
``k``-dimensional subspaces, as in [OBr90](@cite).
The orbit algorithm runs on these integers.
Where the solvable radical of ``G`` is known, its orbits are computed along
a polycyclic generating sequence, see [LNS84](@cite); this takes about one
image per subspace instead of one per generator of ``G``.

## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what this means.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Max Horn](https://math.rptu.de/en/wgs/agag/people/head/prof-dr-max-horn)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

## Functions

```@docs
subspace_orbit_representatives_and_lengths
subspace_orbit_representatives_and_stabilizers
```
