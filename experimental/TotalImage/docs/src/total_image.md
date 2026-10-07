```@meta
CurrentModule = Oscar.TotalImage
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Total images of rational maps

This experimental module computes the image of a rational map as a
constructible set. It translates the algorithms from the
[Macaulay2 TotalImage package](https://github.com/coreysharris/TotalImage) by
Corey Harris, Mateusz Michałek, and Emre Sertöz.

For a rational map from a projective variety to projective space, the module
first computes the Zariski closure of the image. It then resolves the base
locus through the graph closure, projects the exceptional divisor, and repeats
the calculation on the relevant pullbacks. The result is an alternating tree:
the root is included, its children are excluded, their children are included,
and so on.

## Basic use

The Cremona transformation illustrates the representation.

```julia
using Oscar.TotalImage

R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
tree = total_image([y*z, x*z, x*y])
```

Tree vertices are ideals in an automatically generated target ring with
coordinates `b[0], b[1], ...`. Use `vertices(tree)` and `edges(tree)` to inspect
the data directly. Edge indices follow Julia's one-based convention.
Use `is_contained(x, tree)` to test whether a projective or affine point lies
in the represented constructible set.

Nonhomogeneous coordinates, or a nonhomogeneous domain ideal, are interpreted
as affine input and homogenized automatically. The resulting tree is pruned at
infinity and dehomogenized before it is returned.

Set the `:TotalImage` verbosity level to `1` to print progress information.

## Public interface

```@docs
ConstructibleSetTree
is_contained
partial_image
total_image
is_closed_image
```

## Status

This code is experimental. In particular, the tree API and the choice of
linear sections used for maps with positive-dimensional generic fibers may
change.
