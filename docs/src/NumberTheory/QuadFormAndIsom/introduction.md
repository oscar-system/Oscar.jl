```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

This project is a complement to the code about *hermitian lattices* available
in Hecke. We aim here to connect Hecke and GAP to handle some algorithmic
methods regarding quadratic forms with their isometries. In particular,
the integration of this code within OSCAR is necessary to benefit from all the
performance of GAP with respect to computations with groups and automorphisms
in general.

For now, the project covers methods regarding rational and integral quadratic
forms.

## Content

We introduce two new structures
* [`QuadSpaceWithIsom`](@ref)
* [`ZZLatWithIsom`](@ref)

The former parametrizes pairs $(V, f)$ where $V$ is a rational quadratic form
and $f$ is an isometry of $V$. The latter parametrizes pairs $(L, f)$ where
$L$ is an integer quadratic form, also known as $\mathbb Z$-lattice and $f$
is an isometry of $L$. One of the main features of this project is the
enumeration of isomorphism classes of pairs $(L, f)$, where $f$ is an isometry
of finite order. The methods we resort to for this purpose are developed
in the paper [BH23](@cite).

We also provide some algorithms computing isomorphism classes of primitive
embeddings of integral lattices following Nikulin's theory. More precisely, the
function [`primitive_embeddings`](@ref) offers, under certain conditions,
the possibility to compute representatives of primitive embeddings and classify
them in different ways. Note nonetheless that these functions are not efficient
in the case were the discriminant groups have a large number of subgroups.

## Status

Currently, the project features the following:

* enumeration of conjugacy classes of isometries of finite order for even
  lattices;
* enumeration of conjugacy classes of isometries with irreducible and
  reciprocal minimal polynomial for integral lattices (with maximal equation
  order);
* primitive embeddings/extensions for integral lattices;
* equivariant primitive extensions for integral lattices;
* miscellaneous operations on integral/rational quadratic form endowed with
  an isometry.

## Current applications of this project

The project was initiated by Brandhorst and Hofmann [BH23](@cite) for
classifying finite subgroups of automorphisms of K3 surfaces. Our current goal
is to use this code, and further extensions of it, to classify finite subgroups
of bimeromorphic self-maps of *hyperkähler manifolds*, which are higher
dimensional analogues of K3 surface.

## Tutorials

No tutorials available at the moment.

## Examples

No examples available at the moment.

## Notice to the user

### Report an issue

If you are working with some objects of type `QuadSpaceWithIsom` or
`ZZLatWithIsom` and you need to report an issue, you can produce directly some
lines of codes helping to reconstruct your example. This can help the reviewers
to understand your issue and assist you. We have implemented a method
`to_oscar` which prints few lines of codes for reconstructing your example.

```jldoctest
julia> V = quadratic_space(QQ, 2);

julia> Vf = quadratic_space_with_isometry(V, neg = true)
Quadratic space of dimension 2
  with isometry of finite order 2
  given by
  [-1    0]
  [ 0   -1]

julia> Oscar.to_oscar(Vf)
G = matrix(QQ, 2, 2, [1 0; 0 1]);
V = quadratic_space(QQ, G);
f = matrix(QQ, 2, 2, [-1 0; 0 -1]);
Vf = quadratic_space_with_isometry(V, f);

julia> Lf = lattice(Vf)
Integer lattice of rank 2 and degree 2
  with isometry of finite order 2
  given by
  [-1    0]
  [ 0   -1]

julia> Oscar.to_oscar(Lf)
B = matrix(QQ, 2, 2, [1 0; 0 1]);
G = matrix(QQ, 2, 2, [1 0; 0 1]);
L = integer_lattice(B; gram=G);
f = matrix(QQ, 2, 2, [-1 0; 0 -1]);
Lf = integer_lattice_with_isometry(L, f);
```

## Make the code more talkative

Within the code, there are more hidden messages and testing which are disabled
by default. If you plan to experiment with the codes with your favourite
examples, you may want to be able to detect some issues to be reported, as well
as knowing what the code is doing. Indeed, some functions might take time in
term of compilation but also computations. For this, you can enable these extra
tests and printings by setting:
                                                                                                     
```julia
Oscar.set_lwi_level(2)
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en),
* [Tommy Hofmann](https://www.thofma.com/),
* [Stevell Muller](https://sites.google.com/view/stevellmuller/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on GitHub](https://github.com/oscar-system/Oscar.jl).
