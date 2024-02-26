# Quadratic forms and isometries

This project is a complement to the code about *hermitian lattices* available
in Hecke. We aim here to connect Hecke and GAP to handle some algorithmic
methods regarding quadratic forms with their isometries. In particular,
the integration of this code within Oscar is necessary to benefit from all the
performance of GAP with respect to computations with groups and automorphisms in
general.

For now, the project covers methods regarding rational and integral quadratic
forms.

## Content

We introduce two new structures
* `QuadSpaceWithIsom`
* `ZZLatWithIsom`

The former parametrizes pairs $(V, f)$ where $V$ is a rational quadratic form
and $f$ is an isometry of $V$. The latter parametrizes pairs $(L, f)$ where
$L$ is an integral quadratic form, also known as $\mathbb Z$-lattice and $f$
is an isometry of $L$. One of the main features of this project is the
enumeration of isomorphism classes of pairs $(L, f)$, where $f$ is an isometry
of finite order with at most two prime divisors. The methods we resort to
for this purpose are developed in the paper [BH23](@cite).

We also provide some algorithms computing isomorphism classes of primitive
embeddings of integral lattices following Nikulin's theory. More precisely, the
function `primitive_embeddings` offers, under certain conditions,
the possibility to compute representatives of primitive embeddings and classify
them in different ways. Note nonetheless that these functions are not efficient
in the case were the discriminant groups have a large number of subgroups.

## Status

Currently, the project features the following:

- enumeration of conjugacy classes of isometries of finite order for even lattices (in the case of at most 2 prime divisors);
- enumeration of conjugacy classes of isometries with irreducible and reciprocal minimal polynomial for integral lattices (with maximal equation order);
- primitive embeddings/extensions for integral lattices;
- equivariant primitive extensions for integral lattices;
- miscellaneous operations on integral/rational quadratic form endowed with an isometry.

## Current applications of this project

The project was initiated by S. Brandhorst and T. Hofmann for classifying
finite subgroups of automorphisms of K3 surfaces. Our current goal is to use
this code, and further extensions of it, to classify finite subgroups of
bimeromorphic self-maps of *hyperkaehler manifolds*, which are a higher
dimensional analogues of K3 surface.

## Notice to the user

Since this project is still under development, feel free to try any feature and
report all the bugs you may have found. Any suggestions for improvements or
extensions are more than welcome. Refer to the next section to know who you
should contact and how. Do not hesitate either to ask for new features - we
will be glad to add anything you may need for your research.

One may expect many things to vary within the next months: name of the
functions, available features, performance. This is due to the fact that the
current version of the code is still at an experimental stage.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en),
* [Tommy Hofmann](https://www.thofma.com/)
* [Stevell Muller](https://www.math.uni-sb.de/ag/brandhorst/index.php?option=com_content&view=article&id=26:muller-en-1&catid=18&lang=en&Itemid=114).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on GitHub](https://github.com/oscar-system/Oscar.jl).
