# Integer lattices with isometry

This project is a complement to the code about *hermitian lattices* available
on Hecke. We aim here to connect Hecke and GAP to handle some algorithmic
methods regarding integer lattices with their isometries. In particular,
the integration of this code to Oscar is necessary to benefit all the
performance of GAP with respect to computations with groups and automorphisms in
general.

## Content

We introduce the new type `ZZLatWithIsom` which parametrizes pairs $(L, f)$ where
$L$ is a non-denegerate $\mathbb{Z}$-lattice and $f$ is an isometry of $L$. One
of the main feature of this project is the enumeration of isomorphism classes of
pairs $(L, f)$ where $f$ is an isometry of finite order with at most two prime
divisors. The methods we resort to for this purpose are developed in the paper
[BH23].

We also provide some algorithms computing isomorphism classes of primitive
embeddings of even lattices following [Nikulin]. More precisely, the two
functions `primitive_embeddings_in_primary_lattice` and
`primitive_embeddings_of_primary_lattice` offer, under certain conditions,
the possibility of obtaining representatives of such classes of primitive
embeddings. Note nonetheless that these functions are not efficient in the case
were the discriminant groups are large.

## Status

This project has been slightly tested on simple and known examples. It is
currently being tested on a larger scale to test its reliability. Moreover,
there are still computational bottlenecks due to non optimized algorithms.

Among the possible improvements and extensions:
* Implement methods about for lattices with isometries of infinite order;
* Extend the methods for classification of primitive embeddings for the more
  general case (knowing that we lose efficiency for large discriminant groups);
* Implement methods for all kinds of equivariant primitive extensions (not
  necessarily admissibles);
* Import the methods for extending trivial discriminant actions on lattice whose
  discriminant group is an abelian $p$-group.

## Currently application of this project

The project was initiated by S. Brandhorst and T. Hofmann for classifying
finite subgroups of automorphisms of K3 surfaces. Our current goal is to use
this code, and further extension of it, to classify finite subgroups of
automorphisms and birational transformations on *hyperkaehler manifolds*, which
are a higher dimensional analog of K3 surface.

## Tutorials


## Notice to the user

Since this project is still under development, feel free to try any feature and
report all the bugs you may have found. Any suggestions for improvements or
extensions are more than welcome. Refer to the next section to know who you
should contact and how.

One may expect many things to vary within the next months: name of the functions,
available features, performance. This is due to the fact that the current
version of the code is still at an experimental stage.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en),
* [Tommy Hofmann](https://www.thofma.com/)
* [Stevell Muller](https://www.math.uni-sb.de/ag/brandhorst/index.php?option=com_content&view=article&id=26:muller-en-1&catid=18&lang=en&Itemid=114).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on GitHub](https://github.com/oscar-system/Oscar.jl).
