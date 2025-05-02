# Symmetric intersections

The main purpose of this code is to systematize the computation of
ideals defining intersections of hypersurfaces of the same degree in a given
projective space preserved under a faithful linear action of a finite group.

Note that the code does not require to know explicitly the action of the group
on the projective space. However, one requires to be able to construct the group
in Oscar.

## Warning

The project is experimental. The names of the functions, the functionalities
available, the extent of the content as well as the existence of this project
may vary in the coming months.

Moreover, the current implementation of representation theory of finite groups
available throughout *should not be taken as a standard* for Oscar: it is a
temporary solution for the purpose of this project. This part of the code was
not written by an expert on representation theory, and some of the expensive
functions might not be computationally optimized. This is also why the content
is limited to "what was needed". Nevertheless, it was made to be usable by any
user.

## Content

Here is a list of what kind of functionalities are available:

### Linear representations of finite groups

Basics functions for linear representations of finite groups over fields of
characteristic zero, as well as some complements to character theory, are
available. Namely, one should be able to create and manipulate linear
representations of a given finite group over a big enough field *of
characteristic zero*.

For instance, one can compute the direct sum of representations, the dual
representation, the symmetric/exterior powers of a representation, but also
bases for the *isotypical components*. Finally, few things have been added for
characters, namely the possibility of switching between characters and
representations, and enumerating constituents of a given degree.

### Projective representations

These representations are implemented to allow the user to classify, up to
similarity, faithful projective representations of a finite group over a
field of characteristic zero. Such representations represent faithful
linear actions of groups on projective spaces.

The main function from this part is called `faithful_projective_representations`,
which returns representatives of similarity classes of faithful projective
representations of a given group, of a given dimension.

*Almost* all the functionalities available for linear representations are
available for projective representations.

### Symmetric Grassmannians

We call here *symmetric Grassmannians* spaces which parametrize submodules of a
given group algebra modules. In some sense, we add a group action to the usual
Grassmannian variety, and look for invariant subspaces of a fixed vector space
on which a group acts.

Even though some functionalities about them are available, the naming
"symmetric Grassmannians" is not standard and the usage is quite limited. This
is why the content related to them is not directly exported to the user:
one may use "Oscar.bla(arguments, of, bla)" to call the corresponding function
"bla".

Given a group algebra module $M$, represented by a linear representation of a
finite group, one can create the space parametrising submodules $M$ of a given
dimension, with given character or with given determinant character. Once
created, it is possible to ask for a *standard element*, some parametrization
data, a defining ideal of the corresponding variety (these spaces of modules are
projective varieties) as well as the dimension as a projective variety.

Here by standard element we mean the submodule obtained by setting
all the parameters to 1.

### Parameter space for symmetric intersections

The ultimate goal of this project is to compute defining ideals for
intersections of hypersurfaces preserved under the action of a group on a
projective space. Up to now, only the case of *faithful* actions and intersections
of hypersurfaces of the *same degree* are available.

The main function is `symmetric_intersections`.

## Future of the project

Any idea of improvement or complement is welcome. The possible new directions
for this project will be to extend to non-faithful projective actions, as well
as treating the case of intersections of hypersurfaces of different degrees. A bit
more ambitiously would be to use tools from representation theory in order to
compute sections of vector bundles on other homogeneous varieties whose zero
locus is preserved under the action of a given group.

As of now, the principal use of this code is to compute defining ideals of some
K3 surfaces.

## Contact

Please direct questions about this part of Oscar to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?option=com_content&view=article&id=6&Itemid=107&lang=en)
* [Stevell Muller](https://www.math.uni-sb.de/ag/brandhorst/index.php?option=com_content&view=article&id=26:muller-en-1&catid=18&lang=en&Itemid=114)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on GitHub](https://github.com/oscar-system/Oscar.jl).
