# Tableaux

A **Young diagram** is a diagram of finitely many empty "boxes" arranged
in left-justified rows, with the row lengths in non-increasing order. The
box in row $i$ and and column $j$ has the **coordinates** $(i, j)$. Listing
the number of boxes in each row gives a partition $\lambda$ of a non-negative
integer $n$ (the total number of boxes of the diagram). The diagram is
then said to be of **shape** $\lambda$. Conversely, one can associate to any
partition $\lambda$ a Young diagram in the obvious way, so Young diagrams are
just another way to look at partitions.

A **Young tableau** of shape $\lambda$ is a filling of the boxes of the Young
diagram of $\lambda$ with elements from some set. After relabeling we can (and
will) assume that we fill from a set of integers from $1$ up to some number,
which in applications is often equal to $n$.

In OSCAR, a tableau is internally stored as an array of arrays and is
represented by the type `YoungTableau{T}` which is a subtype of
`AbstractVector{AbstractVector{T}}`, where `T` is the integer type of the
filling. As for partitions, one may increase performance by casting into
smaller integer types, e.g. `Int8`.

```@docs
young_tableau
```

## Operations

```@docs
hook_length
hook_lengths
shape
weight
reading_word
```

## Semistandard tableaux

```@docs
is_semistandard
semistandard_tableaux
```

## Standard tableaux

```@docs
is_standard
standard_tableaux
```
```@docs
number_of_standard_tableaux
```
The number $f^\lambda$ of standard Young tableaux of shape $\lambda$ is computed
using the *hook length formula*

$$f^\lambda = \frac{n!}{\prod_{i, j} h_\lambda(i, j)},$$

where the product is taken over all boxes in the Young diagram of $\lambda$ and
$h_\lambda$ denotes the hook length of the box $(i, j)$.

```@docs
schensted
bump!
```
