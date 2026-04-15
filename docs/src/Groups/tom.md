```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Tables of Marks

The concept of a *Table of Marks* was introduced by W. Burnside in his book
Theory  of Groups of Finite Order [Bur11](@cite).
Therefore a table of marks is sometimes called a *Burnside matrix*.
The table of marks of a finite group ``G`` is a matrix whose rows and columns
are labelled by the conjugacy classes of subgroups of ``G`` and where for
two subgroups ``H`` and ``K`` the ``(H, K)``-entry is the number of
fixed points of ``K`` in the transitive action of ``G`` on the right cosets
of ``H`` in ``G``.
So the table of marks characterizes the set of all permutation representations
of ``G``.
Moreover, the table of marks gives a compact description of the
subgroup lattice of ``G``, since from the numbers of fixed points
the numbers of conjugates of a subgroup ``K`` contained in a subgroup ``H``
can be derived.

For small groups the table of marks of ``G`` can be constructed directly
by first computing the entire subgroup lattice of ``G``,
see [`table_of_marks(G::Union{GAPGroup, FinGenAbGroup})`](@ref).
Besides that, the Table of Marks library [MNP24](@cite) provides access to
several hundred tables of marks of simple groups and maximal subgroups
of simple groups.
These tables of marks can be fetched via the names of these groups,
which coincide with the names of the character tables of these groups
in the Character Table Library, see [`table_of_marks(id::String)`](@ref).

Like the library of character tables, the library of tables of marks
can be used similar to group libraries (see [Group libraries](@ref))
in the sense that [`all_table_of_marks_names`](@ref) returns descriptions
of all those available tables of marks that have certain properties.

## Construct tables of marks

```@docs
table_of_marks(G::Union{GAPGroup, FinGenAbGroup})
table_of_marks(id::String)
Base.show(io::IO, ::MIME"text/plain", tom::GAPGroupTableOfMarks)
all_table_of_marks_names
is_table_of_marks_name
```

## Attributes and operations for tables of marks

```@docs
class_lengths(tom::GAPGroupTableOfMarks)
identifier(tom::GAPGroupTableOfMarks)
order(tom::GAPGroupTableOfMarks)
orders_class_representatives(tom::GAPGroupTableOfMarks)
representative(tom::GAPGroupTableOfMarks, i::Int)
```

## Marks vectors

The ``\mathbb Z``-linear combinations of the rows of the table of marks
of the group ``G`` can be interpreted as elements of
the integral *Burnside ring* of ``G``:
The rows of the table represent the isomorphism classes of
transitive ``G``-sets, the sum of two rows represents the isomorphism class
of the disjoint union of the two ``G``-sets,
and the pointwise product of two rows represents the isomorphism class
of the Cartesian product of the two ``G``-sets.
The coefficients of the decomposition of a linear combination of rows
can be computed by [`coordinates(chi::GAPGroupMarksVector)`](@ref).

The rows of a table of marks and their ``\mathbb Z``-linear combinations
are implemented as marks vector objects, with `parent` the table of marks.

length and iteration:

The length of a marks vector is the number of columns of the table of marks.
iteration is defined w.r.t. the ordering of columns.

arithmetic operations:

- `chi == psi`:
  two marks vectors are equal if and only if they belong to the same
  table of marks and have the same values,
- `chi + psi` and `chi - psi` are the pointwise sum and difference,
  respectively, of the two marks vectors `chi`, `psi`,
- `n * chi` is the pointwise `n`-fold sum of `chi`, for an integer `n`,
- `chi * psi` is the pointwise product of `chi` and `psi`,
- `zero(chi)` is the marks vector that is zero on all classes,
- `one(chi)` is the all-one marks vector,
  corresponding to the ``G``-set that consists of one point,
- `chi^n` is the `n`-th power of `chi`, for positive integers `n`.

```@docs
coordinates(chi::GAPGroupMarksVector)
restrict(chi::GAPGroupMarksVector, tbl::GAPGroupCharacterTable)
```

## The interface between tables of marks and character tables

```@docs
character_table(tom::GAPGroupTableOfMarks)
table_of_marks(tbl::GAPGroupCharacterTable)
```

## Technicalities

```@docs
GAPGroupTableOfMarks
```
