```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["group_characters.md"]
```

# Group characters

Let ``G`` be a finite group, and let ``\rho: G \rightarrow GL(n, R)``
be a group homomorphism, for some ring ``R``.
We call ``\chi: G \rightarrow R``, defined by ``\chi(g) = Trace(\rho(g))``,
the character afforded by ``\rho``.

Since ``\chi`` is constant on conjugacy classes of ``G``,
it can be represented by an array ``l`` of values such that
the value on the ``i``-th conjugacy class of ``G``
(see [`conjugacy_classes`](@ref)) is stored at ``l[i]``.
Note that this makes sense only if we assume that the ordering of
conjugacy classes of ``G`` is fixed once the classes have been computed.

We deal only with the cases that either ``R`` can be embedded into
some number field, or that ``R`` is a finite field.

In the former case, the eigenvalues of the matrix ``\rho(g)``,
for ``g \in G``,
are ``k``-th roots of unity, where ``k`` is the order of ``g``,
thus all values of ``\chi`` can be represented by elements in the
abelian closure of the field of rational numbers,
see [`abelian_closure`](@ref).
The characters obtained this way are called *ordinary characters*.

In the latter case, the list of traces of ``\rho(g)``
(the so-called *Frobenius character* of ``\rho``) is often not so interesting;
instead, one considers the *Brauer character* of ``\rho``,
which is defined on (conjugacy classes of) elements ``g`` whose order is
coprime to the characteristic of ``R``
(the so-called *``p``-regular* elements resp. classes),
by first lifting the eigenvalues of ``\rho(g)`` to complex roots of unity
and then summing up these roots;
this way, one gets again a list of values in the abelian closure of the
field of rationals.

The pointwise sum and product of two characters are again characters,
they are afforded by the direct sum and the tensor product of the
underlying representations.
A character that is not the sum of two characters is called
*absolutely irreducible*.

## Character tables

Putting the values of the absolutely irreducible ordinary characters
of a group ``G`` into an array such that the rows correspond to the characters
and the columns correspond to the conjugacy classes yields the
*ordinary character table* of ``G``, which is in fact a square matrix.
Analogously, the absolutely irreducible Brauer characters of ``G``, for a given
characteristic ``p``, yield a square matrix,
the ``p``-modular *Brauer character table*.

Ordinary character tables can be computed with [`character_table`](@ref)
from a given group.
The computation of ``p``-modular Brauer tables is currently restricted to
the case of ``p``-solvable groups.

Character tables contain a lot of information about their groups,
many questions about a finite group can be answered by computations only with
its characters.
Thus it makes sense to deal also with character tables *without* an
explicit labeling of the columns of the table by conjugacy classes of a group.
For example, the character tables shown in
Atlas of Finite Groups [CCNPW85](@cite) and from the
Atlas of Brauer Characters [JLPW95](@cite) are available in OSCAR.
Such character tables can be fetched with `character_table`
from the database, via their names.

In OSCAR, a character table `t` is identified with the array of absolutely
irreducible characters of ``G``, in the sense that `t[i]` yields the
`i`-th irreducible character of ``G``,
and `t[i, j]` is the value of this character on the `j`-th conjugacy class
of ``G`` (or the `j`-th conjugacy class of ``p``-regular elements
in the case of Brauer tables).

Ordinary and ``p``-modular Brauer tables in OSCAR are distinguished by
the field `characteristic`; its value is `0` for ordinary tables and
``p`` otherwise.

```@docs
GAPGroupCharacterTable
character_table
Base.mod(tbl::GAPGroupCharacterTable, p::Int)
all_character_table_names
```

## Attributes of group characters

```@docs
character_field
conj(chi::GAPGroupClassFunction)
Nemo.degree(chi::GAPGroupClassFunction)
indicator
is_irreducible(chi::GAPGroupClassFunction)
schur_index
```

## Attributes of character tables

```@docs
character_parameters
class_parameters
decomposition_matrix
identifier
induced_cyclic(tbl::GAPGroupCharacterTable)
is_duplicate_table
maxes
names_of_fusion_sources
class_lengths
orders_centralizers
orders_class_representatives
trivial_character(tbl::GAPGroupCharacterTable)
```

## Construct group characters from groups

```@docs
natural_character(G::PermGroup)
natural_character(G::Union{MatrixGroup{fmpq}, MatrixGroup{nf_elem}})
trivial_character(G::GAPGroup)
```

## Operations for group characters

length and iteration:

The length of a class function is the number of conjugacy classes of its
group, iteration is defined w.r.t. the ordering of conjugacy classes.

arithmetic operations:

- `chi == psi`:
  two class functions are equal if and only if they belong to the same
  character table and have the same values,
- `chi + psi` and `chi - psi` are the pointwise sum and difference,
  respectively, of the two class functions `chi`, `psi`,
- `n*chi` is the pointwise `n`-fold sum of `chi`, for an integer `n`,
- `chi*psi` is the pointwise (tensor) product of `chi` and `psi`,
- `zero(chi)` is the class function that is zero on all classes,
- `one(chi)` is the trivial character of the character table of `chi`,
- `chi^n` is the `n`-th tensor power of `chi`, for positive integers `n`,
- `chi(g)` is the value of `chi` at the element `g` of the group of `chi`,
- `chi^g` is the conjugate character of `chi` under the action of
  a group element `g` that normalizes the group ``G`` of `chi`;
  we have `chi^g(x) == chi(g*x*g^-1)` for all `x` in ``G``,
- `chi^galaut` is the Galois conjugate character of `chi` under the
  pointwise action of the field automorphism `galaut`
  (If `galaut` was created as `QQAbAutomorphism(k)` then the action raises
  each root of unity to its `k`-th power;
  this action defines a field automorphism of the `n`-th cyclotomic field
  whenever `n` and `k` are coprime.)
- `chi^tbl` is the character of the character table `tbl`
  that is induced from `chi`,
  where the group of `chi` is a subgroup of the group of `tbl`.

```@docs
scalar_product
induced_class_function
```

## Symmetrizations of group characters

```@docs
symmetrizations(characters::Vector{GAPGroupClassFunction}, n::Int)
symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)
anti_symmetric_parts(characters::Vector{GAPGroupClassFunction}, n::Int)
exterior_power(chi::GAPGroupClassFunction, n::Int)
symmetric_power(chi::GAPGroupClassFunction, n::Int)
orthogonal_components(characters::Vector{GAPGroupClassFunction}, n::Int)
symplectic_components(characters::Vector{GAPGroupClassFunction}, n::Int)
```

## Operations for character tables

```@docs
class_multiplication_coefficient
known_class_fusion
order(tbl::GAPGroupCharacterTable)
possible_class_fusions
```

## Character tables and normal subgroups

Normal subgroups of a group ``G`` are unions of conjugacy classes of elements
of ``G``.
Thus one can often turn questions about a normal subgroup ``N`` of ``G``
into questions about the array of those positions in the list of
conjugacy classes of ``G`` that contain the elements of ``N``.

```@docs
class_positions_of_kernel
class_positions_of_pcore
```
