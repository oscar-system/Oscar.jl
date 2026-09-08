# Notes for GAP users

OSCAR uses GAP for most of its group theory,
via the Julia package [GAP.jl](https://github.com/oscar-system/GAP.jl).
The GAP session runs inside the Julia session, and all GAP functions and
packages are available from it.
Only a part of GAP's functionality has a native OSCAR interface,
but the rest can be called directly, as described in the first section.
The remaining sections describe differences between GAP and OSCAR,
and list GAP functions together with their OSCAR counterparts.

## Calling GAP functions from OSCAR

### Global GAP variables and functions

Every global GAP variable is available as `GAP.Globals.<name>`,
and GAP functions can be called like Julia functions.
Julia integers and booleans can be passed as arguments directly,
all other arguments must be converted with `GapObj` first, see below.
Small integers, booleans, and elements of small finite fields that GAP
returns are converted to Julia automatically.
Most other results are `GapObj`s and get printed with the prefix `GAP:`
(but a GAP function may also return any Julia object, for example one
that was passed to it as an argument).
```jldoctest
julia> GAP.Globals.Factorial(5)
120

julia> GAP.Globals.IsPrimeInt(97)
true

julia> x = GAP.Globals.Factorial(30)
GAP: 265252859812191058636308480000000
```

### Converting GAP objects to OSCAR objects

For numbers, lists, and strings, use the constructor of the Julia type
that you want to get.
Nested lists need the `recursive` keyword argument, and OSCAR matrices are
constructed from GAP matrices with `matrix`.
```jldoctest
julia> ZZ(GAP.Globals.Factorial(30))
265252859812191058636308480000000

julia> Vector{Int}(GAP.Globals.DivisorsInt(12))
6-element Vector{Int64}:
  1
  2
  3
  4
  6
 12

julia> Vector{Vector{Int}}(GAP.Globals.Partitions(4))
5-element Vector{Vector{Int64}}:
 [1, 1, 1, 1]
 [2, 1, 1]
 [2, 2]
 [3, 1]
 [4]

julia> String(GAP.Globals.StructureDescription(GAP.Globals.SymmetricGroup(4)))
"S4"

julia> matrix(ZZ, GAP.Globals.IdentityMat(2))
[1   0]
[0   1]
```
For example, here is how to use GAP's `Derangements`,
for which OSCAR provides no counterpart.
```jldoctest
julia> x = GAP.Globals.Derangements(GapObj(1:4))
GAP: [ [ 2, 1, 4, 3 ], [ 2, 3, 4, 1 ], [ 2, 4, 1, 3 ], [ 3, 1, 4, 2 ], [ 3, 4, 1, 2 ], [ 3, 4, 2, 1 ], [ 4, 1, 2, 3 ], [ 4, 3, 1, 2 ], [ 4, 3, 2, 1 ] ]

julia> Vector{Vector{Int}}(x)
9-element Vector{Vector{Int64}}:
 [2, 1, 4, 3]
 [2, 3, 4, 1]
 [2, 4, 1, 3]
 [3, 1, 4, 2]
 [3, 4, 1, 2]
 [3, 4, 2, 1]
 [4, 1, 2, 3]
 [4, 3, 1, 2]
 [4, 3, 2, 1]
```
GAP groups are turned into OSCAR groups by the constructors `PermGroup`,
`PcGroup`, `FPGroup`, and `matrix_group`.
```jldoctest
julia> G = PermGroup(GAP.Globals.MathieuGroup(11))
Permutation group of degree 11 and order 7920

julia> H = PcGroup(GAP.Globals.SmallGroup(8, 3))
Pc group of order 8
```
Elements of GAP's finite fields and cyclotomic fields are converted by
calling the OSCAR field on them, and the same works for matrices.
```jldoctest
julia> F = GF(3);

julia> F(GAP.Globals.Z(3))
2

julia> K, z = cyclotomic_field(5);

julia> K(GAP.Globals.E(5))
z_5

julia> matrix(K, GAP.evalstr("[[E(5), 1], [0, E(5)^2]]"))
[z_5       1]
[  0   z_5^2]
```
For more complex situations, `Oscar.iso_gap_oscar` and
`Oscar.iso_oscar_gap` construct an isomorphism between a GAP ring or
field and its OSCAR counterpart,
see the section on [GAP Integration](@ref).
More conversions are described in the
[GAP.jl manual](https://oscar-system.github.io/GAP.jl/stable/conversion/).

### Passing OSCAR objects to GAP

The function `GapObj` converts Julia and OSCAR objects into GAP objects.
For an OSCAR group `G`, `GapObj(G)` returns the GAP group that it wraps,
and similarly for group elements, character tables,
and other objects that are backed by GAP.
```jldoctest
julia> G = symmetric_group(4);

julia> GAP.Globals.StructureDescription(GapObj(G))
GAP: "S4"

julia> GAP.Globals.IsAbelian(GapObj(G))
false

julia> GapObj(ZZ(2)^100)
GAP: 1267650600228229401496703205376

julia> GapObj([1, 2, 3])
GAP: [ 1, 2, 3 ]

julia> GapObj([[1, 2], [3, 4]]; recursive = true)
GAP: [ [ 1, 2 ], [ 3, 4 ] ]

julia> GapObj(matrix(QQ, [1 2; 3 4]))
GAP: [ [ 1, 2 ], [ 3, 4 ] ]
```
Julia functions can be passed to GAP functions that expect a function
argument.
```jldoctest
julia> GAP.Globals.Filtered(GapObj(1:10), is_prime)
GAP: [ 2, 3, 5, 7 ]
```

### Running GAP code

`GAP.evalstr` evaluates a string containing GAP code and returns the
result.
```jldoctest
julia> GAP.evalstr("List([1..5], x -> x^2)")
GAP: [ 1, 4, 9, 16, 25 ]

julia> r = GAP.evalstr("rec(a := 1, b := [1, 2])");

julia> r.a
1
```
`GAP.prompt()` opens a GAP prompt inside the Julia session.
In it, the variables of the Julia session are available via `Julia.<name>`,
and the OSCAR functions via `Oscar_jl.<name>`.
Enter `quit;` to get back to the Julia prompt.
GAP's help system is available at the Julia prompt via
`?GAP.Globals.Size` or `GAP.show_gap_help("Size", true)`.

GAP packages that are not loaded automatically can be loaded with
`GAP.Packages.load`,
and packages that are not shipped with OSCAR can be installed with
`GAP.Packages.install`.
```jldoctest
julia> GAP.Packages.load("ctbllib")
true
```

### When OSCAR does not have the function you need

If you cannot find a counterpart of a GAP function in OSCAR
(the [tables below](@ref "Common GAP functions and their OSCAR counterparts")
may help), calling the GAP function as described above will usually
work as long as the objects involved can be converted.
In this case, please
[open an issue](https://github.com/oscar-system/Oscar.jl/issues)
to request a proper OSCAR interface for it.

The developer documentation mentions `Oscar.GAPWrap`;
this is intended for OSCAR's own code, not for interactive use.


## Differences in syntax

- In GAP, equality of two objects is checked with `=`,
  and one assigns a value to a variable with `:=`.
  In Julia, equality is checked with `==`,
  and `=` denotes assignment.
  Similarly, inequality of objects is checked with `<>` in GAP
  and with `!=` in Julia.

- In GAP, the operator `not` is used to negate boolean expressions,
  whereas `!` is used in Julia.
  The operators `and` and `or` are `&&` and `||` in Julia.

- In GAP, object identity is checked with the function `IsIdenticalObj`,
  whereas the infix operator `===` (with negation `!==`)
  is used in Julia.

- In GAP, `if` statements have the form
  ```gap
  if condition1 then
    statements1
  elif condition2 then
    statements2
  else
    statements3
  fi;
  ```
  whereas the Julia syntax is
  ```julia
  if condition1
    statements1
  elseif condition2
    statements2
  else
    statements3
  end
  ```
  Similarly, GAP's `for` loops have the form
  ```gap
  for var in list do
    statements
  od;
  ```
  whereas the Julia syntax is
  ```julia
  for var in list
    statements
  end
  ```
  (The situation with `while` loops is analogous.)

- GAP functions are defined as `f := function(x) ... end;`,
  the Julia equivalents are `function f(x) ... end`
  and, for one-liners, `f(x) = ...`.
  Anonymous functions look the same in both languages: `x -> x^2`.

- GAP's `G.1` for the first generator of a group is `G[1]` in OSCAR.

- The `;` at the end of a statement is mandatory in GAP,
  and two semicolons suppress the output of the value.
  In Julia, no semicolon is needed, and a trailing `;` suppresses the
  output of the value in an interactive session,
  see [Semicolons and output](@ref).
  So in both languages, an extra semicolon suppresses output.

## Differences in semantics

- **Integer literals are machine integers.**
  In GAP, all integers are arbitrary precision integers.
  In Julia, `2^100` evaluates to `0`,
  see [Integers and rational numbers](@ref other_integers).

- **The sum of a matrix and a scalar is different.**
  In GAP, the sum of a matrix (a list of lists) and a scalar is defined
  recursively as the pointwise sum.
  ```gap-repl
  gap> [ [ 1, 2 ], [ 3, 4 ] ] + 2;
  [ [ 3, 4 ], [ 5, 6 ] ]
  ```
  In OSCAR, the sum of a matrix and a scalar is defined as the sum of
  the given matrix and the multiple of the identity matrix that is given
  by the scalar.
  ```jldoctest
  julia> matrix(ZZ, [1 2; 3 4]) + 2
  [3   2]
  [3   6]
  ```

- **There are no natural embeddings.**
  GAP provides natural embeddings of many algebraic structures.
  For example, two finite fields of the same characteristic are embedded
  into each other whenever this makes sense, and the elements of the smaller
  field are regarded also as elements of the larger field.
  Analogously, subfields of cyclotomic fields are naturally embedded
  into each other, and in fact their elements are internally represented
  w.r.t. the smallest possible cyclotomic field.

  In OSCAR, this is not the case.
  Each element of an algebraic structure has a parent,
  and operations involving several elements (such as arithmetic operations)
  are usually restricted to the situation that their parents coincide.
  One has to explicitly coerce a given element into a different parent
  if necessary, see [Every object has a parent](@ref).

- **Permutations belong to a symmetric group of fixed degree.**
  This is a consequence of the previous point.
  `cperm([1, 2, 3])` creates the permutation `(1,2,3)` as an element of
  the symmetric group of degree 3,
  and it cannot be multiplied with an element of the symmetric group of
  degree 4.
  Create permutations as elements of the group `G` you want to work in,
  for example with `cperm(G, [1, 2, 3])` or `perm(G, [2, 3, 1])`.

  Similarly, each permutation group in OSCAR has a fixed degree,
  and the function `is_transitive` checks whether its argument is transitive
  on the points from 1 to the degree.
  In GAP, however, the function `IsTransitive`, called with a permutation
  group, checks whether this group is transitive on the points which are
  moved by it.
  Thus the group generated by the permutation `(1, 2, 4)` is regarded as
  transitive in GAP but as intransitive in OSCAR.

  For the same reason, `transitive_group` in OSCAR supports degree 1,
  where the trivial group is the unique transitive group,
  whereas GAP's library of transitive groups starts at degree 2.

- **Subgroups come together with embeddings.**
  Functions that return a subgroup, such as `sub`, `center`,
  `derived_subgroup`, `sylow_subgroup`, `stabilizer`, and `kernel`,
  return a tuple consisting of the subgroup and its embedding into the
  given group, and `quo` returns the quotient group together with the
  natural projection,
  see [Many constructors return more than one object](@ref).

- **The argument order of subset tests is reversed.**
  GAP's `IsSubset(G, H)` and `IsSubgroup(G, H)` put the larger object
  first.
  OSCAR follows the Julia convention that the arguments of `issubset`,
  `is_subset`, `is_subgroup`, and `is_normal_subgroup` appear in the
  same order as in the mathematical notation ``H \subseteq G``,
  that is, `is_subgroup(H, G)`.

- **`Size` has several counterparts.**
  Depending on the object, GAP's `Size` corresponds to `order(G)` for a
  group, `length(l)` for a list or another collection,
  and `number_of_rows(M)`, `number_of_columns(M)` for a matrix.
  Julia's `size(M)` returns the tuple of dimensions of a matrix.

- **Global variables are not protected.**
  Global OSCAR variables are not write protected,
  contrary to most global GAP variables.
  Thus there is always the danger that assignments overwrite Julia functions.
  For example, it is tempting to use `gens`, `hom`, and `map` as names for
  variables, but Julia or OSCAR define them already.

  (Also copying some lines of code from an OSCAR function into a Julia session
  can be dangerous in this sense,
  because some names of local variables of the function may coincide with the
  names of global variables.)

## Interactive sessions

When an error occurs or when the user hits ctrl-C in a GAP session,
usually a break loop is entered,
from which one can either try to continue the computations,
by entering `return`, or return to the GAP prompt, by entering `quit`;
in the latter case, some objects may be corrupted afterwards.

In a Julia session, one gets automatically back to the Julia prompt
when an error occurs or when the user hits ctrl-C,
and again some objects may be corrupted afterwards.

## Names of functions and variables

Variable names in GAP and Julia are recommended to be written in
camel case and snake case, respectively, see [Naming conventions](@ref).
For example, the GAP function `SylowSubgroup` corresponds to
OSCAR's `sylow_subgroup`.
Guessing the OSCAR name this way works surprisingly often;
the tables below list cases where it does not.

The GAP rule that the names of user variables should start with a
lowercase letter, in order to avoid clashes with system variables,
does not make sense in Julia.

## Common GAP functions and their OSCAR counterparts

### Lists and functional programming

| GAP | OSCAR |
|:----|:------|
| `Length(l)`, `Size(l)` | `length(l)` |
| `[1..10]` | `1:10` (`collect(1:10)` for a vector) |
| `l[i]`, `l{[i, j]}` | `l[i]`, `l[[i, j]]` |
| `List(l, f)` | `map(f, l)` or `[f(x) for x in l]` |
| `Filtered(l, f)` | `filter(f, l)` or `[x for x in l if f(x)]` |
| `ForAll(l, f)`, `ForAny(l, f)` | `all(f, l)`, `any(f, l)` |
| `Number(l, f)` | `count(f, l)` |
| `Position(l, x)`, `PositionProperty(l, f)` | `findfirst(==(x), l)`, `findfirst(f, l)` |
| `Sum(l)`, `Product(l)` | `sum(l)`, `prod(l)` |
| `Maximum(l)`, `Minimum(l)` | `maximum(l)`, `minimum(l)` |
| `Add(l, x)`, `Append(l, m)` | `push!(l, x)`, `append!(l, m)` |
| `Concatenation(l, m)` | `vcat(l, m)` |
| `Reversed(l)` | `reverse(l)` |
| `Sort(l)`, `SortedList(l)` | `sort!(l)`, `sort(l)` |
| `Set(l)` | `sort(unique(l))` (`Set(l)` creates a Julia set) |
| `Union(a, b)`, `Intersection(a, b)`, `Difference(a, b)` | `union(a, b)`, `intersect(a, b)`, `setdiff(a, b)` |
| `IsSubset(a, b)` | `issubset(b, a)` |
| `IsEmpty(l)` | `isempty(l)` |
| `Cartesian(a, b)` | `Iterators.product(a, b)` |
| `Combinations(l, k)` | `combinations(l, k)` |
| `Partitions(n)` | `partitions(n)` |
| `ShallowCopy(x)`, `StructuralCopy(x)` | `copy(x)`, `deepcopy(x)` |
| `rec(a := 1)`, `r.a` | `Dict(:a => 1)`, `r[:a]` |
| `IsBound(x)` | `@isdefined x` |
| `Print(x)`, `Display(x)` | `print(x)`, `display(x)` |
| `Read("file.g")` | `include("file.jl")` |
| `Random(l)` | `rand(l)` |

### Integers

| GAP | OSCAR |
|:----|:------|
| `Factorial(n)` | `factorial(ZZ(n))` |
| `Binomial(n, k)` | `binomial(n, k)` |
| `Gcd(a, b)`, `Lcm(a, b)`, `Gcdex(a, b)` | `gcd(a, b)`, `lcm(a, b)`, `gcdx(a, b)` |
| `QuoInt(a, b)`, `RemInt(a, b)`, `a mod b` | `div(a, b)`, `rem(a, b)`, `mod(a, b)` |
| `PowerMod(a, e, m)` | `powermod(a, e, m)` |
| `IsPrimeInt(n)`, `NextPrimeInt(n)` | `is_prime(n)`, `next_prime(n)` |
| `Factors(n)` | `factor(n)` |
| `DivisorsInt(n)` | `divisors(n)` |
| `Phi(n)` | `euler_phi(n)` |
| `Jacobi(a, n)` | `jacobi_symbol(a, n)` |
| `ChineseRem([m1, m2], [r1, r2])` | `crt([r1, r2], [m1, m2])` |
| `RootInt(n, k)` | `iroot(n, k)` |
| `Fibonacci(n)`, `Bell(n)` | `fibonacci(n)`, `bell(n)` |

### Groups

| GAP | OSCAR |
|:----|:------|
| `SymmetricGroup(n)`, `AlternatingGroup(n)` | `symmetric_group(n)`, `alternating_group(n)` |
| `CyclicGroup(n)`, `DihedralGroup(n)` | `cyclic_group(n)`, `dihedral_group(n)` |
| `AbelianGroup([2, 4])` | `abelian_group([2, 4])` |
| `SmallGroup(n, i)`, `IdGroup(G)` | `small_group(n, i)`, `small_group_identification(G)` |
| `Group(g, h)` (permutations) | `permutation_group(n, [g, h])` |
| `Group(m1, m2)` (matrices) | `matrix_group([m1, m2])` |
| `GL(n, q)`, `SL(n, q)` | `GL(n, q)`, `SL(n, q)` |
| `FreeGroup(2)`, `F / rels` | `free_group(2)`, `quo(F, rels)` |
| `DirectProduct(G, H)` | `direct_product(G, H)` |
| `Subgroup(G, gens)` | `sub(G, gens)` |
| `GeneratorsOfGroup(G)`, `G.1` | `gens(G)`, `G[1]` |
| `One(G)` | `one(G)` |
| `Size(G)`, `Order(g)` | `order(G)`, `order(g)` |
| `Elements(G)` | `collect(G)` |
| `Random(G)` | `rand(G)` |
| `Exponent(G)` | `exponent(G)` |
| `Centre(G)`, `DerivedSubgroup(G)` | `center(G)`, `derived_subgroup(G)` |
| `FittingSubgroup(G)`, `FrattiniSubgroup(G)`, `Socle(G)` | `fitting_subgroup(G)`, `frattini_subgroup(G)`, `socle(G)` |
| `Centralizer(G, x)`, `Normalizer(G, H)` | `centralizer(G, x)`, `normalizer(G, H)` |
| `SylowSubgroup(G, p)` | `sylow_subgroup(G, p)` |
| `Intersection(G, H)` | `intersect(G, H)` |
| `Index(G, H)` | `index(G, H)` |
| `IsSubgroup(G, H)`, `IsNormal(G, H)` | `is_subgroup(H, G)`, `is_normal_subgroup(H, G)` |
| `NormalSubgroups(G)`, `MaximalSubgroups(G)` | `normal_subgroups(G)`, `maximal_subgroups(G)` |
| `ConjugacyClassesSubgroups(G)` | `subgroup_classes(G)` |
| `IsAbelian(G)`, `IsSolvable(G)`, `IsNilpotent(G)`, `IsSimple(G)`, `IsPerfect(G)` | `is_abelian(G)`, `is_solvable(G)`, `is_nilpotent(G)`, `is_simple(G)`, `is_perfect(G)` |
| `ConjugacyClasses(G)`, `NrConjugacyClasses(G)` | `conjugacy_classes(G)`, `number_of_conjugacy_classes(G)` |
| `Representative(C)` | `representative(C)` |
| `IsConjugate(G, x, y)` | `is_conjugate(G, x, y)` |
| `CharacterTable(G)`, `CharacterTable("M11")` | `character_table(G)`, `character_table("M11")` |
| `TableOfMarks(G)` | `table_of_marks(G)` |
| `AutomorphismGroup(G)` | `automorphism_group(G)` |
| `StructureDescription(G)` | `describe(G)` |
| `AbelianInvariants(G)` | `abelian_invariants(G)` |
| `IsomorphismGroups(G, H)` | `isomorphism(G, H)` |
| `IsomorphismPermGroup(G)`, `IsomorphismPcGroup(G)`, `IsomorphismFpGroup(G)` | `isomorphism(PermGroup, G)`, `isomorphism(PcGroup, G)`, `isomorphism(FPGroup, G)` |
| `GroupHomomorphismByImages(G, H, gens, imgs)` | `hom(G, H, gens, imgs)` |
| `Image(f, x)`, `x^f` | `f(x)` |
| `Image(f)`, `Kernel(f)` | `image(f)`, `kernel(f)` |
| `PreImagesRepresentative(f, y)` | `preimage(f, y)` |
| `NaturalHomomorphismByNormalSubgroup(G, N)`, `FactorGroup(G, N)` | `quo(G, N)` |
| `Orbit(G, x)`, `Orbits(G)` | `orbit(G, x)`, `orbits(G)` |
| `Stabilizer(G, x)` | `stabilizer(G, x)` |
| `IsTransitive(G)`, `IsPrimitive(G)` | `is_transitive(G)`, `is_primitive(G)` |

### Permutations

| GAP | OSCAR |
|:----|:------|
| `(1,2,3)` | `cperm(G, [1, 2, 3])` |
| `PermList([2, 3, 1])` | `perm(G, [2, 3, 1])` |
| `ListPerm(g)` | `Vector(g)` |
| `i^g`, `OnPoints(i, g)` | `i^g` |
| `g^h`, `Comm(g, h)` | `g^h`, `comm(g, h)` |
| `SignPerm(g)`, `CycleStructurePerm(g)` | `sign(g)`, `cycle_structure(g)` |
| `MovedPoints(g)`, `NrMovedPoints(g)` | `moved_points(g)`, `number_of_moved_points(g)` |

### Rings, fields, and polynomials

| GAP | OSCAR |
|:----|:------|
| `Integers`, `Rationals` | `ZZ`, `QQ` |
| `GF(q)`, `GF(p, n)` | `GF(q)`, `GF(p, n)` |
| `Integers mod n` | `residue_ring(ZZ, n)` |
| `CF(n)` | `cyclotomic_field(n)` |
| `E(n)` | `K, z = abelian_closure(QQ); z(n)` |
| `PolynomialRing(Rationals, ["x", "y"])` | `polynomial_ring(QQ, [:x, :y])` |
| `Indeterminate(Rationals, "x")` | `polynomial_ring(QQ, :x)` |
| `Value(f, x)` | `evaluate(f, x)` or `f(x)` |
| `Degree(f)`, `Derivative(f)` | `degree(f)`, `derivative(f)` |
| `Factors(f)`, `IsIrreducible(f)` | `factor(f)`, `is_irreducible(f)` |
| `RootsOfUPol(f)` | `roots(f)` |
| `CoefficientsOfUnivariatePolynomial(f)` | `coefficients(f)` |
| `Gcd(f, g)`, `Resultant(f, g)`, `Discriminant(f)` | `gcd(f, g)`, `resultant(f, g)`, `discriminant(f)` |
| `AlgebraicExtension(Rationals, f)` | `number_field(f)` |
| `Ideal(R, [f, g])` | `ideal(R, [f, g])` |
| `GroebnerBasis(I, ord)` | `groebner_basis(I)` |

### Matrices

| GAP | OSCAR |
|:----|:------|
| `[[1, 2], [3, 4]]` | `matrix(ZZ, [1 2; 3 4])` |
| `IdentityMat(n)`, `NullMat(m, n)`, `DiagonalMat(l)` | `identity_matrix(ZZ, n)`, `zero_matrix(ZZ, m, n)`, `diagonal_matrix(l)` |
| `M[i][j]`, `M[i, j]` | `M[i, j]` |
| `M[i]` | `M[i, :]` |
| `NrRows(M)`, `NrCols(M)` | `nrows(M)`, `ncols(M)` |
| `TransposedMat(M)` | `transpose(M)` |
| `DeterminantMat(M)`, `TraceMat(M)`, `RankMat(M)` | `det(M)`, `tr(M)`, `rank(M)` |
| `Inverse(M)`, `M^-1` | `inv(M)`, `M^-1` |
| `NullspaceMat(M)` | `kernel(M)` |
| `SolutionMat(M, v)` | `solve(M, v; side = :left)` |
| `CharacteristicPolynomial(M)`, `MinimalPolynomial(M)` | `charpoly(M)`, `minpoly(M)` |
| `Eigenvalues(F, M)` | `eigenvalues(M)` |
| `SmithNormalFormIntegerMat(M)`, `HermiteNormalFormIntegerMat(M)` | `snf(M)`, `hnf(M)` |
| `TriangulizedMat(M)` | `rref(M)` |
| `KroneckerProduct(M, N)` | `kronecker_product(M, N)` |
