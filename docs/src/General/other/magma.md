# Notes for Magma users

There is no interface between Magma and OSCAR, so this page is only
about translating what you know from Magma to OSCAR.
Many of the number theoretic functions in OSCAR live in the package
[Hecke](https://github.com/thofma/Hecke.jl), whose function names
often are the `snake_case` versions of the Magma names.

!!! note "Help wanted"
    This page is a start. Please tell us what is missing, see
    [Notes for users of other computer algebra systems](@ref).

## Differences in syntax

- Assignment is `:=` in Magma and `=` in Julia.
  The comparison operators `eq`, `ne`, `lt`, `le`, `gt`, `ge` are
  `==`, `!=`, `<`, `<=`, `>`, `>=` in Julia,
  and `and`, `or`, `not` are `&&`, `||`, `!`.

- Every Magma statement ends with `;`.
  In an interactive Julia session, a trailing `;` suppresses the output
  of the value, see [Semicolons and output](@ref).

- Magma's `if ... then ... elif ... else ... end if;` is
  `if ... elseif ... else ... end` in Julia,
  and `for x in L do ... end for;` is `for x in L ... end`.
  The same holds for `while` loops.
  There are no `then` and `do` keywords.

- Functions are defined with `function f(x) ... end` in Julia;
  one-liners can be written as `f(x) = ...`.
  Procedures that modify their arguments (`~x` in Magma) are ordinary
  functions in Julia, by convention their names end with `!`.

- The cardinality `#S` of a set, sequence, or group is `length(S)` for
  collections and `order(G)` for groups.

- Sequence and set constructors translate to comprehensions:
  `[ f(x) : x in L | c(x) ]` is `[f(x) for x in L if c(x)]`,
  and `{ ... }` is `Set(...)` applied to a comprehension.
  The reductions `&+L` and `&*L` are `sum(L)` and `prod(L)`,
  `&cat L` is `reduce(vcat, L)`.
  `forall{ x : x in L | c(x) }` and `exists{ ... }` are
  `all(c, L)` and `any(c, L)`.

- Coercion `K!x` is `K(x)` in OSCAR.

- Strings are concatenated with `*` instead of `cat`.

- Comments start with `#` instead of `//`; multi-line comments are
  `#= ... =#` instead of `/* ... */`.

- `time f(x);` is `@time f(x)`, and `load "file";` is `include("file.jl")`.
  `assigned x` is `@isdefined x`.

- Magma's `?` help syntax is available in Julia as `?name`.

- Many Magma constructors use angle brackets:
  `sub<G | a, b>`, `quo<G | N>`, `hom<G -> H | ...>`, `ideal<R | f, g>`,
  `R<x, y> := PolynomialRing(K, 2)`.
  The OSCAR counterparts are ordinary functions:
  `sub(G, [a, b])`, `quo(G, N)`, `hom(G, H, ...)`, `ideal(R, [f, g])`,
  and `R, (x, y) = polynomial_ring(K, [:x, :y])`.

## Differences in semantics

- **Integer literals are machine integers.**
  In Magma, `2^100` is computed exactly and `3/4` is a rational number.
  In Julia, `2^100` evaluates to `0` and `3/4` to `0.75`.
  Write `ZZ(2)^100` and `QQ(3, 4)` instead,
  see [Integers and rational numbers](@ref other_integers).

- **Multiple return values are tuples.**
  Both Magma and Julia allow functions to return several values,
  but in Julia these form a tuple that is printed as a whole,
  and there are no "optional" return values that are silently dropped.
  Where a Magma predicate returns a witness as a second value,
  OSCAR has a separate function whose name says so:
  ```jldoctest
  julia> is_square(ZZ(16))
  true

  julia> is_square_with_sqrt(ZZ(16))
  (true, 4)
  ```
  For example, `H, f := sub<G | ...>` returns the subgroup together with
  its embedding, and `H := sub<G | ...>` silently drops the embedding.
  OSCAR's `sub` always returns both,
  see [Many constructors return more than one object](@ref).

- **The order of dihedral groups.**
  Magma's `DihedralGroup(n)` is the dihedral group of order ``2n``,
  OSCAR's `dihedral_group(n)` is the dihedral group of order ``n``,
  following GAP.

- **The default monomial ordering differs.**
  `PolynomialRing(K, n)` in Magma uses the lexicographical ordering,
  and Gröbner bases are computed w.r.t. it.
  OSCAR's multivariate polynomial rings use the degree reverse
  lexicographical ordering by default, and `groebner_basis` takes the
  ordering as a keyword argument.
  ```jldoctest
  julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

  julia> I = ideal(R, [x^2 + y^2 - 1, x - y]);

  julia> groebner_basis(I; ordering = lex(R))
  Gröbner basis with elements
    1: 2*y^2 - 1
    2: x - y
  with respect to the ordering
    lex([x, y])
  ```

- **Argument order of `ChangeRing`.**
  Magma's `ChangeRing(M, R)` is `change_base_ring(R, M)` in OSCAR;
  the ring comes first, as in `matrix(R, ...)`.

- **Types are not intrinsics.**
  Magma's `Type(x)` is `typeof(x)` in Julia, `Parent(x)` is `parent(x)`.
  There is no `Category`.
  The signatures of a function `f` are listed by `methods(f)`.

## Common Magma functions and their OSCAR counterparts

### Sequences and sets

| Magma | OSCAR |
|:------|:------|
| `#L` | `length(L)` |
| `[1..10]`, `[1..10 by 2]` | `1:10`, `1:2:10` |
| `Append(~L, x)` | `push!(L, x)` |
| `L cat M` | `vcat(L, M)` |
| `Reverse(L)`, `Sort(L)` | `reverse(L)`, `sort(L)` |
| `Position(L, x)`, `Index(L, x)` | `findfirst(==(x), L)` |
| `Max(L)`, `Min(L)` | `maximum(L)`, `minimum(L)` |
| `IsEmpty(L)` | `isempty(L)` |
| `Seqset(L)`, `Setseq(S)` | `Set(L)`, `collect(S)` |
| `A join B`, `A meet B`, `A diff B` | `union(A, B)`, `intersect(A, B)`, `setdiff(A, B)` |
| `Include(~S, x)`, `Exclude(~S, x)` | `push!(S, x)`, `delete!(S, x)` |
| `x in S`, `x notin S` | `x in S`, `!(x in S)` |
| `Universe(L)` | `eltype(L)` |
| `ChangeUniverse(L, R)` | `R.(L)` |

### Integers and rational numbers

| Magma | OSCAR |
|:------|:------|
| `Integers()`, `Rationals()` | `ZZ`, `QQ` |
| `a div b`, `a mod b` | `div(a, b)`, `mod(a, b)` |
| `Factorization(n)` | `factor(n)` |
| `IsPrime(n)`, `NextPrime(n)` | `is_prime(n)`, `next_prime(n)` |
| `Divisors(n)`, `EulerPhi(n)` | `divisors(n)`, `euler_phi(n)` |
| `GCD(a, b)`, `LCM(a, b)` | `gcd(a, b)`, `lcm(a, b)` |
| `Binomial(n, k)`, `Factorial(n)` | `binomial(n, k)`, `factorial(ZZ(n))` |
| `Numerator(q)`, `Denominator(q)` | `numerator(q)`, `denominator(q)` |
| `Floor(q)`, `Ceiling(q)`, `Round(q)` | `floor(q)`, `ceil(q)`, `round(q)` |
| `Isqrt(n)`, `IsSquare(n)` | `isqrt(n)`, `is_square(n)`, `is_square_with_sqrt(n)` |
| `Integers(n)` | `residue_ring(ZZ, n)` |

### Groups

| Magma | OSCAR |
|:------|:------|
| `Sym(n)`, `Alt(n)` | `symmetric_group(n)`, `alternating_group(n)` |
| `CyclicGroup(n)`, `DihedralGroup(n)` | `cyclic_group(n)`, `dihedral_group(2n)` |
| `AbelianGroup([2, 4])` | `abelian_group([2, 4])` |
| `SmallGroup(n, i)`, `IdentifyGroup(G)` | `small_group(n, i)`, `small_group_identification(G)` |
| `PermutationGroup<n \| g, h>` | `permutation_group(n, [g, h])` |
| `Sym(n)!(1,2,3)` | `cperm(G, [1, 2, 3])` |
| `sub<G \| a, b>` | `sub(G, [a, b])` |
| `quo<G \| N>` | `quo(G, N)` |
| `hom<G -> H \| a :-> x, b :-> y>` | `hom(G, H, [a, b], [x, y])` |
| `f(x)`, `x @ f` | `f(x)` |
| `Kernel(f)`, `Image(f)` | `kernel(f)`, `image(f)` |
| `#G`, `Order(G)`, `Order(g)` | `order(G)`, `order(g)` |
| `Generators(G)`, `G.1` | `gens(G)`, `G[1]` |
| `Random(G)` | `rand(G)` |
| `Centre(G)`, `DerivedSubgroup(G)` | `center(G)`, `derived_subgroup(G)` |
| `Centralizer(G, x)`, `Normalizer(G, H)` | `centralizer(G, x)`, `normalizer(G, H)` |
| `Sylow(G, p)` | `sylow_subgroup(G, p)` |
| `NormalSubgroups(G)` | `normal_subgroups(G)` |
| `IsAbelian(G)`, `IsSoluble(G)`, `IsSimple(G)` | `is_abelian(G)`, `is_solvable(G)`, `is_simple(G)` |
| `IsIsomorphic(G, H)` | `is_isomorphic(G, H)`, `isomorphism(G, H)` |
| `Classes(G)`, `ConjugacyClasses(G)` | `conjugacy_classes(G)` |
| `CharacterTable(G)` | `character_table(G)` |
| `AutomorphismGroup(G)` | `automorphism_group(G)` |
| `GroupName(G)` | `describe(G)` |
| `Orbit(G, x)`, `Stabilizer(G, x)` | `orbit(G, x)`, `stabilizer(G, x)` |

### Polynomials and ideals

| Magma | OSCAR |
|:------|:------|
| `R<x, y> := PolynomialRing(K, 2)` | `R, (x, y) = polynomial_ring(K, [:x, :y])` |
| `P<t> := PolynomialRing(K)` | `P, t = polynomial_ring(K, :t)` |
| `Evaluate(f, [1, 2])` | `evaluate(f, [1, 2])` or `f(1, 2)` |
| `Degree(f)`, `TotalDegree(f)` | `degree(f)`, `total_degree(f)` |
| `Derivative(f, x)` | `derivative(f, x)` |
| `Coefficients(f)`, `Monomials(f)`, `Terms(f)` | `coefficients(f)`, `monomials(f)`, `terms(f)` |
| `LeadingCoefficient(f)`, `LeadingTerm(f)` | `leading_coefficient(f)`, `leading_term(f)` |
| `Factorization(f)`, `IsIrreducible(f)` | `factor(f)`, `is_irreducible(f)` |
| `Roots(f)` | `roots(f)` |
| `GCD(f, g)`, `Resultant(f, g)`, `Discriminant(f)` | `gcd(f, g)`, `resultant(f, g)`, `discriminant(f)` |
| `ideal<R \| f, g>` | `ideal(R, [f, g])` |
| `GroebnerBasis(I)` | `groebner_basis(I)` |
| `Dimension(I)`, `Radical(I)` | `dim(I)`, `radical(I)` |
| `PrimaryDecomposition(I)` | `primary_decomposition(I)` |
| `EliminationIdeal(I, {x})` | `eliminate(I, [x])` |
| `NormalForm(f, I)`, `f in I` | `normal_form(f, I)`, `f in I` |
| `Variety(I)` | `rational_solutions(I)` |
| `quo<R \| I>` | `quo(R, I)` |

### Number fields

| Magma | OSCAR |
|:------|:------|
| `K<a> := NumberField(f)` | `K, a = number_field(f, "a")` |
| `CyclotomicField(n)`, `QuadraticField(d)` | `cyclotomic_field(n)`, `quadratic_field(d)` |
| `MaximalOrder(K)`, `RingOfIntegers(K)` | `maximal_order(K)` |
| `Discriminant(K)`, `Degree(K)` | `discriminant(K)`, `degree(K)` |
| `MinimalPolynomial(a)`, `Norm(a)`, `Trace(a)` | `minpoly(a)`, `norm(a)`, `tr(a)` |
| `ClassGroup(O)`, `ClassNumber(K)` | `class_group(O)`, `class_number(K)` |
| `UnitGroup(O)` | `unit_group(O)` |
| `Basis(O)` | `basis(O)` |
| `ideal<O \| 2>` | `ideal(O, 2)` |
| `Factorization(I)`, `IsPrime(I)` | `factor(I)`, `is_prime(I)` |
| `GaloisGroup(f)` | `galois_group(K)` |

### Matrices

| Magma | OSCAR |
|:------|:------|
| `Matrix(K, 2, 2, [1, 2, 3, 4])` | `matrix(K, 2, 2, [1, 2, 3, 4])` or `matrix(K, [1 2; 3 4])` |
| `IdentityMatrix(K, n)`, `ZeroMatrix(K, m, n)` | `identity_matrix(K, n)`, `zero_matrix(K, m, n)` |
| `Nrows(M)`, `Ncols(M)` | `nrows(M)`, `ncols(M)` |
| `M[i, j]`, `M[i]` | `M[i, j]`, `M[i, :]` |
| `Determinant(M)`, `Rank(M)` | `det(M)`, `rank(M)` |
| `Transpose(M)`, `M^-1` | `transpose(M)`, `M^-1` |
| `ChangeRing(M, R)` | `change_base_ring(R, M)` |
| `CharacteristicPolynomial(M)`, `MinimalPolynomial(M)` | `charpoly(M)`, `minpoly(M)` |
| `Eigenvalues(M)` | `eigenvalues(M)` |
| `Kernel(M)`, `NullspaceMatrix(M)` | `kernel(M)` |
| `HermiteForm(M)`, `SmithForm(M)` | `hnf(M)`, `snf(M)` |
| `Vector(K, [1, 2])` | `K.([1, 2])` or `matrix(K, [1 2])` |
