# Notes for SageMath users

SageMath and OSCAR share several of their components,
for example both use GAP for group theory and Singular for commutative
algebra.
There is no interface between the two systems, though,
so this page is about translating what you know from SageMath to OSCAR.
SageMath is built on Python, OSCAR is built on Julia,
so most of the differences are differences between these two languages.

!!! note "Help wanted"
    This page is a start. Please tell us what is missing, see
    [Notes for users of other computer algebra systems](@ref).

## Differences in syntax

- Blocks are not defined by indentation but end with the keyword `end`:
  `for x in L ... end`, `if c ... elseif d ... else ... end`,
  `function f(x) ... end`, `while c ... end`.
  There is no colon after the condition.

- `True`, `False`, `None` are `true`, `false`, `nothing`,
  and `and`, `or`, `not` are `&&`, `||`, `!`.

- One-line functions can be written as `f(x) = x^2`,
  and `lambda x: x^2` is `x -> x^2`.

- Indexing starts at 1, not at 0, and `L[end]` is the last entry.
  `range(1, 11)` is `1:10`; both endpoints are included.
  Slicing `L[1:3]` (Python) is `L[2:3]` (Julia).

- Functions are not methods of objects:
  `G.order()` is `order(G)`, `f.factor()` is `factor(f)`,
  `I.groebner_basis()` is `groebner_basis(I)`.
  Instead of `G.<TAB>` for exploring what one can do with `G`,
  use `methodswith(typeof(G))`,
  see [Names of functions](@ref).

- Exponentiation is `^` (as in SageMath, but unlike in plain Python).

- Comprehensions look alike: `[f(x) for x in L if c(x)]`.

- `print(x)` is `println(x)`.
  `str(x)` is `string(x)`, `type(x)` is `typeof(x)`,
  and `isinstance(x, T)` is `x isa T`.

- `?f` works in both systems.
  `%time f(x)` is `@time f(x)`, and `load("file.sage")` is
  `include("file.jl")`.

- SageMath's preparser syntax `R.<x, y> = PolynomialRing(QQ)` and
  `K.<a> = NumberField(f)` do not exist.
  The OSCAR constructors return the ring together with its generators:
  ```jldoctest
  julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
  (Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

  julia> P, t = polynomial_ring(QQ, :t);

  julia> K, a = number_field(t^2 + 1, "a")
  (Number field of degree 2 over QQ, a)
  ```

## Differences in semantics

- **Integer literals are machine integers.**
  The SageMath preparser turns `2^100` into an exact integer and `3/4`
  into a rational number.
  In Julia, `2^100` evaluates to `0` and `3/4` to `0.75`,
  see [Integers and rational numbers](@ref other_integers).

- **There is no automatic coercion between parents.**
  SageMath's coercion system finds a common parent for the operands of
  an arithmetic operation.
  OSCAR requires that the parents coincide, except that integers and
  rational numbers can be mixed with elements of any ring;
  see [Every object has a parent](@ref).

- **There is no symbolic ring.**
  `sqrt(2)`, `pi`, and `var('u')` have no counterparts;
  work in an algebraic structure instead, such as a number field
  (`quadratic_field(2)`) or the algebraic closure `algebraic_closure(QQ)`.

- **Constructors return tuples.**
  `Integers(5)` is `residue_ring(ZZ, 5)`, which returns the ring
  together with the projection from `ZZ`.
  Similarly, `G.center()` is `center(G)`, which returns the center together
  with its embedding into `G`, and `G.quotient(N)` is `quo(G, N)`,
  which returns the quotient together with the projection,
  see [Many constructors return more than one object](@ref).

- **The order of dihedral groups.**
  `DihedralGroup(n)` in SageMath is the dihedral group of order ``2n``,
  OSCAR's `dihedral_group(n)` is the dihedral group of order ``n``,
  following GAP.

- **Matrices are not Python lists.**
  SageMath's `matrix(QQ, [[1, 2], [3, 4]])` is `matrix(QQ, [1 2; 3 4])`
  in OSCAR; the entry `M[0, 1]` is `M[1, 2]`,
  and the first row `M[0]` is `M[1, :]`.
  `M.kernel()` is `kernel(M)`, both compute the left kernel by default.

- **Accessing GAP and Singular.**
  Prefer the OSCAR functions (`symmetric_group(4)`, `groebner_basis(I)`)
  whenever they exist, they take care of the conversions.
  If you need something that is only available in GAP,
  `libgap.SymmetricGroup(4)` is `GAP.Globals.SymmetricGroup(4)`,
  and `G.gap()` is `GapObj(G)`, see [Notes for GAP users](@ref).
  The Singular kernel is accessed via the Julia package
  [Singular.jl](https://github.com/oscar-system/Singular.jl),
  available as `Oscar.Singular`.

## Common SageMath functions and their OSCAR counterparts

### Lists

| SageMath | OSCAR |
|:---------|:------|
| `len(L)` | `length(L)` |
| `L.append(x)` | `push!(L, x)` |
| `L + M` | `vcat(L, M)` |
| `L[::-1]`, `sorted(L)` | `reverse(L)`, `sort(L)` |
| `L.index(x)` | `findfirst(==(x), L)` |
| `max(L)`, `min(L)`, `sum(L)`, `prod(L)` | `maximum(L)`, `minimum(L)`, `sum(L)`, `prod(L)` |
| `set(L)` | `Set(L)` or `unique(L)` |
| `A \| B`, `A & B`, `A - B` | `union(A, B)`, `intersect(A, B)`, `setdiff(A, B)` |
| `Combinations(L, k)`, `Partitions(n)` | `combinations(L, k)`, `partitions(n)` |
| `cartesian_product([A, B])` | `Iterators.product(A, B)` |
| `copy(x)`, `deepcopy(x)` | `copy(x)`, `deepcopy(x)` |

### Integers and rational numbers

| SageMath | OSCAR |
|:---------|:------|
| `ZZ`, `QQ` | `ZZ`, `QQ` |
| `factor(n)`, `n.factor()` | `factor(n)` |
| `is_prime(n)`, `next_prime(n)` | `is_prime(n)`, `next_prime(n)` |
| `divisors(n)`, `euler_phi(n)` | `divisors(n)`, `euler_phi(n)` |
| `gcd(a, b)`, `lcm(a, b)` | `gcd(a, b)`, `lcm(a, b)` |
| `binomial(n, k)`, `factorial(n)` | `binomial(n, k)`, `factorial(ZZ(n))` |
| `power_mod(a, e, m)` | `powermod(a, e, m)` |
| `CRT([r1, r2], [m1, m2])` | `crt([r1, r2], [m1, m2])` |
| `q.numerator()`, `q.denominator()` | `numerator(q)`, `denominator(q)` |
| `floor(q)`, `ceil(q)`, `round(q)` | `floor(q)`, `ceil(q)`, `round(q)` |
| `isqrt(n)`, `n.is_square()` | `isqrt(n)`, `is_square(n)` |
| `Integers(n)`, `Zmod(n)` | `residue_ring(ZZ, n)` |
| `GF(q)`, `GF(q, 'a')` | `GF(q)` |
| `Mod(a, n)` | `residue_ring(ZZ, n)[1](a)` |

### Groups

| SageMath | OSCAR |
|:---------|:------|
| `SymmetricGroup(n)`, `AlternatingGroup(n)` | `symmetric_group(n)`, `alternating_group(n)` |
| `CyclicPermutationGroup(n)`, `DihedralGroup(n)` | `cyclic_group(PermGroup, n)`, `dihedral_group(PermGroup, 2n)` |
| `AbelianGroup([2, 4])` | `abelian_group([2, 4])` |
| `PermutationGroup([[(1,2,3,4)], [(1,2)]])` | `permutation_group(4, [cperm([1, 2, 3, 4]), cperm([1, 2])])` |
| `G((1,2,3))` | `cperm(G, [1, 2, 3])` |
| `Permutation([2, 3, 1])` | `perm([2, 3, 1])` |
| `G.order()`, `g.order()` | `order(G)`, `order(g)` |
| `G.gens()`, `G.gen(0)` | `gens(G)`, `G[1]` |
| `G.identity()` | `one(G)` |
| `G.list()` | `collect(G)` |
| `G.random_element()` | `rand(G)` |
| `G.subgroup([a, b])` | `sub(G, [a, b])` |
| `G.center()`, `G.commutator()` | `center(G)`, `derived_subgroup(G)` |
| `G.centralizer(x)`, `G.normalizer(H)` | `centralizer(G, x)`, `normalizer(G, H)` |
| `G.sylow_subgroup(p)` | `sylow_subgroup(G, p)` |
| `G.quotient(N)` | `quo(G, N)` |
| `H.is_subgroup(G)`, `H.is_normal(G)` | `is_subgroup(H, G)`, `is_normal_subgroup(H, G)` |
| `G.is_abelian()`, `G.is_solvable()`, `G.is_simple()` | `is_abelian(G)`, `is_solvable(G)`, `is_simple(G)` |
| `G.conjugacy_classes()` | `conjugacy_classes(G)` |
| `G.character_table()` | `character_table(G)` |
| `G.structure_description()`, `G.group_id()` | `describe(G)`, `small_group_identification(G)` |
| `G.is_isomorphic(H)` | `is_isomorphic(G, H)` |
| `G.orbit(x)`, `G.stabilizer(x)` | `orbit(G, x)`, `stabilizer(G, x)` |
| `G.is_transitive()` | `is_transitive(G)` |
| `g.sign()`, `g.cycle_type()` | `sign(g)`, `cycle_structure(g)` |

### Polynomials and ideals

| SageMath | OSCAR |
|:---------|:------|
| `R.<x, y> = PolynomialRing(QQ)` | `R, (x, y) = polynomial_ring(QQ, [:x, :y])` |
| `S.<t> = QQ[]` | `S, t = polynomial_ring(QQ, :t)` |
| `f.degree()`, `f.total_degree()` | `degree(f)`, `total_degree(f)` |
| `f.derivative(x)` | `derivative(f, x)` |
| `f(1, 2)` | `f(1, 2)` or `evaluate(f, [1, 2])` |
| `f.coefficients()`, `f.monomials()` | `coefficients(f)`, `monomials(f)` |
| `f.lc()`, `f.lm()`, `f.lt()` | `leading_coefficient(f)`, `leading_monomial(f)`, `leading_term(f)` |
| `f.factor()`, `f.is_irreducible()` | `factor(f)`, `is_irreducible(f)` |
| `f.roots()` | `roots(f)` |
| `f.discriminant()`, `f.resultant(g)` | `discriminant(f)`, `resultant(f, g)` |
| `R.ideal([f, g])` | `ideal(R, [f, g])` |
| `I.groebner_basis()` | `groebner_basis(I)` |
| `I.dimension()`, `I.radical()` | `dim(I)`, `radical(I)` |
| `I.primary_decomposition()` | `primary_decomposition(I)` |
| `I.elimination_ideal(x)` | `eliminate(I, [x])` |
| `I.reduce(f)`, `f in I` | `normal_form(f, I)`, `f in I` |
| `I.variety()` | `rational_solutions(I)` |
| `R.quotient(I)` | `quo(R, I)` |

### Number fields

| SageMath | OSCAR |
|:---------|:------|
| `K.<a> = NumberField(f)` | `K, a = number_field(f, "a")` |
| `CyclotomicField(n)`, `QuadraticField(d)` | `cyclotomic_field(n)`, `quadratic_field(d)` |
| `K.ring_of_integers()`, `K.maximal_order()` | `maximal_order(K)` |
| `K.discriminant()`, `K.degree()` | `discriminant(K)`, `degree(K)` |
| `a.minpoly()`, `a.norm()`, `a.trace()` | `minpoly(a)`, `norm(a)`, `tr(a)` |
| `K.class_group()`, `K.class_number()` | `class_group(K)`, `class_number(K)` |
| `K.unit_group()` | `unit_group(maximal_order(K))` |
| `K.galois_group()` | `galois_group(K)` |
| `K.ideal(2)` | `ideal(maximal_order(K), 2)` |
| `I.factor()`, `I.is_prime()` | `factor(I)`, `is_prime(I)` |

### Matrices

| SageMath | OSCAR |
|:---------|:------|
| `matrix(QQ, [[1, 2], [3, 4]])` | `matrix(QQ, [1 2; 3 4])` |
| `identity_matrix(QQ, n)`, `zero_matrix(QQ, m, n)` | `identity_matrix(QQ, n)`, `zero_matrix(QQ, m, n)` |
| `vector(QQ, [1, 2])` | `QQ.([1, 2])` or `matrix(QQ, [1 2])` |
| `M.nrows()`, `M.ncols()` | `nrows(M)`, `ncols(M)` |
| `M[i, j]`, `M[i]` | `M[i + 1, j + 1]`, `M[i + 1, :]` |
| `M.det()`, `M.trace()`, `M.rank()` | `det(M)`, `tr(M)`, `rank(M)` |
| `M.transpose()`, `M.inverse()` | `transpose(M)`, `inv(M)` |
| `M.change_ring(R)` | `change_base_ring(R, M)` |
| `M.charpoly()`, `M.minpoly()` | `charpoly(M)`, `minpoly(M)` |
| `M.eigenvalues()` | `eigenvalues(M)` |
| `M.kernel()`, `M.right_kernel()` | `kernel(M)`, `kernel(M; side = :right)` |
| `M.rref()` | `rref(M)` |
| `M.hermite_form()`, `M.smith_form()` | `hnf(M)`, `snf(M)` |

### Polyhedral geometry and graphs

| SageMath | OSCAR |
|:---------|:------|
| `Polyhedron(vertices = [[0, 0], [1, 0], [0, 1]])` | `convex_hull([0 0; 1 0; 0 1])` |
| `polytopes.cube()` | `cube(3)` |
| `P.vertices()`, `P.facets()` | `vertices(P)`, `facets(P)` |
| `P.f_vector()`, `P.dim()`, `P.volume()` | `f_vector(P)`, `dim(P)`, `volume(P)` |
| `Graph([(1, 2), (2, 3)])` | `graph_from_edges([[1, 2], [2, 3]])` |
| `G.is_connected()`, `G.automorphism_group()` | `is_connected(G)`, `automorphism_group(G)` |
