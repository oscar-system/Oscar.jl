# Notes for Macaulay2 users

OSCAR and Macaulay2 overlap in commutative algebra and algebraic
geometry; OSCAR uses Singular for much of this.
This page describes differences between Macaulay2 and OSCAR,
and lists common Macaulay2 commands together with their OSCAR
counterparts.

!!! note "Help wanted"
    This page is a start. Please tell us what is missing, see
    [Notes for users of other computer algebra systems](@ref).

## Differences in syntax

- Function arguments are always parenthesized.
  Macaulay2 allows `factor 120`, `gens R` and `res I`;
  in Julia one writes `factor(120)`, `gens(R)` and `free_resolution(I)`.

- Lists are written `[1, 2, 3]` instead of `{1, 2, 3}`,
  and **they are indexed from 1**: Macaulay2's `L#0` is `L[1]`,
  and `#L` is `length(L)`.

- Control structures need an `end`:
  `for i from 1 to 3 do ...` becomes `for i in 1:3 ... end`,
  and `if c then A else B` becomes `if c A else B end`.

- Comments start with `#` instead of `--`,
  and `help foo` is `?foo`.

- A polynomial ring is created together with its variables,
  see [Many constructors return more than one object](@ref):
  `R = QQ[x,y]` becomes
  `R, (x, y) = polynomial_ring(QQ, [:x, :y])`,
  and `S = R/I` becomes `S, _ = quo(R, I)`.

## Differences in semantics

- **`/` and `//` mean the opposite of what you expect.**
  In Macaulay2, `7/2` is the rational number ``\frac{7}{2}`` and `7//2`
  is the integer `3`.
  In Julia, `7/2` is the floating point number `3.5`, `7//2` is the
  rational number ``\frac{7}{2}``, and integer division is `div(7, 2)`.
  Beware that `2^100` evaluates to `0`, see
  [Integers and rational numbers](@ref other_integers).

- **There is no current ring.**
  In Macaulay2, the symbols `x` and `y` refer to the ring created most
  recently, and `use R` makes them refer to `R` again.
  In OSCAR, every element knows its parent, that is, the algebraic
  structure it lives in (see [Every object has a parent](@ref)),
  and the variables returned by `polynomial_ring` keep referring to
  their ring, whatever else is created afterwards.
  ```jldoctest
  julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

  julia> S, (u, v) = polynomial_ring(QQ, [:u, :v]);

  julia> parent(x)
  Multivariate polynomial ring in 2 variables x, y
    over rational field

  julia> x + u
  ERROR: parents do not match
  [...]
  ```

- **You rarely need a Gröbner basis yourself.**
  Both systems compute Gröbner bases on demand and cache them:
  `dim(I)`, `normal_form(f, I)` and `f in I` take the ideal,
  as `dim I` and `f % I` do in Macaulay2.
  Call `groebner_basis(I)` only if you want the basis itself,
  where you would write `gens gb I`.
  If the ideal is there to define a variety, ask the variety instead:
  `dim(variety(I))`.
  Both systems use the degree reverse lexicographical ordering by
  default; another one is passed to `groebner_basis` as a keyword
  argument, as in `groebner_basis(I; ordering = lex(R))`.

- **Gradings are explicit.**
  A polynomial ring in OSCAR is not graded unless you say so.
  Functions that need a grading, such as `hilbert_series` and
  `betti_table`, work over `grade(R)[1]` rather than over `R`.

- **A quotient ring does not take over the variables.**
  In Macaulay2, `S = R/I` makes `x` and `y` refer to elements of `S`
  from then on:
  ```m2
  i1 : R = QQ[x,y]; S = R/ideal(x^2-y);

  i2 : x^2 == y

  o2 = true
  ```
  In OSCAR, `quo(R, I)` returns the quotient ring together with the
  projection map, `x` stays an element of `R`, and the map moves
  elements into the quotient.
  ```jldoctest
  julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

  julia> Q, p = quo(R, ideal(R, [x^2 - y]));

  julia> x^2 == y
  false

  julia> p(x)^2 == p(y)
  true
  ```

## Common Macaulay2 commands and their OSCAR counterparts

### Rings and polynomials

| Macaulay2 | OSCAR |
|:----------|:------|
| `R = QQ[x,y]` | `R, (x, y) = polynomial_ring(QQ, [:x, :y])` |
| `R = ZZ/32003[x,y]` | `R, (x, y) = polynomial_ring(GF(32003), [:x, :y])` |
| `QQ[x,y,MonomialOrder=>Lex]` | `groebner_basis(I; ordering = lex(R))` |
| `coefficientRing R`, `numgens R`, `gens R` | `coefficient_ring(R)`, `ngens(R)`, `gens(R)` |
| `ring f` | `parent(f)` |
| `degree f` | `total_degree(f)` |
| `leadTerm f`, `leadCoefficient f`, `leadMonomial f` | `leading_term(f)`, `leading_coefficient(f)`, `leading_monomial(f)` |
| `terms f`, `support f`, `exponents f` | `terms(f)`, `vars(f)`, `exponents(f)` |
| `size f` | `length(collect(terms(f)))` |
| `diff(x, f)` | `derivative(f, x)` |
| `substitute(f, {x=>1})` | `evaluate(f, [x], [1])` |
| `factor f`, `gcd(f, g)` | `factor(f)`, `gcd(f, g)` |
| `resultant(f, g, x)`, `discriminant(f, x)` | `resultant(f, g, i)`, `discriminant(f, i)` for the `i`-th variable |
| `isHomogeneous f` | `is_homogeneous(f)` |
| `homogenize(f, z)` | `h = homogenizer(R, :z); h(f)` |
| `S = R/I` | `S, _ = quo(R, I)` |

### Ideals and modules

| Macaulay2 | OSCAR |
|:----------|:------|
| `I = ideal(f, g)` | `I = ideal(R, [f, g])` |
| `numgens I`, `I_0` | `ngens(I)`, `I[1]` |
| `gens gb I` | `groebner_basis(I)` |
| `dim I`, `codim I`, `degree I` | `dim(I)`, `codim(I)`, `degree(I)` |
| `f % (gb I)` | `normal_form(f, I)` |
| `isSubset(ideal f, I)` | `f in I` |
| `radical I` | `radical(I)` |
| `primaryDecomposition I`, `minimalPrimes I` | `primary_decomposition(I)`, `minimal_primes(I)` |
| `eliminate(I, x)` | `eliminate(I, [x])` |
| `I : J`, `saturate(I, J)` | `quotient(I, J)`, `saturation(I, J)` |
| `intersect(I, J)`, `I + J`, `I * J` | `intersect(I, J)`, `I + J`, `I * J` |
| `syz M` | `syzygy_generators(gens(I))` |
| `res I` | `free_resolution(I)` |
| `betti res I` | `betti_table(free_resolution(quo(grade(R)[1], I)[1]))` |
| `hilbertSeries(R/I)` | `hilbert_series(quo(grade(R)[1], I)[1])` |
| `jacobian f` | `jacobian_matrix(f)` |

### Matrices

| Macaulay2 | OSCAR |
|:----------|:------|
| `matrix{{1,2},{3,4}}` | `matrix(ZZ, [1 2; 3 4])` |
| `id_(ZZ^2)` | `identity_matrix(ZZ, 2)` |
| `M_(0,0)` | `M[1, 1]` |
| `numrows M`, `numcols M` | `nrows(M)`, `ncols(M)` |
| `det M`, `rank M`, `transpose M` | `det(M)`, `rank(M)`, `transpose(M)` |
| `inverse M` | `inv(M)` |
| `ker M` | `kernel(M)` |

### Lists and numbers

| Macaulay2 | OSCAR |
|:----------|:------|
| `{1, 2, 3}`, `#L`, `L#0` | `[1, 2, 3]`, `length(L)`, `L[1]` |
| `append(L, x)`, `join(L, M)` | `push!(L, x)`, `vcat(L, M)` |
| `apply(L, f)`, `select(L, f)` | `map(f, L)`, `filter(f, L)` |
| `position(L, f)`, `member(x, L)` | `findfirst(f, L)`, `x in L` |
| `sum L`, `product L`, `max L` | `sum(L)`, `prod(L)`, `maximum(L)` |
| `sort L`, `reverse L` | `sort(L)`, `reverse(L)` |
| `set L`, `toList S` | `Set(L)`, `collect(S)` |
| `new HashTable from {a=>1}`, `H#a` | `Dict(:a => 1)`, `H[:a]` |
| `ZZ`, `QQ`, `GF 9`, `ZZ/7` | `ZZ`, `QQ`, `GF(9)`, `residue_ring(ZZ, 7)` |
| `isPrime n`, `nextPrime n`, `factor n` | `is_prime(n)`, `next_prime(n)`, `factor(n)` |
| `gcd(a, b)`, `lcm(a, b)`, `binomial(n, k)` | `gcd(a, b)`, `lcm(a, b)`, `binomial(n, k)` |
