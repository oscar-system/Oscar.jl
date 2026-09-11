# Notes for PARI/GP users

OSCAR and PARI/GP overlap in number theory.
Much of OSCAR's number theory lives in the packages
[Nemo](https://github.com/nemocas/Nemo.jl) and
[Hecke](https://github.com/thofma/Hecke.jl), whose functions are
available in OSCAR.
This page describes differences between GP and OSCAR,
and lists common GP functions together with their OSCAR counterparts.

!!! note "Help wanted"
    This page is a start. Please tell us what is missing, see
    [Notes for users of other computer algebra systems](@ref).

## Differences in syntax

- Defining a function looks the same in both languages: `f(x) = x^2`.
  Anonymous functions are `x -> x^2` in Julia,
  and local variables need no `my`.

- Control structures need an `end` and take no parentheses:
  `for(i = 1, 3, ...)` becomes `for i in 1:3 ... end`,
  and `if(c, A, B)` becomes `if c A else B end`.

- Comments start with `#` instead of `\\`,
  and multi-line comments are `#= ... =#` instead of `/* ... */`.

- `#v` is `length(v)`, and `M~` is `transpose(M)`.
  Vectors and matrices are indexed from 1 in both systems.
  GP's `apply`, `select` and `vecsort` are Julia's `map`, `filter` and
  `sort`; `concat` is `vcat`, and `vecsum`, `vecprod`, `vecmax` are
  `sum`, `prod`, `maximum`.

- `?name` shows the documentation in both systems.
  `##` is `@time`, and `\r file` is `include("file.jl")`.

- A trailing `;` suppresses the output of a value in both systems,
  see [Semicolons and output](@ref).

## Differences in semantics

- **Integer literals are machine integers.**
  In GP, integers and rational numbers are exact and of arbitrary size,
  and `7/2` is the rational number `7/2`.
  In Julia, `2^100` evaluates to `0` and `7/2` to the floating point
  number `3.5`; write `ZZ(2)^100` and `QQ(7, 2)`,
  see [Integers and rational numbers](@ref other_integers).
  GP's `7\2` and `7%2` are `div(7, 2)` and `mod(7, 2)`.

- **The ring comes first.**
  In GP, `Mod(3, 7)` and `x^2 + 1` are just values, and the ring they
  live in is implied.
  In OSCAR, one creates the ring first and takes elements from it:
  `Zn, _ = residue_ring(ZZ, 7)` and then `Zn(3)`,
  or `R, x = polynomial_ring(QQ, :x)` and then `x^2 + 1`,
  see [Many constructors return more than one object](@ref).
  Every element knows its ring, and `parent(a)` returns it,
  see [Every object has a parent](@ref).
  OSCAR elements have a Julia type as well, `typeof(a)`, but the type
  does not determine the ring: elements of `ZZ/7` and of `ZZ/11` have
  the same type and different parents.
  `lift` maps an element of `ZZ/n` back to the integers.

- **Real and complex numbers come from a field with a precision.**
  GP computes with a global real precision that `\p` changes, so
  `sqrt(2)` and `Pi` just work.
  OSCAR's `RealField()` and `ComplexField()` behave in the same way:
  their precision is a global setting, `precision(Balls)`, which
  `set_precision!(Balls, 128)` changes.
  `ArbField(prec)` and `AcbField(prec)` fix the precision in the field
  instead.
  All four compute with balls that carry rigorous error bounds,
  `sqrt(RealField()(2))` shows one.
  Exact alternatives are often preferable, such as `quadratic_field(2)`
  or `algebraic_closure(QQ)`.

- **Number fields are objects, not `init` structures.**
  Instead of `nfinit` and `bnfinit`, which precompute data that later
  functions expect, OSCAR has the field itself, and each invariant is a
  function of it, or of its maximal order; what is needed gets computed
  and cached along the way.
  The GP session
  ```gp
  ? K = bnfinit(x^2 + 5);
  ? K.disc
  -20
  ? K.zk
  [1, x]
  ? K.clgp
  [2, [2], [[2, 1; 0, 1]]]
  ```
  corresponds to
  ```jldoctest
  julia> R, x = polynomial_ring(QQ, :x);

  julia> K, a = number_field(x^2 + 5, "a");

  julia> discriminant(K)
  -20

  julia> OK = maximal_order(K);

  julia> basis(OK)
  2-element Vector{AbsSimpleNumFieldOrderElem}:
   1
   a

  julia> class_group(OK)
  (Z/2, Class group map of set of ideals of OK)
  ```

## Common GP functions and their OSCAR counterparts

### Integers and rational numbers

| GP | OSCAR |
|:---|:------|
| `factor(n)`, `isprime(n)`, `nextprime(n)` | `factor(n)`, `is_prime(n)`, `next_prime(n)` |
| `divisors(n)`, `numdiv(n)`, `sigma(n)` | `divisors(n)`, `divisor_sigma(n, 0)`, `divisor_sigma(n, 1)` |
| `eulerphi(n)`, `moebius(n)` | `euler_phi(n)`, `moebius_mu(n)` |
| `gcd(a, b)`, `lcm(a, b)`, `bezout(a, b)` | `gcd(a, b)`, `lcm(a, b)`, `gcdx(a, b)` |
| `binomial(n, k)`, `fibonacci(n)` | `binomial(n, k)`, `fibonacci(n)` |
| `issquare(n)`, `sqrtint(n)`, `sqrtnint(n, k)` | `is_square(n)`, `isqrt(n)`, `iroot(n, k)` |
| `kronecker(a, n)` | `kronecker_symbol(a, n)`, `jacobi_symbol(a, n)` |
| `numerator(q)`, `denominator(q)` | `numerator(q)`, `denominator(q)` |
| `floor(q)`, `ceil(q)`, `round(q)`, `truncate(q)` | `floor(ZZRingElem, q)`, `ceil(ZZRingElem, q)`, `round(ZZRingElem, q)`, `trunc(ZZRingElem, q)` |
| `Mod(a, n)`, `lift(x)` | `Zn, _ = residue_ring(ZZ, n); Zn(a)`, `lift(x)` |
| `Mod(a, n)^-1` | `inv(Zn(a))` |
| `chinese(Mod(a, m), Mod(b, n))` | `crt([a, b], [m, n])` |
| `ffgen(ffinit(p, n))` | `gen(GF(p, n))` |
| `O(p^k)` | `padic_field(p, precision = k)` |

### Polynomials

| GP | OSCAR |
|:---|:------|
| `x` | `R, x = polynomial_ring(QQ, :x)` |
| `poldegree(f)`, `polcoef(f, i)` | `degree(f)`, `coeff(f, i)` |
| `Vec(f)` | `collect(coefficients(f))` |
| `deriv(f)`, `subst(f, x, 2)` | `derivative(f)`, `evaluate(f, 2)` |
| `factor(f)`, `polisirreducible(f)` | `factor(f)`, `is_irreducible(f)` |
| `polroots(f)` | `roots(AcbField(64), f)` |
| `polrootsreal(f)` | `roots(ArbField(64), f)` |
| `polresultant(f, g)`, `poldisc(f)` | `resultant(f, g)`, `discriminant(f)` |
| `polgalois(f)`, `nfsplitting(f)` | `galois_group(number_field(f)[1])`, `splitting_field(f)` |

### Number fields

| GP | OSCAR |
|:---|:------|
| `nfinit(f)`, `bnfinit(f)` | `K, a = number_field(f, "a")` |
| `nf.pol`, `nf.disc`, `poldegree(nf.pol)` | `defining_polynomial(K)`, `discriminant(K)`, `degree(K)` |
| `nf.zk`, `nfbasis(f)` | `basis(maximal_order(K))` |
| `bnf.clgp`, `bnf.no`, `bnf.cyc` | `class_group(maximal_order(K))`, `class_number(K)` |
| `bnf.fu`, `bnfunit(bnf)` | `unit_group(maximal_order(K))` |
| `idealprimedec(nf, p)` | `prime_decomposition(maximal_order(K), p)` |
| `idealfactor(nf, x)`, `idealnorm(nf, I)` | `factor(I)`, `norm(I)` |
| `nfeltnorm(nf, a)`, `nfelttrace(nf, a)` | `norm(a)`, `tr(a)` |
| `minpoly(a)`, `charpoly(a)` | `minpoly(a)`, `charpoly(a)` |
| `quadgen(d)`, `qfbclassno(d)` | `quadratic_field(d)`, `class_number(quadratic_field(d)[1])` |
| `polcyclo(n)` | `cyclotomic_field(n)` |

### Elliptic curves

| GP | OSCAR |
|:---|:------|
| `ellinit([a4, a6])` | `E = elliptic_curve(QQ, [a4, a6])` |
| `e.disc`, `e.j` | `discriminant(E)`, `j_invariant(E)` |
| `ellglobalred(e)[1]` | `conductor(E)` |
| `elltors(e)` | `torsion_structure(E)`, `torsion_points(E)` |
| `elladd(e, P, Q)`, `ellneg(e, P)` | `P + Q`, `-P` |
| `ellorder(e, P)` | `order(P)` |
| `ellisoncurve(e, P)` | `P in E` |

### Matrices and vectors

| GP | OSCAR |
|:---|:------|
| `[1, 2; 3, 4]` | `matrix(ZZ, [1 2; 3 4])` |
| `matid(n)`, `matsize(M)` | `identity_matrix(ZZ, n)`, `size(M)` |
| `matdet(M)`, `matrank(M)`, `M~` | `det(M)`, `rank(M)`, `transpose(M)` |
| `1/M` | `inv(M)` |
| `matsolve(M, b)` | `solve(M, b; side = :right)` |
| `matker(M)` | `kernel(M; side = :right)` |
| `mathnf(M)`, `matsnf(M)` | `hnf(M)`, `snf(M)` |
| `charpoly(M)`, `mateigen(M)` | `charpoly(M)`, `eigenvalues(M)` |
| `[1, 2, 3]`, `#v`, `v[1]` | `[1, 2, 3]`, `length(v)`, `v[1]` |
| `concat(v, w)`, `vecsort(v)` | `vcat(v, w)`, `sort(v)` |
| `apply(f, v)`, `select(f, v)` | `map(f, v)`, `filter(f, v)` |
| `vecsum(v)`, `vecprod(v)`, `vecmax(v)` | `sum(v)`, `prod(v)`, `maximum(v)` |
| `Set(v)` | `Set(v)` |
