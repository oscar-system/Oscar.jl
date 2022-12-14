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
Pages = ["groebner_bases_integers.md"]
```

# Gröbner/Standard Bases Over $\mathbb Z$

In this section, we consider a polynomial ring
$\mathbb Z[x] = \mathbb Z[x_1, \dots, x_n]$ over the integers. As in the previous section
on Gröbner/standard bases over fields, let $>$ be a monomial ordering on $\text{Mon}_n(x)$.
With respect to this ordering, the localization $\mathbb Z[x]_>$ and, given a nonzero element
$f \in \mathbb Z[x]_>$, the notions *leading term*, *leading monomial*, *leading exponent*,
*leading coefficient*, and *tail*  of $f$ are defined as before.

!!! note
    Over $\mathbb Z$, the basic idea of multivariate polynomial division with remainder in OSCAR is as follows:
    If $\text{LT}(\widetilde{g}) = ax^\alpha$ is the leading term of the intermediate dividend, $f_i$
    is *some* divisor whose leading monomial equals $x^\alpha$, say
    $\text{LT}(f_i) = bx^\alpha$, and $r$ is the remainder of $a$ on division by $b$ in
    $\mathbb Z$, then $ax^\alpha$ is replaced by $rx^\alpha$.

##### Examples

```jldoctest
julia> R, (x, y) = PolynomialRing(ZZ, ["x", "y"]);

julia> reduce(3*x, [2*x])
x

julia> reduce(6*x, [5*x, 2*x])
0
```

The notion of *leading ideals*  as formulated in the previous section and the definitions of
standard bases (Gröbner bases) carry over: A *standard basis* for an ideal $I\subset K[x]_>$
with respect to $>$ is a finite subset $G$ of $I$ such that $\text{L}_>(G) = \text{L}_>(I)$ (a
standard basis with respect to a global monomial ordering is also called a *Gröbner basis*).

There is, however, a sublety: Over a field, the defining condition of a standard basis as stated
above is equivalent to saying that the $\text{LT}_>(g)$, $g\in G\setminus\{0\}$ generate
$\text{L}_>(I)$. Over $\mathbb Z$, the latter condition implies the former one, but
not vice versa. Consequently, over $\mathbb Z$, a finite subset $G$ of $I$ satisfying the
latter condition is called a *strong standard basis* for $I$ (with respect to $>$).

We refer to the  textbook [AL94](@cite) for more on this.

!!! note
    Over $\mathbb Z$, the standard bases returned by OSCAR are strong in the sense above.

##### Examples

```jldoctest
julia> R, (x,y) = PolynomialRing(ZZ, ["x","y"])
(Multivariate Polynomial Ring in x, y over Integer Ring, fmpz_mpoly[x, y])

julia> I = ideal(R, [3*x^2*y+7*y, 4*x*y^2-5*x])
ideal(3*x^2*y + 7*y, 4*x*y^2 - 5*x)

julia> G = groebner_basis(I, ordering = lex(R))
Gröbner basis with elements
1 -> 28*y^3 - 35*y
2 -> 4*x*y^2 - 5*x
3 -> 15*x^2 + 28*y^2
4 -> 3*x^2*y + 7*y
5 -> x^2*y^2 - 5*x^2 - 7*y^2
with respect to the ordering
lex([x, y])
```

