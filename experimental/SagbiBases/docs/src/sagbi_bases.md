# [Subalgbera (SAGBI) Bases](@id sagbi_bases)

A SAGBI (*Subalgebra Basis Analogue to Gröbner Bases*) basis is, as the
acronym suggests, a basis of a subalgebra analogous to a Gröbner 
basis of an ideal. For the theory of SAGBI basis, see [RS90](@cite).

Let $R$ be a polynomial ring with a specified 
[monomial order](@ref monomial_orderings). Given a subalgebra $A$ of $R$, a 
SAGBI basis for $A$ is a set of generators $B=\{b_i\}_{i\in I}$ of $A$ whose
initial terms generate the initial algebra of $A$. This means that for any 
$f\in A$, there are non-negative integers $\{\alpha_i\}_{i\in I}$ such that

$\text{LM}(f) = \prod_{i\in I} \text{LM}(b_i)^{\alpha_i},$

where $\text{LM}(f)$ is the leading monomial of $f$.

## Subduction of an Element

Suppose $R$ has variables $x=\{x_1,...,x_n\}$, and let $B=\{b_i\}_{i\in I}$ 
be a generating set for a subalgebra $A$. Let 
$\text{LM}(b_i) = x^{\beta_i} := x_1^{\beta_{i,1}}\cdot ...\cdot x_n^{\beta_{i,n}}$, which defines exponent
vectors $\{\beta_i\}_{i\in I}$. Given any other monomial $x^\alpha$, the 
function `monoid_representation` searches for coefficients 
$c_i \in \mathbb{N}^n$ such that $\alpha = \sum c_i \beta_i$. 

```@docs
monoid_representation(
m::MPolyRingElem, generators::Vector{<:MPolyRingElem}
)
```
Given any polynomial $f\in R$, a *subduction step* of $f$ with respect to $B$ 
is computed as follows:
1. Let $\text{LM}(f) = x^\alpha$ and compute $c_i$ such that $\alpha = \sum c_i \beta_i$.
2. Let $f' = f - \prod_{i\in I} b_i^{\beta_i}$ in $\mathbb{N}^n$.

As $f$ and $\prod_{i\in I} b_i^{\beta_i}$ have the same leading monomial, their 
difference $f'$ has a strictly smaller leading monomial. The *subduction* of
$f$ by $B$ is the remainder $f'$ computed by repeating this process until either 
$f'=0$ or there do not exist $c_i\in \mathbb{N}^n$ such that
 $\text{LM}(f) = \sum c_i \beta_i$. In either case, the result  $f'$ is called
 the *subduction of $f$ by $B$*.

```@docs
subduct(
f::MPolyRingElem, generators::Vector{<:MPolyRingElem};
ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
```


## Checking the SAGBI Criterion

A *tête-à-tête* $t\in A$ for the generating set $B=\{b_i\}_{i\in I}$ is a
polynomial of the form

$t = \prod_{i\in I} \left(b_i^{\alpha_i} - b_i^{\beta_i}\right)$

for some non-negative integer exponent vectors $\alpha$ and $\beta$ such that 
$\text{LM}(\prod_{i\in I} b_i^{\alpha_i}) = \text{LM}(\prod_{i\in I} b_i^{\beta_i}).$ 
A generating set $B$ is a SAGBI basis if and only if all its tête-à-têtes 
subduct to 0.

```@docs
tete_a_tetes(
generators::Vector{<:MPolyRingElem};
ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
```

To check if a set of generators is a SAGBI basis, call `is_sagbi`:

```@docs
is_sagbi(
generators::Vector{<:MPolyRingElem};
ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
```
This function computes the subduction of all tête-à-têtes of `generators`, 
verifying if they all return `zero` (returning true) or not (returning false). 



## Computing a SAGBI Basis for a Subalgbera

If a polynomial $f$ subducts to a remainder $r$ by a generating set $B$, then
$f$ will subduct to $0$ by $B\cup\{r\}$. This can be used to complete $B$ to
a SAGBI basis as follows:
1. Compute all the tête-à-têtes of $B$ and their remainders $r$ under subduction.
2. Add these remainders to $B$. This may cause new tête-à-têtes to exist.
3. Repeat steps 1 and 2 until all tête-à-têtes subduct to 0.

However SAGBI bases can be countably infinite, and a finite SAGBI basis need not
exist for a given subalgebra $A$. For this reason, we say the set $B$ has 
*SAGBI degree* $d$ if the SAGBI criterion holds for all monomials in $A$
of degree at most $d$. That is, for all $x^\alpha$ with $|\alpha|\leq d$, 
there exist $\beta_i \in \mathbb{N}^{|B|}$ such that

$x^\alpha = \prod_{i\in I} \text{LM}(b_i)^{\beta_i}.$

A SAGBI basis has SAGBI degree $\infty$. This lets us stop the SAGBI completion
algorithm once we've achieved a SAGBI basis of a given degree bound:

```@docs
sagbi(
    generating_set::Vector{<:MPolyRingElem};
    degree_bound::Integer,
    ordering::MonomialOrdering = default_ordering(parent(generating_set[1]))
)
```
The return type of `sagbi` is `SubalgebraBasis`. The elements of a
SubalgebraBasis $B$ can be accessed with `elements`:

```@docs
elements(
B::SubalgebraBasis
)
```

A SubalgebraBasis $B$ can also be indexed like a vector: `B[i]` returns the
`i`-th element, `length(B)` returns the number of elements, and iterating
over $B$ yields its elements. The polynomial ring containing the elements of
$B$ is available via `base_ring`:

```@docs
base_ring(
B::SubalgebraBasis
)
```

The SAGBI degree of the basis is recorded by `sagbi_degree`. If a
SubalgebraBasis $B$ is a SAGBI basis, then `sagbi_degree(B)` is $-1$:

```@docs
sagbi_degree(
B::SubalgebraBasis
)
```

Given a generating set $B$, you may also compute its SAGBI degree without 
adding new elements to $B$.

```@docs
sagbi_degree(
generators::Vector{<:MPolyRingElem};
ordering=default_ordering(parent(generators[1]))
)
```
