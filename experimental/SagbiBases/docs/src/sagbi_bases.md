# [Subalgbera (SAGBI) Bases](@id sagbi_bases)

!!! note
    If this code is moved from `experimental` to `src`, I suggest that this
    documentation page be placed under `Commutative Algebra/Grobner Bases/` in the
    documentation tree.

A SAGBI (*Subalgebra Basis Analogue to Gr&ouml;bner Bases*) basis is, as the
acronym suggests, a basis of a subalgebra which is analogous to a Gr&ouml;bner 
basis of an ideal. 
See: [Robbiano and Sweedler 1990](https://doi.org/10.1007/BFb0085537)

Let $R$ be a polynomial ring with a specified 
[monomial order](@ref monomial_orderings). Given a subalgebra $A$ of $R$, a 
SAGBI basis for $A$ is a set of generators $B=\{b_i\}_{i\in I}$ of $A$ whose
initial terms generate the initial algebra of $A$. This means that for any 
$f\in A$, there are non-negative integers $\{\alpha_i\}_{i\in I}$ such that

$ \text{LT}(f) = \prod_{i\in I} \text{LT}(b_i)^{\alpha_i}.$


## Subduction of an Element

lorem ipsum (todo)

```@docs
subduct(
    f::MPolyRingElem, generators::Vector{<:MPolyRingElem}
)-> MPolyRingElem
```

```@docs
monoid_representation(
    m::MPolyRingElem, 
    B::SAGBICandidate,
)-> Vector{Int}
```

## Checking the SAGBI Criterion

lorem ipsum (todo)

```@docs
is_sagbi(
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
) -> Bool
```

*t&ecirc;te-&agrave;-t&ecirc;te*

```@docs
function tete_a_tetes(
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(generators[1]))
)
```

## Computing a SAGBI basis for a Subalgbera

lorem ipsum (todo)

```@docs
is_sagbi(
    generators::Vector{<:MPolyRingElem};
    ordering::MonomialOrdering = default_ordering(parent(B[1]))
) -> Bool
```

```@docs
compute_sagbi_degree!(B::SAGBICandidate) -> SAGBICandidate
```
