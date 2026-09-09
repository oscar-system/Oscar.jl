```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Blow-ups

The *blow-up* of a smooth variety $X$ along a smooth subvariety $Z$ replaces $Z$ by the
exceptional divisor $E = \mathbb P(N_{Z/X})$, the projectivization of the normal bundle of
$Z$ in $X$. The resulting variety $\widetilde{X}$ comes equipped with a map $\pi\colon \widetilde{X}\to X$
that is an isomorphism away from $E$.

In OSCAR, a blow-up is constructed from an inclusion map $i\colon Z \hookrightarrow X$ (an
`AbstractVarietyMap`). The function `blowup` returns the triple `(Bl, E, j)` where
`Bl` is the blow-up, `E` is the exceptional divisor, and `j: E → Bl` is the inclusion.
The structure maps `Bl → X` and `E → Z` are obtained via `structure_map`.

## Constructors

```@docs
blowup(i::AbstractVarietyMap; symbol::String = "e")
```

```@docs
blowup_points(X::AbstractVariety, n::Int; symbol::String = "e")
```

## Worked examples

### Blow-up of the Veronese surface

Blow up $\mathbb P^5$ along the Veronese surface (the image of $\mathbb P^2$
under the degree-2 map). The Chow ring generators, Betti numbers, tangent
bundle, and intersection numbers can be read off:

```jldoctest blowup_veronese
julia> P2 = abstract_projective_space(2);

julia> P5 = abstract_projective_space(5);

julia> i = map(P2, P5, [2*P2.O1]);

julia> Bl, E, j = blowup(i);

julia> betti_numbers(Bl)
6-element Vector{Int64}:
 1
 2
 3
 3
 2
 1

julia> euler_number(Bl)
12

julia> e, H = gens(chow_ring(Bl)) # exceptional divisor class and pullback of hyperplane class
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 e
 h

julia> integral(top_chern_class(tangent_bundle(Bl)))
12

```

We can verify that the proper transform of a quadric through the Veronese
has degree 1 (it is a linear space):

```jldoctest blowup_veronese
julia> e, H = gens(chow_ring(Bl));

julia> quad = pullback(structure_map(Bl), 2*P5.O1) - e;  # proper transform of a quadric

julia> integral(quad^5)
1

```

### Blow-up of a twisted cubic

Blow up $\mathbb P^3$ along a twisted cubic (the image of $\mathbb P^1$ via
$\mathcal{O}(3)$). The proper transform of a quadric containing the cubic
has self-intersection 0, while a cubic surface meets it in a single residual point:

```jldoctest
julia> P1 = abstract_projective_space(1);

julia> P3 = abstract_projective_space(3);

julia> i = map(P1, P3, [3*P1.O1]);

julia> Bl, E, j = blowup(i);

julia> e = pushforward(j, E(1));

julia> quad = pullback(structure_map(Bl), 2*P3.O1) - e;

julia> integral(quad^3)   # quadrics through the cubic form a net
0

julia> cubic = pullback(structure_map(Bl), 3*P3.O1) - e;

julia> integral(quad^2 * cubic)   # residual intersection
1

```

