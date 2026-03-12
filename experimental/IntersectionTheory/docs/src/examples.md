```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Illustrating examples from enumerative geometry

## How many lines in $\mathbb P^3$ meet four general lines in $\mathbb P^3$?

This is Example 14.7.2 in [Ful98](@cite).

```jldoctest
julia> G = abstract_grassmannian(2, 4)
AbstractVariety of dim 4

julia> schubert_classes(G)
5-element Vector{Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 [1]
 [-c[1]]
 [c[1]^2 - c[2], c[2]]
 [-c[1]*c[2]]
 [c[2]^2]

julia> s1 = schubert_class(G, 1)
-c[1]

julia> integral(s1^4)
2

```

## How many lines in $\mathbb P^3$ are secant to two general twisted cubic curves in $\mathbb P^3$?

This is discussed in Section 3.4 of [EH16](@cite).

```jldoctest
julia> G = abstract_grassmannian(2, 4)
AbstractVariety of dim 4

julia> s2 = schubert_class(G, 2)
c[1]^2 - c[2]

julia> s11 = schubert_class(G, [1, 1])
c[2]

julia> integral((s2 + 3*s11)^2)
10

```

## How many lines in $\mathbb P^4$ meet six general planes in $\mathbb P^4$?

```jldoctest
julia> G = abstract_grassmannian(2, 5)
AbstractVariety of dim 6

julia> s1 = schubert_class(G, 1)
-c[1]

julia> integral(s1^6)
5

```

## How many planes in $\mathbb P^4$ meet six general lines in $\mathbb P^4$?

```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> s1 = schubert_class(G, 1)
-c[1]

julia> integral(s1^6)
5

```

## Steiner's problem: How many conics are tangent to 5 general conics in $\mathbb P^2$?

This is the question that is the explanation for the name of [EH16](@cite),
and corresponds to Theorem 8.9.

```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P5 = abstract_projective_space(5, symbol = "H")
AbstractVariety of dim 5

julia> h = gens(P2)[1]
h

julia> H = gens(P5)[1]
H

julia> i = map(P2, P5, [2*h])
AbstractVarietyMap from AbstractVariety of dim 2 to AbstractVariety of dim 5

julia> Bl, E, j = blowup(i)
(AbstractVariety of dim 5, AbstractVariety of dim 4, AbstractVarietyMap from AbstractVariety of dim 4 to AbstractVariety of dim 5)

julia> e, HBl = gens(chow_ring(Bl))
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 e
 H

julia> integral((6*HBl - 2*e)^5)
3264

```

## How many conics meet eight general lines in $\mathbb P^3$?

This is Example 14.7.12 in [Ful98](@cite).

```jldoctest
julia> G = abstract_grassmannian(3, 4)
AbstractVariety of dim 3

julia> USBd = dual(tautological_bundles(G)[1])
AbstractBundle of rank 3 on AbstractVariety of dim 3

julia> F = symmetric_power(USBd, 2)
AbstractBundle of rank 6 on AbstractVariety of dim 3

julia> PF = projective_bundle(F) # the parameter space of conics in P3
AbstractVariety of dim 8

julia> UQB = tautological_bundles(G)[2]
AbstractBundle of rank 1 on AbstractVariety of dim 3

julia> p = pullback(structure_map(PF), chern_class(UQB, 1))
-c[1]

julia> Z = dual(tautological_bundles(PF)[1])
AbstractBundle of rank 1 on AbstractVariety of dim 8

julia> z = chern_class(Z, 1)
z

julia> integral((2*p + z)^8)
92

```

## How many conics lie on the general cubic surface and meet a general line in $\mathbb P^3$?

```jldoctest
julia> G = abstract_grassmannian(3, 4)
AbstractVariety of dim 3

julia> USBd = dual(tautological_bundles(G)[1])
AbstractBundle of rank 3 on AbstractVariety of dim 3

julia> F = symmetric_power(USBd, 2)
AbstractBundle of rank 6 on AbstractVariety of dim 3

julia> PF = projective_bundle(F)
AbstractVariety of dim 8

julia> UQB = tautological_bundles(G)[2]
AbstractBundle of rank 1 on AbstractVariety of dim 3

julia> p = pullback(structure_map(PF), chern_class(UQB, 1))
-c[1]

julia> Z = dual(tautological_bundles(PF)[1])
AbstractBundle of rank 1 on AbstractVariety of dim 8

julia> z = chern_class(Z, 1)
z

julia> E = symmetric_power(USBd, 3)-USBd*OO(PF, -1)
AbstractBundle of rank 7 on AbstractVariety of dim 8

julia> C = top_chern_class(E)
27*z^5*c[2] - 135*z^4*c[3]

julia> integral(C*(2*p + z))
81

```

## How many conics lie on the general quintic hypersurface in $\mathbb P^4$?

This is a famous question, relevant for the development of mirror symmetry.
See Theorem 3.1 in [K86](@cite).

```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> USBd = dual(tautological_bundles(G)[1])
AbstractBundle of rank 3 on AbstractVariety of dim 6

julia> F = symmetric_power(USBd, 2)
AbstractBundle of rank 6 on AbstractVariety of dim 6

julia> PF = projective_bundle(F) # the parameter space of conics in P4
AbstractVariety of dim 11

julia> A = symmetric_power(USBd, 5) - symmetric_power(USBd, 3)*OO(PF, -1)
AbstractBundle of rank 11 on AbstractVariety of dim 11

julia> integral(top_chern_class(A))
609250

```

## How many elliptic cubics lie on the general sextic hypersurface in $\mathbb P^5$?

See Table 2 in [KP07](@cite).

```jldoctest
julia> G = abstract_grassmannian(3, 6);

julia> USBd = dual(tautological_bundles(G)[1])
AbstractBundle of rank 3 on AbstractVariety of dim 9

julia> F = symmetric_power(USBd, 3)
AbstractBundle of rank 10 on AbstractVariety of dim 9

julia> PF = projective_bundle(F)
AbstractVariety of dim 18

julia> E = symmetric_power(USBd, 6) - F*OO(PF, -1)
AbstractBundle of rank 18 on AbstractVariety of dim 18

julia> integral(top_chern_class(E))
2734099200

```

## Euler characteristics of $\mathcal{O}_{\mathbb P^3}(n)$

The Euler characteristic $\chi(\mathcal{O}_{\mathbb P^n}(d)) = \binom{n+d}{n}$
is computed by Hirzebruch–Riemann–Roch. Here we tabulate it for $\mathbb P^3$
over a range of twists, mirroring the Sage `Proj(3).o(n).euler_characteristic()` example:

```jldoctest
julia> P3 = abstract_projective_space(3);

julia> [euler_characteristic(OO(P3, n)) for n in 0:6]
7-element Vector{QQFieldElem}:
 1
 4
 10
 20
 35
 56
 84

```

## Grassmannian duality: $\mathrm{Gr}(2,5) \cong \mathrm{Gr}(3,5)$

The Grassmannians $\mathrm{Gr}(k,n)$ and $\mathrm{Gr}(n-k,n)$ are isomorphic.
In particular they have the same Betti numbers, Euler number, and basis
structure. The universal sub-bundle of one corresponds to the dual of the
universal quotient bundle of the other:

```jldoctest
julia> G1 = abstract_grassmannian(2, 5);

julia> G2 = abstract_grassmannian(3, 5);

julia> betti_numbers(G1) == betti_numbers(G2)
true

julia> euler_number(G1)
10

julia> euler_number(G2)
10

julia> integral(top_chern_class(symmetric_power(dual(tautological_bundles(G1)[1]), 3)))
27

julia> integral(top_chern_class(symmetric_power(dual(tautological_bundles(G2)[1]), 3)))
27

```

## Blow-up of the Veronese surface: invariants

Blow up $\mathbb P^5$ along the Veronese image of $\mathbb P^2$ (degree-2 embedding).
We compute the Chow ring generators, Betti numbers, tangent bundle Chern class,
and several intersection numbers matching the Sage `blowup` examples:

```jldoctest
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

julia> integral(top_chern_class(tangent_bundle(Bl)))  # Euler number via top Chern class
12

julia> e = pushforward(j, E(1)); quad = pullback(structure_map(Bl), 2*P5.O1) - e;

julia> integral(quad^5)  # proper transform of a quadric through the Veronese
1

julia> sext = pullback(structure_map(Bl), 6*P5.O1) - 2*e;

julia> integral(sext^5)  # Steiner's 3264
3264

```

## Blow-up of a twisted cubic: quadric and cubic divisors

Blow up $\mathbb P^3$ along a twisted cubic. The proper transform of a quadric
containing the cubic has zero self-intersection, and a cubic surface meets
the proper transform of a pair of quadrics in exactly one residual point:

```jldoctest
julia> P1 = abstract_projective_space(1);

julia> P3 = abstract_projective_space(3);

julia> i = map(P1, P3, [3*P1.O1]);

julia> Bl, E, j = blowup(i);

julia> e = pushforward(j, E(1));

julia> quad = pullback(structure_map(Bl), 2*P3.O1) - e;

julia> integral(quad^3)
0

julia> cubic = pullback(structure_map(Bl), 3*P3.O1) - e;

julia> integral(quad^2 * cubic)
1

```

With parameters $r,s,t$ for the degrees of three divisors, the triple
intersection number is a polynomial in $r,s,t$:

```jldoctest
julia> T, (r, s, t) = polynomial_ring(QQ, [:r, :s, :t]);

julia> F = fraction_field(T);

julia> (r, s, t) = gens(F);

julia> P1 = abstract_projective_space(1; base=F);

julia> P3 = abstract_projective_space(3; base=F);

julia> i = map(P1, P3, [3*P1.O1]);

julia> Bl, E, j = blowup(i);

julia> e = pushforward(j, E(1));

julia> rH, sH, tH = [pullback(structure_map(Bl), x * P3.O1) - e for x in [r, s, t]];

julia> integral(rH * sH * tH)
r*s*t - 3*r - 3*s - 3*t + 10

```

## Segre blow-up: $\mathbb P^2 \times \mathbb P^2 \hookrightarrow \mathbb P^8$

Blow up $\mathbb P^8$ along the Segre image of $\mathbb P^2 \times \mathbb P^2$.
We verify the Betti numbers, Euler number, and the generators of the Chow ring:

```jldoctest
julia> P2xP2 = abstract_projective_space(2, symbol = "k") * abstract_projective_space(2, symbol = "l");

julia> P8 = abstract_projective_space(8);

julia> k, l = gens(P2xP2);

julia> Se = map(P2xP2, P8, [k + l]);

julia> Bl, E, j = blowup(Se);

julia> betti_numbers(Bl)
9-element Vector{Int64}:
 1
 2
 4
 7
 8
 7
 4
 2
 1

julia> euler_number(Bl)
36

julia> length(gens(chow_ring(Bl)))  # number of generators
5

```
