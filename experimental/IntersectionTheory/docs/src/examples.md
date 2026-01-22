```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Illustrating Examples From Enumerative Geometry

#### How Many Lines in $\mathbb P^3$ Meet Four General Lines in $\mathbb P^3$?

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

#### How many lines in $\mathbb P^3$ are secant to two General Twisted Cubic Curves in $\mathbb P^3$?

```jldoctest
julia> G = abstract_grassmannian(2, 4)
AbstractVariety of dim 4

julia> s2 = schubert_class(G, 2)
c[1]^2 - c[2]

julia> s11 = schubert_class(G, [1, 1])
c[2]

julia> integral((s2+3*s11)^2)
10

```

#### How Many Lines in $\mathbb P^4$ Meet six General Planes in $\mathbb P^4$?

```jldoctest
julia> G = abstract_grassmannian(2, 5)
AbstractVariety of dim 6

julia> s1 = schubert_class(G, 1)
-c[1]

julia> integral(s1^6)
5

```

#### How Many Planes in $\mathbb P^4$ Meet six General Lines in $\mathbb P^4$?

```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> s1 = schubert_class(G, 1)
-c[1]

julia> integral(s1^6)
5

```

#### Steiner's Problem: How Many Conics are Tangent to 5 General Conics in $\mathbb P^2$?

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

julia> Bl, E, j = blow_up(i)
(AbstractVariety of dim 5, AbstractVariety of dim 4, AbstractVarietyMap from AbstractVariety of dim 4 to AbstractVariety of dim 5)

julia> e, HBl = gens(chow_ring(Bl))
2-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 e
 H

julia> integral((6*HBl-2*e)^5)
3264

```

#### How Many Conics Meet Eight General Lines in $\mathbb P^3$?

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

#### How Many Conics lie on the General Cubic Surface and meet a General Line in $\mathbb P^3$?

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

julia> E = symmetric_power(USBd , 3)-USBd*OO(PF, -1)
AbstractBundle of rank 7 on AbstractVariety of dim 8

julia> C = top_chern_class(E)
27*z^5*c[2] - 135*z^4*c[3]

julia> integral(C*(2*p + z))
81

```

#### How Many Conics lie on the General Quintic Hypersurface in $\mathbb P^4$?

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

#### How Many Elliptic Cubics lie on the General Sextic Hypersurface in $\mathbb P^5$?

See [KP07](@cite).

```jldoctest
julia> G = abstract_grassmannian(3, 6);

julia> USBd = dual(tautological_bundles(G)[1])
AbstractBundle of rank 3 on AbstractVariety of dim 9

julia> F = symmetric_power(USBd, 3)
AbstractBundle of rank 10 on AbstractVariety of dim 9

julia> PF = projective_bundle(F)
AbstractVariety of dim 18

julia> E = symmetric_power(USBd, 6)-F*OO(PF, -1)
AbstractBundle of rank 18 on AbstractVariety of dim 18

julia> integral(top_chern_class(E))
2734099200

```

