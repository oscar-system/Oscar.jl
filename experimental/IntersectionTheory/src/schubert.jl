@doc raw"""
    schubert_class(G::AbstractVariety, λ::Int...)
    schubert_class(G::AbstractVariety, λ::Vector{Int})
    schubert_class(G::AbstractVariety, λ::Partition)

Return the Schubert class $\sigma_\lambda$ on a (relative) Grassmannian `G`.

# Examples

```jldoctest
julia> G = abstract_grassmannian(2,4)
AbstractVariety of dim 4

julia> s0 = schubert_class(G, 0)
1

julia> s1 = schubert_class(G, 1)
-c[1]

julia> s2 = schubert_class(G, 2)
c[1]^2 - c[2]

julia> s11 = schubert_class(G, [1, 1])
c[2]

julia> s21 = schubert_class(G, [2, 1])
-c[1]*c[2]

julia> s22 = schubert_class(G, [2, 2])
c[2]^2

julia> s1*s1 == s2+s11
true

julia> s1*s2 == s1*s11 == s21
true

julia> s1*s21 == s2*s2 == s22
true

```

```jldoctest
julia> G = abstract_grassmannian(2,5)
AbstractVariety of dim 6

julia> s3 = schubert_class(G, 5-2)
-c[1]^3 + 2*c[1]*c[2]

julia> s3^2 == point_class(G)
true

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 3 on AbstractVariety of dim 6

julia> chern_class(Q, 3)
-c[1]^3 + 2*c[1]*c[2]

```

```jldoctest
julia> G = abstract_grassmannian(2,4)
AbstractVariety of dim 4

julia> s1 = schubert_class(G, 1)
-c[1]

julia> s1 == schubert_class(G, [1, 0])
true

julia> integral(s1^4)
2

```
"""
function schubert_class(G::AbstractVariety, λ::Int...) schubert_class(G, collect(λ)) end
function schubert_class(G::AbstractVariety, λ::Partition) schubert_class(G, Vector(λ)) end
function schubert_class(G::AbstractVariety, λ::Vector{Int})
  get_attribute(G, :grassmannian) === nothing && error("the abstract_variety is not a Grassmannian")
  (length(λ) > rank(G.bundles[1]) || sort(λ, rev=true) != λ) && error("the Schubert input is not well-formed")
  giambelli(G.bundles[2], λ)
end

@doc raw"""
    schubert_classes(G::AbstractVariety, m::Int)

Return all Schubert classes in codimension `m` on a (relative) Grassmannian `G`.
"""
function schubert_classes(G::AbstractVariety, m::Int)
  get_attribute(G, :grassmannian) === nothing && error("the abstract_variety is not a Grassmannian")
  S, Q = G.bundles
  res = elem_type(G.ring)[]
  for i in 0:rank(S)
    append!(res, [schubert_class(G, l) for l in partitions(m, i, 1, rank(Q))])
  end
  return res
end

@doc raw"""
    schubert_classes(G::AbstractVariety)

Return all Schubert classes on a (relative) Grassmannian `G`.

# Examples

```jldoctest
julia> G = abstract_grassmannian(2,4)
AbstractVariety of dim 4

julia> schubert_classes(G)
5-element Vector{Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}:
 [1]
 [-c[1]]
 [c[1]^2 - c[2], c[2]]
 [-c[1]*c[2]]
 [c[2]^2]

julia> basis(G)
5-element Vector{Vector{MPolyQuoRingElem}}:
 [1]
 [c[1]]
 [c[2], c[1]^2]
 [c[1]*c[2]]
 [c[2]^2]
 
```
"""
function schubert_classes(G::AbstractVariety)
   get_attribute(G, :grassmannian) === nothing && error("the abstract_variety is not a Grassmannian")
   S, Q = G.bundles
   return [schubert_classes(G, i) for i = 0:rank(S)*rank(Q)]
end
