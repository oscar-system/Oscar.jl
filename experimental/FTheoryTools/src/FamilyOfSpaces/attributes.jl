#####################################################
# 1 Basic attributes
#####################################################

@doc raw"""
    coordinate_ring(f::FamilyOfSpaces)

Return the coordinate ring of a generic member of the family of spaces.

```jldoctest
julia> ring, (f, g, Kbar, u) = QQ["f", "g", "Kbar", "u"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[f, g, Kbar, u])

julia> grading = [4 6 1 0]
1×4 Matrix{Int64}:
 4  6  1  0

julia> d = 3
3

julia> f = family_of_spaces(ring, grading, d)
A family of spaces of dimension d = 3

julia> coordinate_ring(f)
Multivariate polynomial ring in 4 variables f, g, Kbar, u
  over rational field
```
"""
coordinate_ring(f::FamilyOfSpaces) = f.coordinate_ring


@doc raw"""
    weights(f::FamilyOfSpaces)

Return the grading of the coordinate ring of a generic member of the family of spaces.

```jldoctest
julia> ring, (f, g, Kbar, u) = QQ["f", "g", "Kbar", "u"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[f, g, Kbar, u])

julia> grading = [4 6 1 0]
1×4 Matrix{Int64}:
 4  6  1  0

julia> d = 3
3

julia> f = family_of_spaces(ring, grading, d)
A family of spaces of dimension d = 3

julia> weights(f)
1×4 Matrix{Int64}:
 4  6  1  0
```
"""
weights(f::FamilyOfSpaces) = f.grading


@doc raw"""
    dim(f::FamilyOfSpaces)

Return the dimension of the generic member of the family of spaces.

```jldoctest
julia> ring, (f, g, Kbar, u) = QQ["f", "g", "Kbar", "u"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[f, g, Kbar, u])

julia> grading = [4 6 1 0]
1×4 Matrix{Int64}:
 4  6  1  0

julia> d = 3
3

julia> f = family_of_spaces(ring, grading, d)
A family of spaces of dimension d = 3

julia> dim(f)
3
```
"""
dim(f::FamilyOfSpaces) = f.dim


#####################################################
# 2 Advanced attributes
#####################################################

@doc raw"""
    stanley_reisner_ideal(f::FamilyOfSpaces)

Return the equivalent of the Stanley-Reisner ideal
for the generic member of the family of spaces.

```jldoctest
julia> coordinate_ring, (f, g, Kbar, u) = QQ["f", "g", "Kbar", "u"]
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[f, g, Kbar, u])

julia> grading = [4 6 1 0]
1×4 Matrix{Int64}:
 4  6  1  0

julia> d = 3
3

julia> f = family_of_spaces(coordinate_ring, grading, d)
A family of spaces of dimension d = 3

julia> stanley_reisner_ideal(f)
ideal(f*g*Kbar, f*g*u, f*Kbar*u, g*Kbar*u)
```
"""
@attr MPolyIdeal{QQMPolyRingElem} function stanley_reisner_ideal(f::FamilyOfSpaces)
  ring = coordinate_ring(f)
  variables = gens(ring)
  combis = Oscar.combinations(length(variables), dim(f))
  ideal_generators = [prod([variables[i] for i in c]) for c in combis]
  return ideal(ideal_generators)
end
