###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################
@doc raw"""
    intersect(C0::Cone{T}, C1::Cone{T}) where T<:scalar_types

Return the intersection $C0 \cap C1$ of `C0` and `C1`.

# Examples
```jldoctest
julia> C0 = positive_hull([1 0])
Polyhedral cone in ambient dimension 2

julia> C1 = positive_hull([0 1])
Polyhedral cone in ambient dimension 2

julia> C01 = intersect(C0, C1)
Polyhedral cone in ambient dimension 2

julia> rays(C01)
0-element SubObjectIterator{RayVector{QQFieldElem}}

julia> dim(C01)
0
```
"""
function intersect(C::Cone{T}...) where T<:scalar_types
    pmo = [pm_object(c) for c in C]
    return Cone{T}(Polymake.polytope.intersection(pmo...), get_parent_field(C[1]))
end
intersect(C::AbstractVector{Cone{T}}) where T<:scalar_types = intersect(C...)


@doc raw"""
    polarize(C::Cone)

Return the polar dual of `C`, the cone consisting of all those linear functions
that evaluate positively on all of `C`.

# Examples
```jldoctest
julia> C = positive_hull([1 0; -1 2])
Polyhedral cone in ambient dimension 2

julia> Cv = polarize(C)
Polyhedral cone in ambient dimension 2

julia> rays(Cv)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 1//2]
 [0, 1]
```
"""
function polarize(C::Cone{T}) where T<:scalar_types
    return Cone{T}(Polymake.polytope.polarize(pm_object(C)), get_parent_field(C))
end
