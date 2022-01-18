###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################
@doc Markdown.doc"""
    intersect(C0::Cone, C1::Cone)

Return the intersection $C0 \cap C1$ of `C0` and `C1`.

# Examples
```jldoctest
julia> C0 = positive_hull([1 0])
A polyhedral cone in ambient dimension 2

julia> C1 = positive_hull([0 1])
A polyhedral cone in ambient dimension 2

julia> C01 = intersect(C0, C1)
A polyhedral cone in ambient dimension 2

julia> rays(C01)
0-element SubObjectIterator{RayVector{Polymake.Rational}}

julia> dim(C01)
0
```
"""
function intersect(C0::Cone, C1::Cone)
   return Cone(Polymake.polytope.intersection(pm_object(C0), pm_object(C1)))
end


@doc Markdown.doc"""
    polarize(C::Cone)

Return the dual cone of `C` consisting of all those linear functions that
evaluate positively on all of `C`.

# Examples
```jldoctest
julia> C = positive_hull([1 0; -1 2])
A polyhedral cone in ambient dimension 2

julia> Cv = polarize(C)
A polyhedral cone in ambient dimension 2

julia> rays(Cv)
2-element SubObjectIterator{RayVector{Polymake.Rational}}:
 [1, 1/2]
 [0, 1]
```
"""
function polarize(C::Cone)
    return Cone(Polymake.polytope.polarize(pm_object(C)))
end
