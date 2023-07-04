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
function intersect(C::Cone...)
    T, f = _promote_scalar_field((get_parent_field(c) for c in C)...)
    pmo = [pm_object(c) for c in C]
    return Cone{T}(Polymake.polytope.intersection(pmo...), f)
end
intersect(C::AbstractVector{<:Cone}) = intersect(C...)


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


@doc raw"""
    transform(C::Cone{T}, A::AbstractMatrix) where T<:scalar_types

Return the cone $A\cdot C$. The matrix $A$ must be square, full-rank, and
unimodular.

# Examples
Inversion at the origin:
```jldoctest
julia> c = positive_hull([1 0 0; 1 1 0; 1 1 1; 1 0 1])
Polyhedral cone in ambient dimension 3

julia> A = [-1 0 0; 0 -1 0; 0 0 -1]
3×3 Matrix{Int64}:
 -1   0   0
  0  -1   0
  0   0  -1

julia> ct = transform(c, A)
Polyhedral cone in ambient dimension 3

julia> rays(ct)
4-element SubObjectIterator{RayVector{QQFieldElem}}:
 [-1, 0, 0]
 [-1, -1, 0]
 [-1, -1, -1]
 [-1, 0, -1]

julia> ctt = transform(ct, A)
Polyhedral cone in ambient dimension 3

julia> c == ctt
true
```
"""
function transform(C::Cone{T}, A::Union{AbstractMatrix, MatElem{U}}) where {T<:scalar_types, U<:FieldElem}
  @assert ambient_dim(C) == nrows(A) "Incompatible dimension of cone and transformation matrix"
  @assert nrows(A) == ncols(A) "Transformation matrix must be square"
  @assert Polymake.common.rank(A) == nrows(A) "Transformation matrix must have full rank."
  return _transform(C, A)
end
function _transform(C::Cone{T}, A::AbstractMatrix) where T<:scalar_types
  OT = _scalar_type_to_polymake(T)
  raymod = Polymake.Matrix{OT}(permutedims(A))
  facetmod = Polymake.Matrix{OT}(Polymake.common.inv(permutedims(raymod)))
  return _transform(C, raymod, facetmod)
end
function _transform(C::Cone{T}, A::MatElem{U}) where {T<:scalar_types, U<:scalar_types}
  OT = _scalar_type_to_polymake(T)
  raymod = Polymake.Matrix{OT}(transpose(A))
  facetmod = Polymake.Matrix{OT}(inv(A))
  return _transform(C, raymod, facetmod)
end
function _transform(C::Cone{T}, raymod, facetmod) where T<:scalar_types
  OT = _scalar_type_to_polymake(T)
  result = Polymake.polytope.Cone{OT}()
  for prop in ("RAYS", "INPUT_RAYS", "LINEALITY_SPACE", "INPUT_LINEALITY")
    if Polymake.exists(pm_object(C), prop)
      resultprop = Polymake.Matrix{OT}(Polymake.give(pm_object(C), prop) * raymod)
      Polymake.take(result, prop, resultprop)
    end
  end
  for prop in ("INEQUALITIES", "EQUATIONS", "LINEAR_SPAN", "FACETS")
    if Polymake.exists(pm_object(C), prop)
      resultprop = Polymake.Matrix{OT}(Polymake.give(pm_object(C), prop) * facetmod)
      Polymake.take(result, prop, resultprop)
    end
  end
  return Cone{T}(result, get_parent_field(C))
end
