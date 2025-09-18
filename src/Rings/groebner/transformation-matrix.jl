@doc raw"""
    _compute_standard_basis_with_transform(B::IdealGens, ordering::MonomialOrdering, complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `IdealGens` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a standard basis (if `ordering` is a global monomial ordering and `complete_reduction = true`
the reduced Gröbner basis) of the ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering`
and the transformation matrix from the ideal to the standard basis. Return value is a IdealGens together with a map.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> A = Oscar.IdealGens([x*y-3*x,y^3-2*x^2*y])
Ideal generating system with elements
  1: x*y - 3*x
  2: -2*x^2*y + y^3

julia> B,m = Oscar._compute_standard_basis_with_transform(A, degrevlex(R))
(Ideal generating system with 3 elements with associated ordering degrevlex([x, y]), [1 2*x -2*x^2+y^2+3*y+9; 0 1 -x])
```
"""
function _compute_standard_basis_with_transform(B::IdealGens, ordering::MonomialOrdering, complete_reduction::Bool = false)
  istd, m = Singular.lift_std(singular_generators(B, ordering), complete_reduction = complete_reduction)
  R = base_ring(B)
  return IdealGens(R, istd), map_entries(R, m)
end

@doc raw"""
    standard_basis_with_transformation_matrix(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool=false)

Return a pair `G`, `T`, say, where `G` is a standard basis of `I` with respect to `ordering`, and `T` 
is a transformation matrix from `gens(I)` to `G`. That is, `gens(I)*T == G`.

!!! note
    The returned Gröbner basis is reduced if `ordering` is a global monomial ordering and `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ,[:x,:y]);

julia> I = ideal([x*y^2-1,x^3+y^2+x*y]);

julia> G, T = standard_basis_with_transformation_matrix(I, ordering=neglex(R))
(Standard basis with 1 element w.r.t. neglex([x, y]), [-1; 0])

julia> gens(I)*T == gens(G)
true
```
"""
function standard_basis_with_transformation_matrix(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool = false)
  complete_reduction && @assert is_global(ordering)
  G, m = _compute_standard_basis_with_transform(I.gens, ordering, complete_reduction)
  G.isGB = true
  I.gb[ordering]  = G
  return G, m
end

@doc raw"""
    groebner_basis_with_transformation_matrix(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool=false)

Return a pair `G`, `T`, say, where `G` is a Gröbner basis of `I` with respect to `ordering`, and `T` 
is a transformation matrix from `gens(I)` to `G`. That is, `gens(I)*T == G`.

!!! note
    The returned Gröbner basis is reduced if `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ,[:x,:y]);

julia> I = ideal([x*y^2-1,x^3+y^2+x*y]);

julia> G, T = groebner_basis_with_transformation_matrix(I)
(Gröbner basis with 3 elements w.r.t. degrevlex([x, y]), [1 0 -x^2-y; 0 1 y^2])

julia> gens(I)*T == gens(G)
true
```
"""
function groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool = false)
    is_global(ordering) || error("Ordering must be global")
    return standard_basis_with_transformation_matrix(I, ordering=ordering, complete_reduction=complete_reduction)
end
