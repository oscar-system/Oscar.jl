@doc Markdown.doc"""
    underlying_scheme(X::ToricSpec)

For an affine toric scheme ``X``, this returns
the underlying scheme. In other words, by applying
this method, you obtain a scheme that has forgotten
its toric origin.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[1, 0], [0, 1]]

julia> underlying_scheme(affine_toric_scheme)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2 over Rational Field by ideal()
```
"""
underlying_scheme(X::ToricSpec) = X.X
export underlying_scheme


@doc Markdown.doc"""
    affine_normal_toric_variety(X::ToricSpec)

For an affine toric scheme ``X``, this returns
the underlying affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[1, 0], [0, 1]]

julia> affine_normal_toric_variety(affine_toric_scheme)
A normal, affine toric variety
```
"""
affine_normal_toric_variety(X::ToricSpec) = X.antv
export affine_normal_toric_variety


@doc Markdown.doc"""
    cone(X::ToricSpec)

For an affine toric scheme ``X``, this returns
the cone of the underlying affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[1, 0], [0, 1]]

julia> cone(affine_toric_scheme)
A polyhedral cone in ambient dimension 2
```
"""
cone(X::ToricSpec) = cone(affine_normal_toric_variety(X))
export cone


@doc Markdown.doc"""
    dual_cone(X::ToricSpec)

For an affine toric scheme ``X``, this returns the dual
of the cone of the underlying affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[1, 0], [0, 1]]

julia> dual_cone(affine_toric_scheme)
A polyhedral cone in ambient dimension 2

julia> polarize(cone(affine_toric_scheme)) == dual_cone(affine_toric_scheme)
true
```
"""
dual_cone(X::ToricSpec) = X.dual_cone
export dual_cone


@doc Markdown.doc"""
    hilbert_basis(X::ToricSpec)

For an affine toric scheme ``X``, this returns the Hilbert
basis of the cone dual to the cone of the underlying
affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([-1 1; 1 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[-1, 1], [1, 1]]

julia> dc = dual_cone(affine_toric_scheme)
A polyhedral cone in ambient dimension 2

julia> rays(dc)
2-element SubObjectIterator{RayVector{fmpq}}:
 [1, 1]
 [-1, 1]

julia> hilbert_basis(affine_toric_scheme)
[-1   1]
[ 1   1]
[ 0   1]
```
"""
hilbert_basis(X::ToricSpec) = X.hb
export hilbert_basis
