@doc Markdown.doc"""
    underlying_scheme(X::ToricSpec)

For an affine toric scheme ``X``, this returns
the underlying scheme. In other words, by applying
this method, you obtain a scheme that has forgotten
its toric origin.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{QQFieldElem}[[1, 0], [0, 1]]

julia> underlying_scheme(affine_toric_scheme)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2 over Rational Field by ideal(0)
```
"""
function underlying_scheme(X::ToricSpec)
  if !isdefined(X, :X)
    if length(X.var_names) == 0
      I = toric_ideal(affine_normal_toric_variety(X))
      X.X = Spec(base_ring(I), I)
    else
      R, _ = polynomial_ring(QQ, X.var_names)
      I = toric_ideal(R, hilbert_basis(X))
      X.X = Spec(R, I)
    end
  end
  return X.X
end


@doc Markdown.doc"""
    affine_normal_toric_variety(X::ToricSpec)

For an affine toric scheme ``X``, this returns
the underlying affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{QQFieldElem}[[1, 0], [0, 1]]

julia> affine_normal_toric_variety(affine_toric_scheme)
Normal, affine toric variety
```
"""
affine_normal_toric_variety(X::ToricSpec) = X.antv


@doc Markdown.doc"""
    cone(X::ToricSpec)

For an affine toric scheme ``X``, this returns
the cone of the underlying affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{QQFieldElem}[[1, 0], [0, 1]]

julia> cone(affine_toric_scheme)
Polyhedral cone in ambient dimension 2
```
"""
cone(X::ToricSpec) = cone(affine_normal_toric_variety(X))


@doc Markdown.doc"""
    dual_cone(X::ToricSpec)

For an affine toric scheme ``X``, this returns the dual
of the cone of the underlying affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([1 0; 0 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{QQFieldElem}[[1, 0], [0, 1]]

julia> dual_cone(affine_toric_scheme)
Polyhedral cone in ambient dimension 2

julia> polarize(cone(affine_toric_scheme)) == dual_cone(affine_toric_scheme)
true
```
"""
function dual_cone(X::ToricSpec) 
  if !isdefined(X, :dual_cone)
    X.dual_cone = polarize(cone(X))
  end
  return X.dual_cone
end


@doc Markdown.doc"""
    hilbert_basis(X::ToricSpec)

For an affine toric scheme ``X``, this returns the Hilbert
basis of the cone dual to the cone of the underlying
affine toric variety.

# Examples
```jldoctest
julia> C = positive_hull([-1 1; 1 1])
Polyhedral cone in ambient dimension 2

julia> antv = affine_normal_toric_variety(C)
Normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{QQFieldElem}[[-1, 1], [1, 1]]

julia> dc = dual_cone(affine_toric_scheme)
Polyhedral cone in ambient dimension 2

julia> rays(dc)
2-element SubObjectIterator{RayVector{QQFieldElem}}:
 [1, 1]
 [-1, 1]

julia> hilbert_basis(affine_toric_scheme)
[-1   1]
[ 1   1]
[ 0   1]
```
"""
function hilbert_basis(X::ToricSpec)
  if !isdefined(X, :hb)
    X.hb = matrix(ZZ, hilbert_basis(dual_cone(X)))
  end
  return X.hb
end

