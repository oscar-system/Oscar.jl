@doc raw"""
      generalized_cyclic_torus_cover(n::Int64, d::Int64, vslits::Vector, hslits::Vector)

Given integers ``d, n \ge 1``, this function returns an origami of degree ``n^{2d}`` in the following way: First,
for each ``i = 1, \dots, d`` we take ``n^2`` cornerless squares (``s_{i,1}, \dots, s_{i,n^2}``) and arrange them in
an ``n \times n``-grid indexed from left to right and then from the bottom to the top. The gluing is then described
by the two tuples `vslits` and `hslits` in ``(\mathbb{Z} / d \mathbb{Z})^{n^2}``: The right edge of a square
``s_{i,j}`` will be glued to the left edge of the right neighbor of ``s_{i+ \texttt{vslits}[j],j}``. Similarly,
`hslits` sets the gluing of the upper edges.

# Examples
```jldoctest
julia> generalized_cyclic_torus_cover(2, 2, [1,0,0,0], [0,0,0,0])
GAP: Origami((1,6,5,2)(3,4)(7,8), (1,3)(2,4)(5,7)(6,8), 8)
```
"""
function generalized_cyclic_torus_cover(n::Int64, d::Int64, vslits::Vector, hslits::Vector)
  return GAP.Globals.GeneralizedCyclicTorusCover(n, d, GapObj(vslits), GapObj(hslits))
end

@doc raw"""
      comb_origami(n::Int64, x::Int64, y::Int64)

Returns a cyclic torus cover of degree ``d = 2`` whose critical points are exactly the point ``P = (x,y)`` and the
three copies of ``P`` that are obtained by rotating ``P`` around the center of the ``n``-torus. The coordinates are
given in the range ``0, \dots, (n - 1)^2``, where the point ``(0,0)`` is located in the lower left corner.
``P = (x,y)`` must not be a 2-torsion point, that is, it must not be ``(0,0)``, ``(n/2, n/2)``, ``(n/2,0)`` or
``(0,n/2)``. The coordinates are considered modulo ``n``. See [HS07](@cite) for more details.

# Examples
```jldoctest
julia> comb_origami(3,0,1)
GAP: Origami((1,2,3)(4,5,6,13,14,15)(7,8,9)(10,11,12)(16,17,18), (1,4,7)(2,5,8,11,14,17)(3,6,9)(10,13,16)(12,15,18), 18)
```
"""
function comb_origami(n::Int64, x::Int64, y::Int64)
  return GAP.Globals.CombOrigami(n, x, y)
end

@doc raw"""
      cyclic_torus_cover_origamiS(n::Int64, d::Int64, v::Vector)

Returns: a cyclic torus cover origami whose monodromy vector with respect to the basis ``S`` is `v`.

# Examples
```jldoctest
julia> cyclic_torus_cover_origamiS(2,2,[1,0,1,0,0])
GAP: Origami((1,2,5,6)(3,4)(7,8), (1,3,5,7)(2,4)(6,8), 8)
```
"""
function cyclic_torus_cover_origamiS(n::Int64, d::Int64, v::Vector)
  return GAP.Globals.CyclicTorusCoverOrigamiS(n, d, GapObj(v))
end

@doc raw"""
      cyclic_torus_cover_origamiL(n::Int64, d::Int64, v::Vector)

Returns: a cyclic torus cover origami whose monodromy vector with respect to the basis ``L`` is `v`.

# Examples
```jldoctest
julia> cyclic_torus_cover_origamiL(2,2,[1,0,1,0,0])
GAP: Origami((1,2,5,6)(3,4)(7,8), (1,7)(2,4,6,8)(3,5), 8)
```
"""
function cyclic_torus_cover_origamiL(n::Int64, d::Int64, v::Vector)
  return GAP.Globals.CyclicTorusCoverOrigamiL(n, d, GapObj(v))
end

@doc raw"""
      base_change_l_to_s(n::Int64)

Returns a matrix corresponding to the change of basis from ``L`` to ``S`` on the homology of ``T_n``. The matrix has
the following property: given any cyclic torus cover origami as a monodromy vector ``v`` with respect to the basis
``S``, you may obtain the corresponding monodromy vector with respect to basis ``L`` using ``v \cdot D_{SL}``.

# Examples
```jldoctest
julia> base_change_l_to_s(2)
GAP: [ [ 0, 1, 1, -1, 0 ], [ 0, 0, -1, 1, 0 ], [ 1, 0, -1, 0, 1 ], [ 0, 0, 1, 0, -1 ], [ 0, 1, 0, 0, 1 ] ]
```
"""
function base_change_l_to_s(n::Int64)
  return GAP.Globals.BaseChangeLToS(n)
end

function translation_group_on_homology_of_tn(n::Int64)
  return GAP.Globals.TranslationGroupOnHomologyOfTn(n)
end

function action_of_t_on_homology_of_tn(n::Int64)
  return GAP.Globals.ActionOfTOnHomologyOfTn(n)
end

function action_of_s_on_homology_of_tn(n::Int64)
  return GAP.Globals.ActionOfSOnHomologyOfTn(n)
end

function action_of_matrix_on_homology_of_tn(n::Int64, A::Matrix)
  return GAP.Globals.ActionOfMatrixOnHomologyOfTn(n, GapObj(A))
end
