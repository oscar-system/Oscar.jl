# see https://arxiv.org/pdf/1809.10327.pdf Lemma 5.4

@doc raw"""
    action_s(o::Origami)

For a given origami ``O`` this methods computes the origami ``S \cdot O``, where
``S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}`` is one of the two generators of ``{\rm SL}_2(\mathbb{Z})``. 

# Examples
```jldoctest
julia> o = origami(cperm([1,6,4,7,5,3], [2,8]), cperm([1,4,5,3,8,2,6]))
Origami ((1,6,4,7,5,3)(2,8),(1,4,5,3,8,2,6), 8)

julia> action_s(o)
Origami ((1,6,2,8,3,5,4),(1,6,4,7,5,3)(2,8), 8)
```
"""
function action_s(o::Origami)
  h = horizontal_perm(o)
  v = vertical_perm(o)
  return origami_disconnected(v^-1, h, degree(o))
end

@doc raw"""
    action_t(o::Origami)

For a given origami ``O`` this methods computes the origami ``T \cdot O``, where
``T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}`` is one of the two generators of ``{\rm SL}_2(\mathbb{Z})``. 

# Examples
```jldoctest
julia> o = origami(cperm([1,6,4,7,5,3], [2,8]), cperm([1,4,5,3,8,2,6]))
Origami ((1,6,4,7,5,3)(2,8),(1,4,5,3,8,2,6), 8)

julia> action_t(o)
Origami ((1,6,4,7,5,3)(2,8),(1,8,6,4)(5,7), 8)
```
"""
function action_t(o::Origami)
  h = horizontal_perm(o)
  return origami_disconnected(h, h^-1 * vertical_perm(o), degree(o))
end

@doc raw"""
    action_t_inv(o::Origami)

For a given origami ``O`` this methods computes the origami ``T^{-1} \cdot O``, where
``T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}`` is one of the two generators of ``{\rm SL}_2(\mathbb{Z})``.

# Examples
```jldoctest
julia> o = origami(cperm([1,6,4,7,5,3], [2,8]), cperm([1,4,5,3,8,2,6]))
Origami ((1,6,4,7,5,3)(2,8),(1,4,5,3,8,2,6), 8)

julia> action_t_inv(o)
Origami ((1,6,4,7,5,3)(2,8),(3,4,7)(5,8,6), 8)
```
"""
function action_t_inv(o::Origami)
  h = horizontal_perm(o)
  v = vertical_perm(o)
  return origami_disconnected(h, h * v, degree(o))
end

@doc raw"""
    action_s_inv(o::Origami)

For a given origami ``O`` this methods computes the origami ``S^{-1} \cdot O``, where
``S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}`` is one of the two generators of ``{\rm SL}_2(\mathbb{Z})``. 

# Examples
```jldoctest
julia> o = origami(cperm([1,6,4,7,5,3], [2,8]), cperm([1,4,5,3,8,2,6]))
Origami ((1,6,4,7,5,3)(2,8),(1,4,5,3,8,2,6), 8)

julia> action_s_inv(o)
Origami ((1,4,5,3,8,2,6),(1,3,5,7,4,6)(2,8), 8)
```
"""
function action_s_inv(o::Origami)
  h = horizontal_perm(o)
  v = vertical_perm(o)
  return origami_disconnected(v, h^-1, degree(o))
end

@doc raw"""
    action_sl2(A::ZZMatrix,o::Origami)

This function computes ``A \cdot O``, where ``A`` is a matrix in ``{\rm SL}_2(\mathbb{Z})`` and ``O`` is an origami.

# Examples
```jldoctest
julia> o = origami(cperm([1,6,4,7,5,3], [2,8]), cperm([1,4,5,3,8,2,6]))
Origami ((1,6,4,7,5,3)(2,8),(1,4,5,3,8,2,6), 8)

julia> action_sl2(ZZ[0 -1; 1 1], o)
Origami ((1,4,6,8)(5,7),(1,6,4,7,5,3)(2,8), 8)
```
"""
function action_sl2(A::ZZMatrix, o::Origami)
  @req is_one(det(A)) "matrix has to be an element of SL2(Z)"

  # TODO implement this in Oscar when implementing the modular subgroup pkg
  # because this is really slow and ugly
  gap_origami = GAP.Globals.ActionOfSL2(GapObj(A), GapObj(o))
  d = GAP.Globals.DegreeOrigami(gap_origami)::Int
  G = symmetric_group(d)
  h = PermGroupElem(G, GAP.Globals.HorizontalPerm(gap_origami))
  v = PermGroupElem(G, GAP.Globals.VerticalPerm(gap_origami))
  return origami_disconnected(h, v, degree(o))
end
