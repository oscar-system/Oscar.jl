# see https://arxiv.org/pdf/1809.10327.pdf Lemma 5.4

@doc raw"""
    action_s(o::Origami)

For a given origami `o` this methods computes the origami ``S \cdot \texttt{o}``, where
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
  return origami(v^-1, h; check=false)
end

@doc raw"""
    action_t(o::Origami)

For a given origami `o` this methods computes the origami ``T \cdot \texttt{o}``, where
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
  v = vertical_perm(o)
  return origami(h, h^-1 * v; check=false)
end

@doc raw"""
    action_t_inv(o::Origami)

For a given origami `o` this methods computes the origami ``T^{-1} \cdot \texttt{o}``, where
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
  return origami(h, h * v; check=false)
end

@doc raw"""
    action_s_inv(o::Origami)

For a given origami `o` this methods computes the origami ``S^{-1} \cdot \texttt{o}``, where
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
  return origami(v, h^-1; check=false)
end

@doc raw"""
    action_sl2(A::ZZMatrix, o::Origami)

This function computes ``\texttt{A} \cdot \texttt{o}``, where `A` is a matrix in ``{\rm SL}_2(\mathbb{Z})`` and `o` is an origami.
The matrix is decomposed into a word in the generators
```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad
T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix},
```
as by [`s_t_decomposition`](@ref), and the letters of that word are applied from the right.

# Examples
```jldoctest
julia> o = origami(cperm([1,6,4,7,5,3], [2,8]), cperm([1,4,5,3,8,2,6]))
Origami ((1,6,4,7,5,3)(2,8),(1,4,5,3,8,2,6), 8)

julia> action_sl2(ZZ[0 -1; 1 1], o)
Origami ((1,4,6,8)(5,7),(1,6,4,7,5,3)(2,8), 8)
```
"""
function action_sl2(A::ZZMatrix, o::Origami)
  @req size(A) == (2, 2) && isone(det(A)) "matrix has to be an element of SL2(Z)"

  # Syllables of `S` (generator 1) and `T` (generator 2), read from the right,
  # so that a product acts on the left. This is GAP's `ActionOfSL2`.
  for (gen, exp) in Iterators.reverse(syllables(s_t_decomposition(A)))
    o = _action_sl2_syllable(o, gen, exp)
  end
  return o
end

# One syllable `gen^exp`. `T^e` sends `(h, v)` to `(h, h^{-e} * v)` and `S` has
# order 4, which is the same action as one letter of the syllable at a time.
function _action_sl2_syllable(o::Origami, gen::Int, exp::ZZRingElem)
  if gen == 2
    iszero(exp) && return o
    h = horizontal_perm(o)
    v = vertical_perm(o)
    return origami(h, _origami_perm_pow(h, -exp) * v; check=false)
  end
  @req gen == 1 "word must be a word in two generators"
  r = mod(exp, 4)
  if iszero(r)
    return o
  elseif isone(r)
    return action_s(o)
  elseif r == 2
    return action_s(action_s(o))
  else
    return action_s_inv(o)
  end
end

# `^(::GAPGroupElem, ::Int)` only accepts GAP immediate integers, and the
# generic `ZZRingElem` power needs `mul!`, which permutations do not support.
function _origami_perm_pow(p::PermGroupElem, e::ZZRingElem)
  if fits(Int, e) && abs(e) < (ZZ(1) << 60)
    return p^Int(e)
  end
  return group_element(parent(p), (GapObj(p)^GapObj(BigInt(e)))::GapObj)
end
