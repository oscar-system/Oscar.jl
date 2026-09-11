include("exports.jl")
include("types.jl")

# Degree of the coset action, i.e. the index of the subgroup in SL(2,Z).
# The trivial action (index 1) moves no points.
function _index(s::PermGroupElem, t::PermGroupElem)
  return max(1, largest_moved_point(s), largest_moved_point(t))
end

function modular_subgroup(s::PermGroupElem, t::PermGroupElem, check::Bool)
  # Coerce both permutations into the same symmetric group, so that
  # the caller may pass permutations of different degrees.
  Sym = symmetric_group(_index(s, t))
  s = Sym(s)
  t = Sym(t)
  return ModularGroup(s, t, check)
end

@doc raw"""
    modular_subgroup_via_right_action(s::PermGroupElem, t::PermGroupElem; check=true)

Construct a `ModularGroup` object corresponding to the finite-index subgroup
of ``{\rm SL}_2(\mathbb{Z})`` described by the permutations ``s`` and ``t``. For `check = true`, this
constructor tests if the given permutations actually describe the (right) coset action of the matrices
```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}, \qquad T= \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix},
```
by checking that they act transitively and satisfy the relations
``s^4 = (s^3 t)^3 = s^2 t s^{-2} t^{-1} = 1``.
The point ``1`` corresponds to the coset of the identity matrix.

# Examples
```jldoctest
julia> s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
(1,2)(3,4)(5,6)(7,8)(9,10)

julia> t = cperm([1,4], [2,5,9,10,8], [3,7,6])
(1,4)(2,5,9,10,8)(3,7,6)

julia> G = modular_subgroup_via_right_action(s, t)
Modular subgroup of index 10
```
"""
function modular_subgroup_via_right_action(s::PermGroupElem, t::PermGroupElem; check=true)
  return modular_subgroup(s, t, check)
end

@doc raw"""
    modular_subgroup_via_left_action(s::PermGroupElem, t::PermGroupElem; check=true)

Same as [`modular_subgroup_via_right_action`](@ref), but now the permutations describe the action
by left multiplication on the left cosets. Under the bijection ``gH \mapsto Hg^{-1}``
this is the right action of the inverse matrices, hence the subgroup is stored via
the inverse permutations.

# Examples
```jldoctest
julia> s = cperm([1,2], [3,4], [5,6], [7,8], [9,10])
(1,2)(3,4)(5,6)(7,8)(9,10)

julia> t = cperm([1,4], [2,5,9,10,8], [3,7,6])
(1,4)(2,5,9,10,8)(3,7,6)

julia> G = modular_subgroup_via_left_action(s, t)
Modular subgroup of index 10
```
"""
function modular_subgroup_via_left_action(s::PermGroupElem, t::PermGroupElem; check=true)
  return modular_subgroup(s^-1, t^-1, check)
end

function Base.:(==)(G::ModularGroup, H::ModularGroup)
  return index(G) == index(H) && issubset(G, H)
end

@attr NTuple{3, CycleType} function _cycle_structures(G::ModularGroup)
  return (cycle_structure(G.s), cycle_structure(G.t), cycle_structure(G.s * G.t))
end

function Base.hash(G::ModularGroup, h::UInt)
  s_struc, t_struc, s_t_struc = _cycle_structures(G)
  return hash((index(G), s_struc, t_struc, s_t_struc), h)
end

function Base.show(io::IO, G::ModularGroup)
  idx = index(G)
  print(io, "Modular subgroup of index $(idx)")
end

function Base.show(io::IO, x::ModularGroupElem)
  A = matrix(x)
  idx = index(parent(x))
  print(io, "$(A) in moodular group of index $(idx)")
end

elem_type(::Type{ModularGroup}) = ModularGroupElem
parent_type(::Type{ModularGroupElem}) = ModularGroup

parent(x::ModularGroupElem) = x.parent
matrix(x::ModularGroupElem) = x.mat

one(G::ModularGroup) = ModularGroupElem(G, identity_matrix(ZZ, 2))
one(x::ModularGroupElem) = one(parent(x))

function order(::Type{T}, x::ModularGroupElem) where {T<:IntegerUnion}
  A = matrix(x)
  t = A[1, 1] + A[2, 2]

  if t == 2
    isone(A) || throw(InfiniteOrderError(x))
    return T(1)
  elseif t == -2
    isone(-A) || throw(InfiniteOrderError(x))
    return T(2)
  elseif t == -1
    return T(3)
  elseif t == 0
    return T(4)
  elseif t == 1
    return T(6)
  else
    throw(InfiniteOrderError(x))
  end
end

order(x::ModularGroupElem) = order(ZZRingElem, x)

function is_finite_order(x::ModularGroupElem)
  A = matrix(x)
  t = tr(A)
  return abs(t) < 2 || (t == 2 && isone(A)) || (t == -2 && A == -one(A))
end

@doc raw"""
    s_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}``
on the right cosets of `G`.
"""
function s_right_perm(G::ModularGroup)
  return G.s
end

@doc raw"""
    t_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``T= \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
function t_right_perm(G::ModularGroup)
  return G.t
end

@doc raw"""
    r_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``R= \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
@attr PermGroupElem function r_right_perm(G::ModularGroup)
  return G.s^-1*G.t^-1*G.s
end

@doc raw"""
    j_right_perm(G::ModularGroup)

Return the permutation describing the action of the matrix ``J= \begin{pmatrix} 0 & 1 \\ -1 & 1 \end{pmatrix}``
on the right cosets of `G`.
"""
@attr PermGroupElem function j_right_perm(G::ModularGroup)
  return G.s^-1*G.t^-1
end

function defines_coset_action_s_t(s::PermGroupElem, t::PermGroupElem)
  isone(s^4) || return false
  isone((s^3*t)^3) || return false
  isone(s^2*t*s^-2*t^-1) || return false
  return is_transitive(permutation_group(_index(s, t), [s, t]))
end

@attr Int function index(G::ModularGroup)
  return _index(s_right_perm(G), t_right_perm(G))
end

const _SL2Z_FP_CACHE = Ref{FPGroup}()
const _MATRIX_HOM_CACHE = Ref{GAPGroupHomomorphism{FPGroup, MatGroup{ZZRingElem, ZZMatrix}}}()

# cache the SL2Z presentation so that it can be reused consistently
# (especially in testing)
function _SL2Z_fp()
  if !isassigned(_SL2Z_FP_CACHE)
    # Words like T^k with huge k arise from matrices with large entries;
    # the syllable representation stores them in constant space.
    F = free_group(["S", "T"]; eltype = :syllable)
    S, T = gens(F)

    SL2Z, _ = quo(F, [S^4, (S^3*T)^3, S^2*T*S^-2*T^-1])

    _SL2Z_FP_CACHE[] = SL2Z
  end

  return _SL2Z_FP_CACHE[]::FPGroup
end

_matrix_S() = matrix(ZZ, [0 -1; 1 0])
_matrix_T() = matrix(ZZ, [1 1; 0 1])

# cache the homomorphism from the fp presentation of SL2Z to the matrix group
function _matrix_hom()
  if !isassigned(_MATRIX_HOM_CACHE)
    M = matrix_group([_matrix_S(), _matrix_T()])
    MS, MT = gens(M)

    SL2Z = _SL2Z_fp()
    _MATRIX_HOM_CACHE[] = hom(SL2Z, M, [MS, MT])
  end

  return _MATRIX_HOM_CACHE[]::GAPGroupHomomorphism{FPGroup, MatGroup{ZZRingElem, ZZMatrix}}
end

# cache the permutation group generated by s and t
@attr PermGroup function _perm_group(G::ModularGroup)
  return permutation_group(index(G), [s_right_perm(G), t_right_perm(G)])
end

# cache the homomorphism from SL2Z to P
@attr GAPGroupHomomorphism{FPGroup, PermGroup} function _coset_action_hom(G::ModularGroup)
  SL2Z = _SL2Z_fp()
  P = _perm_group(G)
  return hom(SL2Z, P, [s_right_perm(G), t_right_perm(G)])
end

inv(x::ModularGroupElem) = ModularGroupElem(parent(x), inv(matrix(x)))

*(x::ModularGroupElem, y::ModularGroupElem) = ModularGroupElem(parent(x), matrix(x) * matrix(y))

^(x::ModularGroupElem, n::Integer) =
  ModularGroupElem(parent(x), matrix(x)^n)

==(x::ModularGroupElem, y::ModularGroupElem) =
  parent(x) === parent(y) && matrix(x) == matrix(y)

Base.hash(x::ModularGroupElem, h::UInt) = hash(matrix(x), h)

@doc raw"""
    word_gens(G::ModularGroup)

Return generators of `G` as words in the generators `S` and `T` of the
finitely presented group returned by `Oscar._SL2Z_fp()`.
"""
function word_gens(G::ModularGroup)
  P = _perm_group(G)
  phi = _coset_action_hom(G)

  # TODO: might need to check whether this is efficient enough for large index
  Hperm, _ = stabilizer(P, 1)

  H, inc = preimage(phi, Hperm)

  return [inc(h) for h in gens(H)]
end

@attr Vector{ModularGroupElem} function gens(G::ModularGroup)
  w_gens = word_gens(G)
  phi = _matrix_hom()
  return [ModularGroupElem(G, matrix(phi(w))) for w in w_gens]
end

function gen(G::ModularGroup, i::Int)
  gs = gens(G)
  @req abs(i) <= length(gs) "i must be in the range -$(length(gs)):$(length(gs))"
  i > 0 && return gs[i]
  i < 0 && return inv(gs[-i])
  return one(G)
end

function number_of_generators(G::ModularGroup)
  return length(gens(G))
end

@doc raw"""
    s_t_decomposition(M::ZZMatrix)

Return the matrix `M` in ``{\rm SL}_2(\mathbb{Z})`` as a word in the generators `S` and `T`
of the finitely presented group returned by `Oscar._SL2Z_fp()`.

# Examples
```jldoctest
julia> s_t_decomposition(matrix(ZZ, [1 0; -4 1]))
S*T^4*S^-1
```
"""
function s_t_decomposition(M::ZZMatrix)
  @req size(M) == (2, 2) "Matrix needs to be in SL(2, Z)"
  @req isone(det(M)) "Matrix needs to be in SL(2, Z)"

  MatS = _matrix_S()
  MatT = _matrix_T()
  SL2Z = _SL2Z_fp()
  S, T = gens(SL2Z)

  decomp = one(SL2Z)

  # Euclidean algorithm on the second row: M * T^-k * S has second row
  # (d - k*c, -c), so the lower left entry decreases in absolute value.
  while !iszero(M[2, 1])
    k = div(M[2, 2], M[2, 1])
    decomp = S^-1 * T^k * decomp
    M = M * MatT^-k * MatS
  end

  # now M[2, 1] = 0 and since det(M) == 1, we have M == +-T^r where r = +-M[1, 2]
  if isone(M[1, 1])
    decomp = T^M[1, 2] * decomp
  else
    decomp = S^2 * T^-M[1, 2] * decomp
  end

  return decomp
end

function coset_action_of(A::ZZMatrix, G::ModularGroup)
  w = s_t_decomposition(A)
  phi = _coset_action_hom(G)
  return phi(w)
end

@doc raw"""
    coset_right_action_of(A::ZZMatrix, G::ModularGroup)

Return the permutation describing the action of the matrix `A` in ``{\rm SL}_2(\mathbb{Z})``
on the right cosets of `G` by right multiplication.
"""
function coset_right_action_of(A::ZZMatrix, G::ModularGroup)
  return coset_action_of(A, G)
end

@doc raw"""
    coset_left_action_of(A::ZZMatrix, G::ModularGroup)

Return the permutation describing the action of the matrix `A` in ``{\rm SL}_2(\mathbb{Z})``
on the left cosets of `G` by left multiplication, see [`modular_subgroup_via_left_action`](@ref).
"""
function coset_left_action_of(A::ZZMatrix, G::ModularGroup)
  return coset_action_of(A, G)^-1
end

function _image_of_pt(w::FPGroupElem, G::ModularGroup, pt::Integer)
  @req 1 <= pt <= index(G) "Non-valid point"

  imgs = (s_right_perm(G), t_right_perm(G))

  for (i, e) in syllables(w) # i = 1 or 2, e = exponent
    pt = pt^(imgs[i]^e)
  end
  return pt
end

#TODO: could be optimized by tracing the image of the point 1 directly in the Euclidean loop,
# instead of decomposing A and take the resulting word apart again.
function Base.in(A::ZZMatrix, G::ModularGroup)
  (size(A) == (2, 2) && isone(det(A))) || return false
  w = s_t_decomposition(A)
  return _image_of_pt(w, G, 1) == 1
end

function Base.in(x::ModularGroupElem, G::ModularGroup)
  return parent(x) === G || matrix(x) in G
end

@doc raw"""
    is_word_element_of(w::FPGroupElem, G::ModularGroup)

Return whether the word `w` in the generators `S` and `T` of the finitely
presented group returned by `Oscar._SL2Z_fp()` represents an element of `G`.
"""
function is_word_element_of(w::FPGroupElem, G::ModularGroup)
  return _image_of_pt(w, G, 1) == 1
end

function Base.issubset(H::ModularGroup, G::ModularGroup)
  nH = index(H)
  nG = index(G)

  # Necessary numerical conditions: [SL2Z : H] = [SL2Z : G] * [G : H]
  nH >= nG || return false
  iszero(mod(nH, nG)) || return false

  sH = s_right_perm(H)
  tH = t_right_perm(H)
  sG = s_right_perm(G)
  tG = t_right_perm(G)

  f = zeros(Int, nH) # f[x] = image of H-coset x in the G-cosets
  f[1] = 1 # identity coset |-> identity coset

  queue = Vector{Int}(undef, nH)
  queue[1] = 1
  head = 1
  tail = 1

  @inbounds while head <= tail
    x  = queue[head]
    head += 1
    fx = f[x]

    for (pH, pG) in ((sH, sG), (tH, tG))
      y  = x^pH  # neighbour in the H-coset graph
      fy = fx^pG # forced image under equivariance

      if fy < 1 || fy > nG
        return false
      elseif f[y] == 0
        f[y] = fy
        tail += 1
        queue[tail] = y
      elseif f[y] != fy
        return false # conflicting assignment
      end
    end
  end

  # Transitivity of the H-action (already guaranteed by the constructor,
  # but cheap to assert that the map was defined everywhere).
  return tail == nH
end

is_abelian(::ModularGroup) = false
is_finite(::ModularGroup) = false
is_cyclic(::ModularGroup) = false
is_trivial(::ModularGroup) = false
is_finitely_generated(::ModularGroup) = true
#TODO random elements?
