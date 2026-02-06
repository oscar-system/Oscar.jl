#################################################
# Linear Algebraic Groups
#################################################

@doc raw"""
    linear_algebraic_group(rs::RootSystem, k::Field) -> LinearAlgebraicGroup

Construct the linear algerbaic group of the given type.

Only type :A is implemented so far.

# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> rs = root_system(:A, 3)
Root system of rank 3
  of type A3

julia> LAG = linear_algebraic_group(rs, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)
```
"""
function linear_algebraic_group(rs::RootSystem, k::Field)
  if !is_finite(k)
    error("Currently only finite fields are supported.")
  end
  rst = root_system_type(rs)
  if length(rst) == 1 && rst[1][1] == :A
    n = rank(rs)
    G = special_linear_group(n + 1, k)
    LAG = LinearAlgebraicGroup(rs, G, k)
  else
    error("Only type A is implemented so far.")
  end
  return LAG
end

@doc raw"""
    linear_algebraic_group(rs::RootSystem, k::Field) -> LinearAlgebraicGroup

Construct the linear algerbaic group of the given type.

Only type :A is implemented so far.

# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)
```
"""
function linear_algebraic_group(type::Symbol, n::Int, k::Field)
  return linear_algebraic_group(root_system(type, n), k)
end

#This is here because of https://github.com/oscar-system/Oscar.jl/issues/5661
function gap_likes_the_group(LAG::LinearAlgebraicGroup)
  return is_finite(LAG.k)
end

"""
    has_gens(LAG::LinearAlgebraicGroup)

Return whether generators for the group `LAG` are known.
"""
function has_gens(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return has_gens(LAG.G)
  return false
end

@doc raw"""
    number_of_generators(LAG::LinearAlgebraicGroup) -> Int

Return the number of generators of `LAG`. Errors if `LAG` is not finitely generated.
"""
function number_of_generators(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return number_of_generators(LAG.G)
  error("Group is not finitely generated") # as long as field is QQ
end

@doc raw"""
    gens(LAG::LinearAlgebraicGroup) -> Vector{LinearAlgebraicGroupElem}

Return the generators of `LAG`. Errors if `LAG` is not finitely generated.
"""
function gens(LAG::LinearAlgebraicGroup)
  return [gen(LAG, i) for i in 1:ngens(LAG)]
end

@doc raw"""
    gen(LAG::LinearAlgebraicGroup, i::Int) -> Int

Return the i-th generator of `LAG`. Errors if `LAG` is not finitely generated.
"""
function gen(LAG::LinearAlgebraicGroup, i::Int)
  gap_likes_the_group(LAG) && return linear_algebraic_group_elem(LAG, gen(LAG.G, i))
  error("Group is not finitely generated") # as long as field is QQ
end

@doc raw"""
    is_finite(LAG::LinearAlgebraicGroup) -> Bool

Return whether `LAG` is finite.
"""
function isfinite(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return isfinite(LAG.G)
  if degree(LAG) == 0 || degree(LAG) == 1 #Should not occur
    return true
  else
    return false
  end
end

@doc raw"""
    order(LAG::LinearAlgebraicGroup) -> Int

Return the order of `LAG`. Errors if `LAG` is not finite.
"""
function order(::Type{T}, LAG::LinearAlgebraicGroup) where {T}
  if !is_finite(LAG)
    throw(InfiniteOrderError(LAG))
  else
    return order(T, LAG.G)
  end
end

function Base.rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{LinearAlgebraicGroup})
  LAG = rs[]
  if gap_likes_the_group(LAG)
    return linear_algebraic_group_elem(LAG, rand(LAG.G))
  elseif LAG.k == QQ #pseudo random matrices in SLn(QQ)
    M = matrix_space(QQ, degree(LAG), degree(LAG))
    a = rand(M, 0:(degree(LAG) * 10))
    while det(a) == 0
      a = rand(M, 0:3)
    end
    return linear_algebraic_group_elem(LAG, multiply_row(a, inv(det(a)), rand(1:6)))
  else
    error("Random elements not implemented for this field yet.")
  end
end

function Base.eltype(::Type{LinearAlgebraicGroup})
  return LinearAlgebraicGroupElem
end

function elem_type(::Type{LinearAlgebraicGroup})
  return LinearAlgebraicGroupElem
end

@doc raw"""
    one(LAG::LinearAlgebraicGroup) -> LinearAlgebraicGroupElem

Return the identity element of `LAG`.
"""
function one(LAG::LinearAlgebraicGroup)
  return linear_algebraic_group_elem(LAG, one(LAG.G))
end

@doc raw"""
    is_subgroup(U::MatGroup, LAG::LinearAlgebraicGroup) -> Bool

Return whether the matrix group `U` is a subgroup of `LAG` as a matrix group.
"""
function is_subgroup(U::MatGroup, LAG::LinearAlgebraicGroup)
  return is_subgroup(U, LAG.G)
end

#################################################
# Linear Algebraic Group Elements
#################################################

# function linear_algebraic_group_elem(LAG::LinearAlgebraicGroup, mat::MatElem)
#   #add checks here
#   return LinearAlgebraicGroupElem(LAG, MatGroupElem(LAG.G, mat))
# end

function linear_algebraic_group_elem(LAG::LinearAlgebraicGroup, MGE::MatGroupElem)
  #TODO: add checks here
  return LinearAlgebraicGroupElem(LAG, MGE)
end

function Base.:(==)(a::LinearAlgebraicGroupElem, b::LinearAlgebraicGroupElem)
  check_perent(a, b)
  return a.mat == b.mat
end

function Base.:(*)(a::LinearAlgebraicGroupElem, b::LinearAlgebraicGroupElem)
  check_perent(a, b)
  return linear_algebraic_group_elem(parent(a), a.mat * b.mat)
end

function Base.inv(a::LinearAlgebraicGroupElem)
  return linear_algebraic_group_elem(parent(a), inv(a.mat))
end

function Base.deepcopy_internal(a::LinearAlgebraicGroupElem, dict::IdDict)
  return get!(dict, a) do
    linear_algebraic_group_elem(parent(a), Base.deepcopy_internal(a.mat, dict))
  end
end

function order(::Type{T}, a::LinearAlgebraicGroupElem) where {T}
  return order(T, a.mat)
end

function Base.hash(a::LinearAlgebraicGroupElem, h::UInt)
  b = 0x1df4d55a7b37db2f % UInt
  h = hash(parent(a), h)
  h = hash(a.mat, h)

  return xor(h, b)
end

############# Root Subgroups ############################
#internal function to compute action of root alpha, in case :A return the tuple (i, j) for which alpha acts like e_i-e_j
function _compute_action(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  if root_system_type(root_system(LAG))[1][1] == :A
    c = coefficients(alpha)
    l = number_of_simple_roots(root_system(LAG))
    e = zeros(Int64, l + 1)
    for i in 1:l
      e[i] = e[i] + Int64(c[i])
      e[i + 1] = e[i + 1] - Int64(c[i])
    end
    i = 0
    j = 0
    for k in 1:(l + 1)
      if e[k] == 1
        i = k
      elseif e[k] == -1
        j = k
      end
    end
    return i, j
  else
    error("Only type A is implemented so far.")
  end
end

#internal function to get generators of the unit group of the field
function _generating_set_of_unit_group(k::Field)
  gs = FieldElem[]
  u, f = unit_group(k)
  for i in gens(u)
    push!(gs, f(i))
  end
  return gs
end

@doc raw"""
    root_subgroup_generator(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem) -> Matrix

Return the Matrix that generates the root subgroup of `LAG` corresponding to the root `alpha`.


# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> alpha = simple_root(root_system(LAG),2)
a_2

julia> root_subgroup_generator(LAG, alpha)
[0   0   0   0]
[0   0   1   0]
[0   0   0   0]
[0   0   0   0]
```
"""
function root_subgroup_generator(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  @req is_root(alpha) "The given element is not a root"
  @req root_system(alpha) === root_system(LAG) "parent mismatch"
  i, j = _compute_action(LAG, alpha)
  m = zero_matrix(LAG.k, degree(LAG), degree(LAG))
  m[i, j] = one(LAG.k)
  return m
end

@doc raw"""
    root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem) -> MatGroup

Return the root subgroup of `LAG` corresponding to the root `alpha`.


# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> alpha = simple_root(root_system(LAG),2)
a_2

julia> root_subgroup(LAG, alpha)
Matrix group of degree 4
  over prime field of characteristic 5
```
"""
function root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  if isdefined(LAG, :U_alphas)
    if haskey(LAG.U_alphas, alpha)
      return LAG.U_alphas[alpha]
    end
  else
    LAG.U_alphas = Dict{WeightLatticeElem,MatGroup}()
  end
  G = LAG.G
  I = identity_matrix(LAG.k, degree(LAG))
  m = root_subgroup_generator(LAG, alpha)
  gs = [G(I + lambda * m) for lambda in basis(LAG.k)]
  U, _ = sub(G, gs)
  LAG.U_alphas[alpha] = U
  return U
end

########### Tori ############################
@doc raw"""
    maximal_torus(LAG::LinearAlgebraicGroup) -> MatGroup

Return the standard maximal torus of diagonal elements of `LAG`.


# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> maximal_torus(LAG)
Matrix group of degree 4
  over prime field of characteristic 5
```
"""
function maximal_torus(LAG::LinearAlgebraicGroup)
  isdefined(LAG, :T) && return LAG.T
  G = LAG.G
  gs = MatGroupElem[]
  for t in _generating_set_of_unit_group(LAG.k)
    it = inv(t)
    for i in 1:(degree(LAG) - 1)
      m = identity_matrix(LAG.k, degree(LAG))
      m[i, i] = t
      m[i + 1, i + 1] = it
      push!(gs, G(m))
    end
  end
  T, _ = sub(G, gs)
  LAG.T = T
  return T
end

@doc raw"""
  torus_element(LAG::LinearAlgebraicGroup, diag::Vector{T}) where {T<:FieldElem} -> LinearAlgebraicGroupElem

Return the root element consisting of `diag` as diagonal entries which is an element of the maximal torus of `LAG`.


# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> torus_element(LAG, [F(1),F(2),F(1),F(3)])
LinearAlgebraicGroupElem(LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef), [1 0 0 0; 0 2 0 0; 0 0 1 0; 0 0 0 3], #undef)
```
"""
function torus_element(LAG::LinearAlgebraicGroup, diag::Vector{T}) where {T<:FieldElem}
  @req length(diag) == degree(LAG) "Wrong number of diagonal entries"
  m = diagonal_matrix(LAG.k, diag)
  @req det(m) == one(LAG.k) "Deteminant of torus element must be 1"
  return linear_algebraic_group_elem(LAG, MatGroupElem(LAG.G, m))
end

@doc raw"""
  apply_root_to_torus_element(alpha::RootSpaceElem, t::LinearAlgebraicGroupElem) -> FieldElem

Return the field element obtained by applying the root `alpha` to the torus element `t`.

# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> alpha = simple_root(root_system(LAG),2)
a_2

julia> t = torus_element(LAG, [F(1),F(2),F(1),F(3)])
LinearAlgebraicGroupElem(LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef), [1 0 0 0; 0 2 0 0; 0 0 1 0; 0 0 0 3], #undef)

julia> apply_root_to_torus_element(alpha, t)
2
```
"""
function apply_root_to_torus_element(
  alpha::RootSpaceElem, t::LinearAlgebraicGroupElem
)
  @req is_root(alpha) "The given element is not a root"
  @req in(t.mat, maximal_torus(parent(t))) "The given element is not a torus element"
  i, j = _compute_action(parent(t), alpha)
  return t.mat[i, i] * inv(t.mat[j, j])
end

@doc raw"""
  representative_of_root_in_group(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem) -> LinearAlgebraicGroupElem

Return the linear algerbaic group element corresponding to the root `alpha`.

# Examples
```jldoctest
julia> F, _  = finite_field(4)
(Finite field of degree 2 and characteristic 2, o)

julia> LAG = linear_algebraic_group(:A, 4, F)
LinearAlgebraicGroup(Root system of type A4, SL(5,4), Finite field of degree 2 and characteristic 2, #undef, #undef, #undef)

julia> alpha = simple_root(root_system(LAG),2)
a_2

julia> representative_of_root_in_group(LAG, alpha)
LinearAlgebraicGroupElem(LinearAlgebraicGroup(Root system of type A4, SL(5,4), Finite field of degree 2 and characteristic 2, #undef, #undef, #undef), [1 0 0 0 0; 0 0 1 0 0; 0 1 0 0 0; 0 0 0 1 0; 0 0 0 0 1], #undef)
```
"""
function representative_of_root_in_group(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  @req is_root(alpha) "The given element is not a root"
  i, j = _compute_action(LAG, alpha)
  m = identity_matrix(LAG.k, degree(LAG))
  m[i, i] = zero(LAG.k)
  m[i, j] = -one(LAG.k)
  m[j, i] = one(LAG.k)
  m[j, j] = zero(LAG.k)
  return linear_algebraic_group_elem(LAG, MatGroupElem(LAG.G, m))
end

@doc raw"""
  borel(LAG::LinearAlgebraicGroup) -> MatGroup

Return the standard borel subgroup of the linear algerbaic group `LAG`.

# Examples
```jldoctest
julia> F, _  = finite_field(3)
(Prime field of characteristic 3, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,3), Prime field of characteristic 3, #undef, #undef, #undef)

julia> borel(LAG)
Matrix group of degree 4
  over prime field of characteristic 3
```
"""
function borel(LAG::LinearAlgebraicGroup)
  T = maximal_torus(LAG)
  G = LAG.G
  gs = MatGroupElem[]
  for t in gens(T)
    push!(gs, t)
  end
  for lambda in basis(LAG.k)
    for alpha in simple_roots(root_system(LAG))
      push!(gs, MatGroupElem(G, lambda * root_subgroup_generator(LAG, alpha)))
    end
  end
  B, _ = sub(G, gs)
  return B
end

@doc raw"""
  bruhat_cell_rep(LAG::LinearAlgebraicGroup, w::WeylGroupElem) -> MatGroupElem

Return the representative of the bruhat cell corresponding to the weyl group element `w`.

# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> W = weyl_group(root_system(LAG))
Weyl group
  of root system of rank 3
    of type A3

julia> w = W([1,2])
s1 * s2

julia> bruhat_cell_rep(LAG, w)
[0   0   1   0]
[1   0   0   0]
[0   1   0   0]
[0   0   0   1]
```
"""
function bruhat_cell_rep(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
  @req parent(w) == weyl_group(root_system(LAG)) "parent mismatch"
  rep = identity_matrix(LAG.k, degree(LAG))
  for i in word(w)
    alpha = simple_root(root_system(LAG), Int64(i))
    rep = rep * representative_of_root_in_group(LAG, alpha).mat
  end
  return LAG.G(rep)
end

@doc raw"""
  bruhat_cell(LAG::LinearAlgebraicGroup, w::WeylGroupElem) -> GroupDoubleCoset{MatGroup, MatGroupElem}

Return the bruhat cell corresponding to the weyl group element `w`.

# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> W = weyl_group(root_system(LAG))
Weyl group
  of root system of rank 3
    of type A3

julia> w = W([1,2])
s1 * s2

julia> bruhat_cell(LAG,w)
Double coset of matrix group of degree 4 over F
  and matrix group of degree 4 over F
  with representative [0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1]
  in SL(4,5)
```
"""
function bruhat_cell(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
  rep = bruhat_cell_rep(LAG, w)
  B = borel(LAG)
  return double_coset(B, rep, B)
end

@doc raw"""
  bruhat_decomp(LAG::LinearAlgebraicGroup) -> Vector{GroupDoubleCoset{MatGroup, MatGroupElem}}

Return the bruhat decomposition of the linear algebraic group `LAG`.

# Examples
```jldoctest
julia> F, _  = finite_field(5)
(Prime field of characteristic 5, 0)

julia> LAG = linear_algebraic_group(:A, 3, F)
LinearAlgebraicGroup(Root system of type A3, SL(4,5), Prime field of characteristic 5, #undef, #undef, #undef)

julia> bruhat_decomp(LAG)
24-element Vector{GroupDoubleCoset{MatGroup{FqFieldElem, FqMatrix}, MatGroupElem{FqFieldElem, FqMatrix}}}:
 Double coset of matrix group and matrix group with representative [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
 Double coset of matrix group and matrix group with representative [0 4 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]
 Double coset of matrix group and matrix group with representative [0 4 0 0; 0 0 4 0; 1 0 0 0; 0 0 0 1]
 Double coset of matrix group and matrix group with representative [0 0 1 0; 0 4 0 0; 1 0 0 0; 0 0 0 1]
 Double coset of matrix group and matrix group with representative [0 4 0 0; 0 0 4 0; 0 0 0 4; 1 0 0 0]
 Double coset of matrix group and matrix group with representative [0 0 1 0; 0 4 0 0; 0 0 0 4; 1 0 0 0]
 Double coset of matrix group and matrix group with representative [0 0 1 0; 0 0 0 1; 0 4 0 0; 1 0 0 0]
 Double coset of matrix group and matrix group with representative [0 0 0 4; 0 0 1 0; 0 4 0 0; 1 0 0 0]
 Double coset of matrix group and matrix group with representative [0 4 0 0; 0 0 0 1; 0 0 4 0; 1 0 0 0]
 Double coset of matrix group and matrix group with representative [0 0 0 4; 0 4 0 0; 0 0 4 0; 1 0 0 0]
 ⋮
 Double coset of matrix group and matrix group with representative [0 0 0 4; 0 0 1 0; 1 0 0 0; 0 1 0 0]
 Double coset of matrix group and matrix group with representative [1 0 0 0; 0 0 0 1; 0 0 4 0; 0 1 0 0]
 Double coset of matrix group and matrix group with representative [0 0 0 4; 1 0 0 0; 0 0 4 0; 0 1 0 0]
 Double coset of matrix group and matrix group with representative [1 0 0 0; 0 1 0 0; 0 0 0 4; 0 0 1 0]
 Double coset of matrix group and matrix group with representative [0 4 0 0; 1 0 0 0; 0 0 0 4; 0 0 1 0]
 Double coset of matrix group and matrix group with representative [0 4 0 0; 0 0 0 1; 1 0 0 0; 0 0 1 0]
 Double coset of matrix group and matrix group with representative [0 0 0 4; 0 4 0 0; 1 0 0 0; 0 0 1 0]
 Double coset of matrix group and matrix group with representative [1 0 0 0; 0 0 0 1; 0 1 0 0; 0 0 1 0]
 Double coset of matrix group and matrix group with representative [0 0 0 4; 1 0 0 0; 0 1 0 0; 0 0 1 0]
```
"""
function bruhat_decomp(LAG::LinearAlgebraicGroup)
  return [bruhat_cell(LAG, w) for w in weyl_group(root_system(LAG))]
end
