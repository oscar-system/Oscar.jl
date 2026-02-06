#################################################
# Linear Algebraic Groups
#################################################

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

function linear_algebraic_group(type::Symbol, n::Int, k::Field)
  return linear_algebraic_group(root_system(type, n), k)
end

#This is here because of https://github.com/oscar-system/Oscar.jl/issues/5661
function gap_likes_the_group(LAG::LinearAlgebraicGroup)
  return is_finite(LAG.k)
end

function has_gens(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return has_gens(LAG.G)
  return false
end

function number_of_generators(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return number_of_generators(LAG.G)
  error("Group is not finitely generated") # as long as field is QQ
end

function gens(LAG::LinearAlgebraicGroup)
  return [gen(LAG, i) for i in 1:ngens(LAG)]
end

function gen(LAG::LinearAlgebraicGroup, i::Int)
  gap_likes_the_group(LAG) && return linear_algebraic_group_elem(LAG, gen(LAG.G, i))
  error("Group is not finitely generated") # as long as field is QQ
end

function isfinite(LAG::LinearAlgebraicGroup)
  gap_likes_the_group(LAG) && return isfinite(LAG.G)
  if degree(LAG) == 0 || degree(LAG) == 1 #Should not occur
    return true
  else
    return false
  end
end

function order(::Type{T}, LAG::LinearAlgebraicGroup) where {T}
  if !is_finite(LAG)
    throw(InfiniteOrderError(LAG))
  else
    return order(T, LAG.G)
  end # as long as field is QQ no else is needed
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

function one(LAG::LinearAlgebraicGroup)
  return linear_algebraic_group_elem(LAG, one(LAG.G))
end

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
function _compute_action(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
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
end

function _genrating_set_of_unit_group(k::Field)
  gs = FieldElem[]
  u, f = unit_group(k)
  for i in gens(u)
    push!(gs, f(i))
  end
  return gs
end

function root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  if isdefined(LAG, :U_alphas)
    if haskey(LAG.U_alphas, alpha)
      return LAG.U_alphas[alpha]
    end
  else
    LAG.U_alphas = Dict{WeightLatticeElem,MatGroup}()
  end
  @req is_root(alpha) "The given element is not a root"
  @req root_system(alpha) === root_system(LAG) "parent mismatch"
  G = LAG.G
  i, j = _compute_action(LAG, alpha)
  I = identity_matrix(LAG.k, degree(LAG))
  m = zero_matrix(LAG.k, degree(LAG), degree(LAG))
  m[i, j] = one(LAG.k)
  gs = [G(I + lambda * m) for lambda in _genrating_set_of_unit_group(LAG.k)]
  U, _ = sub(G, gs)
  LAG.U_alphas[alpha] = U
  return U
end

########### Tori ############################
function maximal_torus(LAG::LinearAlgebraicGroup)
  isdefined(LAG, :T) && return LAG.T
  G = LAG.G
  gs = MatGroupElem[]
  for t in _genrating_set_of_unit_group(LAG.k)
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

function torus_element(LAG::LinearAlgebraicGroup, diag::Vector{T}) where {T<:FieldElem}
  @req length(diag) == degree(LAG) "Wrong number of diagonal entries"
  m = diagonal_matrix(LAG.k, diag)
  @req det(m) == one(LAG.k) "Deteminant of torus element must be 1"
  return linear_algebraic_group_elem(LAG, MatGroupElem(LAG.G, m))
end

function apply_root_to_torus_element(
  alpha::RootSpaceElem, t::LinearAlgebraicGroupElem
)
  @req is_root(alpha) "The given element is not a root"
  i, j = _compute_action(parent(t), alpha)
  return t.mat[i, i] * inv(t.mat[j, j])
end

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

function borel(LAG::LinearAlgebraicGroup)
  T = maximal_torus(LAG)
  G = LAG.G
  gs = MatGroupElem[]
  for t in gens(T)
    push!(gs, t)
  end
  n = degree(LAG)
  for i in 1:(n - 1)
    for j in (i + 1):n
      for lambda in LAG.k
        m = identity_matrix(LAG.k, n)
        m[i, j] = lambda
        push!(gs, G(m))
      end
    end
  end
  B, _ = sub(G, gs)
  return B
end

function bruhat_cell(LAG::LinearAlgebraicGroup, w::WeylGroupElem)
  @req parent(w) == weyl_group(root_system(LAG)) "parent mismatch"
  B = borel(LAG)
  rep = identity_matrix(LAG.k, degree(LAG))
  for i in word(w)
    alpha = simple_root(root_system(LAG), Int64(i))
    rep = rep * representative_of_root_in_group(LAG, alpha).mat
  end
  return B * LAG.G(rep) * B
end

function bruhat_decomp(LAG::LinearAlgebraicGroup)
  return [bruhat_cell(LAG, w) for w in weyl_group(root_system(LAG))]
end
