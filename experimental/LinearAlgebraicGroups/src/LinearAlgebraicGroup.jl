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
  if type == :A
    R = root_system(type, n)
    G = special_linear_group(n + 1, k)
    LAG = LinearAlgebraicGroup(R, G, k)
  else
    error("Only type A is implemented so far.")
  end
  return LAG
end

#This is here because of https://github.com/oscar-system/Oscar.jl/issues/5661
function gap_likes_the_group(LAG::LinearAlgebraicGroup)
  return is_finite(LAG.k)
end

function has_gens(LAG::LinearAlgebraicGroup)
  if gap_likes_the_group(LAG)
    return has_gens(LAG.G)
  else
    return false
  end
end

function number_of_generators(LAG::LinearAlgebraicGroup)
  if gap_likes_the_group(LAG)
    return number_of_generators(LAG.G)
  else
    error("Generators not known") # as long as field is QQ
  end
end

function gens(LAG::LinearAlgebraicGroup)
  if gap_likes_the_group(LAG)
    return [linear_algebraic_group_elem(LAG, g) for g in gens(LAG.G)]
  else
    error("Generators not known") # as long as field is QQ
  end
end

function gen(LAG::LinearAlgebraicGroup, i::Int)
  if gap_likes_the_group(LAG)
    return linear_algebraic_group_elem(LAG, gen(LAG.G, i))
  else
    error("Generators not known") # as long as field is QQ    
  end
end

function isfinite(LAG::LinearAlgebraicGroup)
  if gap_likes_the_group(LAG)
    return isfinite(LAG.G)
  else
    if degree(LAG) == 0 || degree(LAG) == 1 #Should not occur
      return true
    else
      return false
    end
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

#################################################
# Linear Algebraic Group Elements
#################################################

# function linear_algebraic_group_elem(LAG::LinearAlgebraicGroup, mat::MatElem)
#   #add checks here
#   return LinearAlgebraicGroupElem(LAG, MatGroupElem(LAG.G, mat))
# end

function linear_algebraic_group_elem(LAG::LinearAlgebraicGroup, MGE::MatGroupElem)
  #add checks here
  return LinearAlgebraicGroupElem(LAG, MGE)
end

function Base.:(==)(a::LinearAlgebraicGroupElem, b::LinearAlgebraicGroupElem)
  return parent(a) == parent(b) && a.mat == b.mat
end

function Base.:(*)(a::LinearAlgebraicGroupElem, b::LinearAlgebraicGroupElem)
  if a.parent != b.parent
    throw(ArgumentError("Elements must be from the same group to be multiplied"))
  end
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
function root_subgroup(LAG::LinearAlgebraicGroup, alpha::RootSpaceElem)
  @req is_root(alpha) "The given element is not a root"
  @req root_system(alpha) === root_system(LAG) "parent mismatch"
  G = LAG.G
  c = coefficients(alpha)
  l = number_of_simple_roots(root_system(LAG))
  e = [0,0,0,0]
  for i in 1:l
    e[i] = e[i] + Int64(c[i])
    e[i+1] = e[i+1] -Int64(c[i])
  end
  i=0
  j=0
  for k in 1:l+1
    if e[k] == 1
      i = k
    elseif e[k] == -1
      j = k
    end
  end
  I = identity_matrix(LAG.k,l+1)
  m = zero_matrix(LAG.k,l+1,l+1)
  m[i,j] = one(LAG.k)
  gens = [G(I + lambda * m) for lambda in LAG.k]
  U = sub(G,gens)
  return U
end
