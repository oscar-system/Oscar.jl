###############################################################################
#
#   Root Systems and Weights
#
###############################################################################

mutable struct RootSystem
  cartan_matrix::ZZMatrix # (generalized) Cartan matrix
  #fw::QQMatrix # fundamental weights as linear combination of simple roots
  positive_roots::Any #::Vector{RootSpaceElem} (cyclic reference)
  weyl_group::Any     #::WeylGroup (cyclic reference)
  type::Vector{Tuple{Symbol,Int}}

  function RootSystem(mat::ZZMatrix)
    roots, refl = positive_roots_and_reflections(mat)
    finite = count(refl .== 0) == nrows(mat)

    R = new(mat)
    R.positive_roots = map(r -> RootSpaceElem(R, r), roots)
    R.weyl_group = WeylGroup(finite, refl, R)

    return R
  end
end

@doc raw"""
    root_system(cartan_matrix::ZZMatrix) -> RootSystem
    root_system(cartan_matrix::Matrix{Int}) -> RootSystem

Constructs the root system defined by the Cartan matrix.
"""
function root_system(cartan_matrix::ZZMatrix)
  return RootSystem(cartan_matrix)
end

function root_system(cartan_matrix::Matrix{<:Integer})
  return RootSystem(matrix(ZZ, cartan_matrix))
end

@doc raw"""
    root_system(fam::Symbol, rk::Int) -> RootSystem

Constructs the root system of the given type. See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.
"""
function root_system(fam::Symbol, rk::Int)
  cartan = cartan_matrix(fam, rk)
  R = root_system(cartan)
  R.type = [(fam, rk)]
  return R
end

function root_system(types::Tuple{Symbol,Int}...)
  cartan = cartan_matrix(types...)
  R = root_system(cartan)
  R.type = types
  return R
end

function Base.show(io::IO, R::RootSystem)
  println(io, "Root system defined by Cartan matrix")
  show(io, cartan_matrix(R))
end

@doc raw"""
    cartan_matrix(R::RootSystem) -> ZZMatrix

Returns the Cartan matrix defining `R`.
"""
function cartan_matrix(R::RootSystem)
  return R.cartan_matrix
end

#function basis(R::RootSystem)
#  rk = rank(R)
#  it = 1:rk
#  return [RootSpaceElem(R, matrix(QQ, 1, rk, i .== it)) for i in 1:rk]
#end

function fundamental_weight(R::RootSystem, i::Int)
  @req 1 <= i <= rank(R) "invalid index"
  return WeightLatticeElem(R, matrix(ZZ, rank(R), 1, i .== 1:rank(R)))
end

@doc raw"""
    fundamental_weights(R::RootSystem) -> Vector{WeightLatticeElem}

Returns the fundamental weights corresponding to the simple roots of `R`.
"""
function fundamental_weights(R::RootSystem)
  return [fundamental_weight(R, i) for i in 1:rank(R)]
end

function has_root_system_type(R::RootSystem)
  return isdefined(R, :type)
end

function is_simple(R::RootSystem)
  if has_root_system_type(R)
    return length(root_system_type(R)) == 1
  end
  error("Not implemented") # TODO: implement is_simple
end

function negative_root(R::RootSystem, i::Int)
  return -R.positive_roots[i]::RootSpaceElem
end

function negative_roots(R::RootSystem)
  return [-r for r in positive_roots(R)]
end

function num_positive_roots(R::RootSystem)
  return length(R.positive_roots)
end

function num_roots(R::RootSystem)
  return 2 * length(R.positive_roots)
end

function num_simple_roots(R::RootSystem)
  return rank(R)
end

function nroots(R::RootSystem)
  return num_roots(R)
end

function positive_root(R::RootSystem, i::Int)
  return R.positive_roots[i]::RootSpaceElem
end

function positive_roots(R::RootSystem)
  return R.positive_roots::Vector{RootSpaceElem}
end

@doc raw"""
    rank(R::RootSystem) -> Int

Returns the rank of `R`.
"""
function rank(R::RootSystem)
  return nrows(R.cartan_matrix)
end

function root_system_type(R::RootSystem)
  @req has_root_system_type(R) "root system type not defined"
  return R.type
end

function root_system_type_string(R::RootSystem)
  @req has_root_system_type(R) "root system type not defined"
  return join([string(t[1]) * string(t[2]) for t in R.type], " x ")
end

function root(R::RootSystem, i::Int)
  if i <= num_positive_roots(R)
    return positive_root(R, i)
  else
    return negative_root(R, i - num_positive_roots(R))
  end
end

function roots(R::RootSystem)
  return [[r for r in positive_roots(R)]; [-r for r in positive_roots(R)]]
end

@doc raw"""
    weyl_group(R::RootSystem) -> WeylGroup

Returns the Weyl group of `R`.
"""
function weyl_group(R::RootSystem)
  return R.weyl_group::WeylGroup
end

@doc raw"""
    weyl_vector(R::RootSystem) -> WeightLatticeElem

Returns the Weyl vector $\rho$ of `R`, which is the sum of all fundamental weights,
or half the sum of all positive roots.
"""
function weyl_vector(R::RootSystem)
  return WeightLatticeElem(R, matrix(ZZ, rank(R), 1, fill(1, rank(R))))
end

###############################################################################
# RootSpaceElem

struct RootSpaceElem
  root_system::RootSystem
  vec::QQMatrix # the coordinate (row) vector with respect to the simple roots
end

function RootSpaceElem(root_system::RootSystem, vec::Vector{<:RationalUnion})
  return RootSpaceElem(root_system, matrix(QQ, 1, length(vec), vec))
end

function Base.:(*)(q::RationalUnion, r::RootSpaceElem)
  return RootSpaceElem(root_system(r), q * r.vec)
end

function Base.:(+)(r::RootSpaceElem, r2::RootSpaceElem)
  @req r.root_system === r2.root_system "$r and $r2 must belong to the same root space"

  return RootSpaceElem(r.root_system, r.vec + r2.vec)
end

function Base.:(-)(r::RootSpaceElem, r2::RootSpaceElem)
  @req r.root_system === r2.root_system "$r and $r2 must belong to the same root space"

  return RootSpaceElem(r.root_system, r.vec - r2.vec)
end

function Base.:(-)(r::RootSpaceElem)
  return RootSpaceElem(r.root_system, -r.vec)
end

function Base.:(==)(r::RootSpaceElem, r2::RootSpaceElem)
  return r.root_system === r2.root_system && r.vec == r2.vec
end

function Base.deepcopy_internal(r::RootSpaceElem, dict::IdDict)
  if haskey(dict, r)
    return dict[r]
  end

  w2 = RootSpaceElem(r.root_system, deepcopy_internal(r.vec, dict))
  dict[r] = w2
  return w2
end

@doc raw"""
    getindex(r::RootSpaceElem, i::Int) -> QQRingElem

Returns the coefficient of the `i`-th simple root in `r`.
"""
function Base.getindex(r::RootSpaceElem, i::Int)
  return coeff(r, i)
end

function Base.hash(r::RootSpaceElem, h::UInt)
  b = 0xbe7603eb38c985ad % UInt
  h = hash(r.root_system, h)
  h = hash(r.vec, h)
  return xor(b, h)
end

function coefficients(r::RootSpaceElem)
  return r.vec
end

function coeff(r::RootSpaceElem, i::Int)
  return r.vec[i]
end

function is_root_with_index(r::RootSpaceElem)
  i = findfirst(==(r), roots(r.root_system))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_positive_root_with_index(r::RootSpaceElem)
  i = findfirst(==(r), positive_roots(r.root_system))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_negative_root_with_index(r::RootSpaceElem)
  i = findfirst(==(r), negative_roots(r.root_system))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function reflect!(r::RootSpaceElem, s::Int)
  addmul!(
    r.vec,
    root_system(r).positive_roots[s],
    dot(view(cartan_matrix(r.root_system), s, :), r.vec),
  )
  return r
end

function root_system(r::RootSpaceElem)
  return r.root_system
end

###############################################################################
# WeightLatticeElem

struct WeightLatticeElem
  root_system::RootSystem
  vec::ZZMatrix # the coordinate (column) vector with respect to the fundamental weights
end

@doc raw"""
    WeightLatticeElem(R::RootSystem, v::Vector{IntegerUnion}) -> WeightLatticeElem

Returns the weight defined by the coefficients `v` of the fundamental weights with respect to the root system `R`.
"""
function WeightLatticeElem(R::RootSystem, v::Vector{<:IntegerUnion})
  return WeightLatticeElem(R, matrix(ZZ, rank(R), 1, v))
end

function Base.:(*)(n::IntegerUnion, w::WeightLatticeElem)
  return WeightLatticeElem(w.root_system, n * w.vec)
end

function Base.:(+)(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req w.root_system === w2.root_system "$w and $w2 must belong to the same weight lattice"

  return RootSpaceElem(w.root_system, w.vec + w2.vec)
end

function Base.:(-)(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req w.root_system === w2.root_system "$w and $w2 must belong to the same weight lattice"

  return WeightLatticeElem(w.root_system, w.vec - w2.vec)
end

function Base.:(-)(w::WeightLatticeElem)
  return WeightLatticeElem(w.root_system, -w.vec)
end

function Base.:(==)(w::WeightLatticeElem, w2::WeightLatticeElem)
  return w.root_system === w2.root_system && w.vec == w2.vec
end

function Base.deepcopy_internal(w::WeightLatticeElem, dict::IdDict)
  if haskey(dict, w)
    return dict[w]
  end

  w2 = WeightLatticeElem(w.root_system, deepcopy_internal(w.vec, dict))
  dict[w] = w2
  return w2
end

@doc raw"""
  getindex(w::WeightLatticeElem, i::Int) -> ZZRingElem

Returns the coefficient of the `i`-th fundamental weight in `w`.
"""
function Base.getindex(w::WeightLatticeElem, i::Int)
  return coeff(w, i)
end

function Base.hash(w::WeightLatticeElem, h::UInt)
  b = 0x7b2fefadacf46f4e % UInt
  h = hash(w.root_system, h)
  h = hash(w.vec, h)
  return xor(b, h)
end

@doc raw"""
    iszero(w::WeightLatticeElem) -> Bool

Returns whether `w` is zero.
"""
function Base.iszero(w::WeightLatticeElem)
  return iszero(w.vec)
end

function coefficients(w::WeightLatticeElem)
  return w.vec
end

function coeff(w::WeightLatticeElem, i::Int)
  return w.vec[i]
end

@doc raw"""
    conjugate_dominant_weight(w::WeightLatticeElem) -> WeightLatticeElem

Returns the unique dominant weight conjugate to `w`.
"""
function conjugate_dominant_weight(w::WeightLatticeElem)
  # conj will be the dominant weight conjugate to w
  conj = deepcopy(w)

  # conj will be dominant once all fundamental weights have a positive coefficient,
  # so search for negative coefficients and make them positive by applying the corresponding reflection.
  s = 1
  while s <= rank(w.root_system)
    if conj.vec[s] < 0
      reflect!(conj, s)
      s = 1
    else
      s += 1
    end
  end

  return conj
end

function expressify(w::WeightLatticeElem, s=:w; context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(w.vec)
    push!(sum.args, Expr(:call, :*, expressify(w.vec[i]; context=context), "$s$i"))
  end
  return sum
end
@enable_all_show_via_expressify WeightLatticeElem

@doc raw"""
    reflect(w::WeightLatticeElem, s::Int) -> WeightLatticeElem
    
Returns the `w` reflected at the `s`-th simple root.
"""
function reflect(w::WeightLatticeElem, s::Int)
  return reflect!(deepcopy(w), s)
end

@doc raw"""
    reflect!(w::WeightLatticeElem, s::Int) -> WeightLatticeElem
    
Reflects the `w` at the `s`-th simple root in place and returns `w`.
"""
function reflect!(w::WeightLatticeElem, s::Int)
  addmul!(w.vec, view(w.root_system.cartan_matrix, :, s), -w.vec[s])
  return w
end

function root_system(w::WeightLatticeElem)
  return w.root_system
end

###############################################################################
# internal helpers

# cartan matrix in the format <a^v, b>
function positive_roots_and_reflections(cartan_matrix::ZZMatrix)
  rank, _ = size(cartan_matrix)

  roots = [[l == s ? 1 : 0 for l in 1:rank] for s in 1:rank]
  coroots = [[l == s ? 1 : 0 for l in 1:rank] for s in 1:rank]
  rootidx = Dict(roots[s] => s for s in 1:rank)
  refl = Dict((s, s) => 0 for s in 1:rank)

  i = 1
  while i <= length(roots)
    for s in 1:rank
      if haskey(refl, (s, i))
        continue
      end

      pairing = sum(roots[i][l] * cartan_matrix[s, l] for l in 1:rank)
      copairing = sum(coroots[i][l] * cartan_matrix[l, s] for l in 1:rank)
      if pairing * copairing >= 4
        refl[s, i] = 0
        continue
      end

      r = copy(roots[i])
      r[s] -= pairing

      # If we've not seen this root before, then add it to our list.
      if !haskey(rootidx, r)
        push!(roots, r)
        rootidx[r] = length(roots)

        new_coroot = copy(coroots[i])
        new_coroot[s] -= copairing
        push!(coroots, new_coroot)
      end

      # record the reflection data
      si = rootidx[r]
      refl[s, i] = si
      refl[s, si] = i
    end
    i += 1
  end

  table = zero_matrix(ZZ, rank, length(roots))
  for i in 1:length(roots), s in 1:rank
    table[s, i] = refl[s, i]
  end

  roots, table
end
