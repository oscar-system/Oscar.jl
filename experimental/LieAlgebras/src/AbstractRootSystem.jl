###############################################################################
#
#   Root Systems and Weights
#
###############################################################################

mutable struct AbstractRootSystem{T<:CoxeterGroup} <: RootSystem
  # familiy (A, ..., G) of the root system
  # fam::Vector{Tuple{Symbol, Int}}

  # (generalized) Cartan matrix
  cartan_matrix::ZZMatrix

  # fundamental weights as linear combination of simple roots
  #fw::QQMatrix

  positive_roots::Vector

  weyl_group::T

  function AbstractRootSystem(mat::ZZMatrix)
    R = new{WeylGroup}()

    roots, refl = positive_roots_and_reflections(mat)
    finite = count(refl .== 0) == nrows(mat)

    R.cartan_matrix = mat
    R.positive_roots = roots
    R.weyl_group = WeylGroup(finite, refl, R)

    return R
  end
end

@doc raw"""
    root_system(cartan_matrix::ZZMatrix) -> AbstractRootSystem

Constructs the root system defined by the Cartan matrix.
"""
function root_system(cartan_matrix::ZZMatrix)
  #fw = solve(change_base_ring(QQ, gcm), identity_matrix(QQ, nrows(gcm)))
  return AbstractRootSystem(ZZMatrix(cartan_matrix))
end

@doc raw"""
    root_system(fam::Symbol, rk::Int) -> AbstractRootSystem

Constructs the root system of the given type. See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.
"""
function root_system(fam::Symbol, rk::Int)
  cartan = cartan_matrix(fam, rk)
  #fw = solve(change_base_ring(QQ, cartan), identity_matrix(QQ, nrows(cartan)))
  return AbstractRootSystem(cartan)
end

function root_system(types::Tuple{Symbol,Int}...)
  return root_system(cartan_matrix(types...))
end

function Base.show(io::IO, R::RootSystem)
  println(io, "root system defined by Cartan matrix")
  show(io, R.cartan_matrix)
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
#  return [RootLatticeElem(R, matrix(QQ, 1, rk, i .== it)) for i in 1:rk]
#end

@doc raw"""
    fundamental_weights(R::RootSystem) -> Vector{WeightLatticeElem}

Returns the fundamental weights corresponding to the simple roots of `R`.
"""
function fundamental_weights(R::RootSystem)
  it = 1:rank(R)
  return [WeightLatticeElem(R, matrix(ZZ, rank(R), 1, i .== it)) for i in 1:rank(R)]
end

#function positive_roots(R::RootSystem) end

@doc raw"""
    rank(R::RootSystem) -> Int

Returns the rank of `R`.
"""
function rank(R::RootSystem)
  return nrows(R.cartan_matrix)
end

@doc raw"""
    weyl_group(R::RootSystem) -> WeylGroup

Returns the Weyl group of `R`.
"""
function weyl_group(R::RootSystem)
  return R.weyl_group
end

@doc raw"""
    weyl_vector(R::RootSystem) -> WeightLatticeElem

Returns the Weyl vector $\rho$ of `R`, which is the sum of all fundamental weights, or half the sum of all positive roots.
"""
function weyl_vector(R::RootSystem)
  rk = rank(R)
  return WeightLatticeElem(R, matrix(ZZ, rk, 1, fill(1, rk)))
end

###############################################################################
# WIP

struct RootSpaceElem
  root_system::AbstractRootSystem

  # vec is the coordinate (row) vector with respect to the simple roots
  vec::QQMatrix
end

function Base.:(*)(q::RationalUnion, r::RootSpaceElem)
  return RootSpaceElem(r.root_system, q * r.vec)
end

@doc raw"""
    getindex(r::RootSpaceElem, i::Int) -> QQRingElem

Returns the coefficient of the `i`-th simple root in `r`.
"""
function Base.getindex(r::RootSpaceElem, i::Int)
  return r.vec[i]
end

function reflect!(r::RootSpaceElem, s::Int)
  addmul!(r.vec, r.root_system.positive_roots[s], dot(view(w.root_system.cartan_matrix, s, :), r.vec))
  return r
end

###############################################################################
# Weights

struct WeightLatticeElem
  root_system::RootSystem

  # vec is the coordinate (column) vector with respect to the fundamental weights
  vec::ZZMatrix
end

@doc raw"""
    weight(R::RootSystem, v::Vector{IntegerUnion}) -> WeightLatticeElem

Returns the weight defined by the coefficients `v` of the fundamental weights with respect to the root system `R`.
"""
function weight(R::RootSystem, v::Vector{IntegerUnion})
  return WeightLatticeElem(R, matrix(ZZ, rank(R), 1, v))
end

function Base.:(*)(n::IntegerUnion, w::WeightLatticeElem)
  return WeightLatticeElem(w.root_system, n * w.vec)
end

function Base.:(+)(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req w.root_system === w2.root_system "$w and $w2 must belong to the same weight lattice"

  return RootLatticeElem(w.root_system, w.vec + w2.vec)
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

function Base.copy(w::WeightLatticeElem)
  return WeightLatticeElem(w.root_system, deepcopy(w.vec))
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
  return w.vec[i]
end

@doc raw"""
    iszero(w::WeightLatticeElem) -> Bool

Returns whether `w` is zero.
"""
function Base.iszero(w::WeightLatticeElem)
  return iszero(w.vec)
end

@doc raw"""
    conjugate_dominant_weight(w::WeightLatticeElem) -> WeightLatticeElem

Returns the unique dominant weight conjugate to `w`.
"""
function conjugate_dominant_weight(w::WeightLatticeElem)
  # conj will be the dominant weight conjugate to w
  conj = copy(w)

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
  return reflect!(WeightLatticeElem(w.root_system, deepcopy(w.vec)), s)
end

@doc raw"""
    reflect!(w::WeightLatticeElem, s::Int) -> WeightLatticeElem
    
Reflects the `w` at the `s`-th simple root in place and returns `w`.
"""
function reflect!(w::WeightLatticeElem, s::Int)
  addmul!(w.vec, view(w.root_system.cartan_matrix, :, s), -w.vec[s])
  return w
end

###############################################################################
# internal helpers

# cartan matrix in the format <a^v, b>
function positive_roots_and_reflections(cartan_matrix::ZZMatrix)
  rank, _ = size(cartan_matrix)

  roots = [[l == s ? 1 : 0 for l in 1:rank] for s in 1:rank]
  rootidx = Dict(roots[s] => s for s in 1:rank)
  refl = Dict((s, s) => 0 for s in 1:rank)

  i = 1
  while i <= length(roots)
    for s in 1:rank
      if haskey(refl, (s, i))
        continue
      end

      pairing = sum(roots[i][l] * cartan_matrix[s, l] for l in 1:rank)
      copairing = sum(roots[i][l] * cartan_matrix[l, s] for l in 1:rank)
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
