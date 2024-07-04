###############################################################################
#
#   Root Systems and Weights
#
###############################################################################

mutable struct RootSystem
  cartan_matrix::ZZMatrix # (generalized) Cartan matrix
  #fw::QQMatrix # fundamental weights as linear combination of simple roots
  positive_roots::Vector #::Vector{RootSpaceElem} (cyclic reference)
  positive_coroots::Vector #::Vector{DualRootSpaceElem} (cyclic reference)
  weyl_group::Any     #::WeylGroup (cyclic reference)

  # optional:
  type::Vector{Tuple{Symbol,Int}}
  type_ordering::Vector{Int}

  function RootSystem(mat::ZZMatrix)
    pos_roots, pos_coroots, refl = positive_roots_and_reflections(mat)
    finite = count(refl .== 0) == nrows(mat)

    R = new(mat)
    R.positive_roots = map(r -> RootSpaceElem(R, r), pos_roots)
    R.positive_coroots = map(r -> DualRootSpaceElem(R, r), pos_coroots)
    R.weyl_group = WeylGroup(finite, refl, R)

    return R
  end
end

@doc raw"""
    root_system(cartan_matrix::ZZMatrix; check::Bool=true, detect_type::Bool=true) -> RootSystem
    root_system(cartan_matrix::Matrix{Int}; check::Bool=true, detect_type::Bool=true) -> RootSystem

Constructs the root system defined by the Cartan matrix.
If `check` is `true`, checks that `cartan_matrix` is a generalized Cartan matrix.
Passing `detect_type=false` will skip the detection of the root system type.
"""
function root_system(cartan_matrix::ZZMatrix; check::Bool=true, detect_type::Bool=true)
  @req !check || is_cartan_matrix(cartan_matrix) "Requires a generalized Cartan matrix"
  R = RootSystem(cartan_matrix)
  detect_type && is_finite(weyl_group(R)) && set_root_system_type(R, cartan_type_with_ordering(cartan_matrix)...)
  return R
end

function root_system(cartan_matrix::Matrix{<:Integer}; kwargs...)
  return root_system(matrix(ZZ, cartan_matrix); kwargs...)
end

@doc raw"""
    root_system(fam::Symbol, rk::Int) -> RootSystem

Constructs the root system of the given type. See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.
"""
function root_system(fam::Symbol, rk::Int)
  cartan = cartan_matrix(fam, rk)
  R = root_system(cartan; check=false, detect_type=false)
  set_root_system_type(R, [(fam, rk)])
  return R
end

function root_system(type::Vector{Tuple{Symbol,Int}})
  cartan = cartan_matrix(type)
  R = root_system(cartan; check=false, detect_type=false)
  set_root_system_type(R, type)
  return R
end

function root_system(type::Tuple{Symbol,Int}...)
  return root_system(collect(type))
end

function Base.show(io::IO, ::MIME"text/plain", R::RootSystem)
  io = pretty(io)
  println(io, "Root system defined by Cartan matrix")
  print(io, Indent())
  show(io, MIME"text/plain"(), cartan_matrix(R))
  print(io, Dedent())
end

function Base.show(io::IO, R::RootSystem)
  if get(io, :supercompact, false)
    print(io, "Root system")
  else
    print(io, "Root system defined by Cartan matrix $(cartan_matrix(R))")
  end
end

@attr ZZMatrix function bilinear_form(R::RootSystem)
  return cartan_bilinear_form(cartan_matrix(R); check=false)
end

@doc raw"""
    cartan_matrix(R::RootSystem) -> ZZMatrix

Returns the Cartan matrix defining `R`.
"""
function cartan_matrix(R::RootSystem)
  return R.cartan_matrix
end

@doc raw"""
    coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the `i`-th coroot of `R`, i.e. the `i`-th root of the dual root system of `R`.
This is a more efficient version for `coroots(R)[i]`.

Also see: `coroots`.
"""
function coroot(R::RootSystem, i::Int)
  if i <= n_positive_roots(R)
    return positive_coroot(R, i)
  else
    return negative_coroot(R, i - n_positive_roots(R))
  end
end

@doc raw"""
    coroots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the coroots of `R`, starting with the coroots of positive roots and then the negative roots,
in the order of `positive_coroots` and `negative_coroots`.

Also see: `coroot`.
"""
function coroots(R::RootSystem)
  return [[r for r in positive_coroots(R)]; [-r for r in positive_coroots(R)]]
end

function fundamental_weight(R::RootSystem, i::Int)
  @req 1 <= i <= rank(R) "invalid index"
  return WeightLatticeElem(R, matrix(ZZ, rank(R), 1, i .== 1:rank(R)))
end

@doc raw"""
    fundamental_weights(R::RootSystem) -> Vector{WeightLatticeElem}

Returns the fundamental weights corresponding to the `simple_roots` of `R`.
"""
function fundamental_weights(R::RootSystem)
  return [fundamental_weight(R, i) for i in 1:rank(R)]
end

function is_simple(R::RootSystem)
  if is_finite(weyl_group(R))
    return length(root_system_type(R)) == 1
  end
  error("Not implemented") # TODO: implement is_simple
end

@doc raw"""
    negative_root(R::RootSystem, i::Int) -> RootSpaceElem

Returns the `i`-th negative root of `R`.
This is a more efficient version for `negative_roots(R)[i]`.

Also see: `negative_roots`.
"""
function negative_root(R::RootSystem, i::Int)
  return -R.positive_roots[i]::RootSpaceElem
end

@doc raw"""
    negative_roots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the negative roots of `R`. The $i$-th element of the returned vector is the negative root corresponding to the $i$-th positive root.

Also see: `negative_root`.
"""
function negative_roots(R::RootSystem)
  return [-r for r in positive_roots(R)]
end

@doc raw"""
    negative_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroot corresponding to the `i`-th negative root of `R`
This is a more efficient version for `negative_coroots(R)[i]`.

Also see: `negative_coroots`.
"""
function negative_coroot(R::RootSystem, i::Int)
  return -R.positive_coroots[i]::DualRootSpaceElem
end

@doc raw"""
    negative_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroots corresponding to the negative roots of `R`

Also see: `negative_coroots`.
"""
function negative_coroots(R::RootSystem)
  return [-r for r in positive_coroots(R)]
end

@doc raw"""
    number_of_positive_roots(R::RootSystem) -> Int

Returns the number of positive roots of `R`. This is the same as the number of negative roots.

Also see: `positive_roots`, `negative_roots`.
"""
function number_of_positive_roots(R::RootSystem)
  return length(positive_roots(R))
end

@doc raw"""
    number_of_roots(R::RootSystem) -> Int

Returns the number of roots of `R`.

Also see: `roots`.
"""
function number_of_roots(R::RootSystem)
  return 2 * number_of_positive_roots(R)
end

@doc raw"""
    number_of_simple_roots(R::RootSystem) -> Int

Returns the number of simple roots of `R`.

Also see: `simple_roots`.
"""
function number_of_simple_roots(R::RootSystem)
  return rank(R)
end

@doc raw"""
    positive_root(R::RootSystem, i::Int) -> RootSpaceElem

Returns the `i`-th positive root of `R`.
This is a more efficient version for `positive_roots(R)[i]`.

Also see: `positive_roots`.
"""
function positive_root(R::RootSystem, i::Int)
  return R.positive_roots[i]::RootSpaceElem
end

@doc raw"""
    positive_roots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the positive roots of `R`, starting with the simple roots in the order of `simple_roots`,
and then increasing in height.

Also see: `positive_root`, `number_of_positive_roots`.
"""
function positive_roots(R::RootSystem)
  return R.positive_roots::Vector{RootSpaceElem}
end

@doc raw"""
    positive_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroot corresponding to the `i`-th positive root of `R`
This is a more efficient version for `positive_coroots(R)[i]`.

Also see: `positive_coroots`.
"""
function positive_coroot(R::RootSystem, i::Int)
  return R.positive_coroots[i]::DualRootSpaceElem
end

@doc raw"""
    positive_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroots corresponding to the positive roots of `R`

Also see: `positive_coroots`.
"""
function positive_coroots(R::RootSystem)
  return R.positive_coroots::Vector{DualRootSpaceElem}
end

@doc raw"""
    rank(R::RootSystem) -> Int

Returns the rank of `R`, i.e. the number of simple roots.
"""
function rank(R::RootSystem)
  return nrows(cartan_matrix(R))
end

function root_system_type(R::RootSystem)
  has_root_system_type(R) || error("Root system type not known and cannot be determined")
  return R.type
end

function root_system_type_with_ordering(R::RootSystem)
  return R.type, R.type_ordering
end

function has_root_system_type(R::RootSystem)
  return isdefined(R, :type) && isdefined(R, :type_ordering)
end

function set_root_system_type(R::RootSystem, type::Vector{Tuple{Symbol,Int}})
  return set_root_system_type(R, type, 1:sum(t[2] for t in type; init=0))
end

function set_root_system_type(
  R::RootSystem, type::Vector{Tuple{Symbol,Int}}, ordering::AbstractVector{Int}
)
  R.type = type
  R.type_ordering = collect(ordering)
  return nothing
end

function root_system_type_string(R::RootSystem)
  return join([string(t[1]) * string(t[2]) for t in root_system_type(R)], " x ")
end

@doc raw"""
    root(R::RootSystem, i::Int) -> RootSpaceElem

Returns the `i`-th root of `R`.
This is a more efficient version for `roots(R)[i]`.

Also see: `roots`.
"""
function root(R::RootSystem, i::Int)
  if i <= n_positive_roots(R)
    return positive_root(R, i)
  else
    return negative_root(R, i - n_positive_roots(R))
  end
end

@doc raw"""
    roots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the roots of `R`, starting with the positive roots and then the negative roots,
in the order of `positive_roots` and `negative_roots`.

Also see: `root`.
"""
function roots(R::RootSystem)
  return [[r for r in positive_roots(R)]; [-r for r in positive_roots(R)]]
end

@doc raw"""
    simple_root(R::RootSystem, i::Int) -> RootSpaceElem

Returns the `i`-th simple root of `R`.
This is a more efficient version for `simple_roots(R)[i]`.

Also see: `simple_roots`.
"""
function simple_root(R::RootSystem, i::Int)
  @req 1 <= i <= rank(R) "Invalid index"
  return positive_root(R, i)
end

@doc raw"""
    simple_roots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the simple roots of `R`.

Also see: `simple_root`.
"""
function simple_roots(R::RootSystem)
  return positive_roots(R)[1:rank(R)]
end

@doc raw"""
    simple_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroot corresponding to the `i`-th simple root of `R`
This is a more efficient version for `simple_coroots(R)[i]`.

Also see: `simple_coroots`.
"""
function simple_coroot(R::RootSystem, i::Int)
  @req 1 <= i <= rank(R) "Invalid index"
  return positive_coroot(R, i)
end

@doc raw"""
    simple_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroots corresponding to the simple roots of `R`

Also see: `simple_coroots`.
"""
function simple_coroots(R::RootSystem)
  return positive_coroots(R)[1:rank(R)]
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

mutable struct RootSpaceElem
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
  @req root_system(r) === root_system(r2) "$r and $r2 must belong to the same root space"

  return RootSpaceElem(root_system(r), r.vec + r2.vec)
end

function Base.:(-)(r::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(r) === root_system(r2) "$r and $r2 must belong to the same root space"

  return RootSpaceElem(root_system(r), r.vec - r2.vec)
end

function Base.:(-)(r::RootSpaceElem)
  return RootSpaceElem(root_system(r), -r.vec)
end

function Base.:(==)(r::RootSpaceElem, r2::RootSpaceElem)
  return r.root_system === r2.root_system && r.vec == r2.vec
end

function Base.deepcopy_internal(r::RootSpaceElem, dict::IdDict)
  if haskey(dict, r)
    return dict[r]
  end

  w2 = RootSpaceElem(root_system(r), deepcopy_internal(r.vec, dict))
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

function dot(r1::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(r1) === root_system(r2) "$r1 and $r2 must belong to the same root space"

  return dot(coefficients(r1) * bilinear_form(root_system(r1)), coefficients(r2))
end

@doc raw"""
    height(r::RootSpaceElem) -> QQFieldElem

For a root `r`, returns the height of `r`, i.e. the sum of the coefficients of the simple roots.
If `r` is not a root, the return value is arbitrary.
"""
function height(r::RootSpaceElem)
  return sum(coefficients(r))
end

function is_root_with_index(r::RootSpaceElem)
  i = findfirst(==(r), roots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_positive_root_with_index(r::RootSpaceElem)
  i = findfirst(==(r), positive_roots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_negative_root_with_index(r::RootSpaceElem)
  i = findfirst(==(r), negative_roots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_simple_root_with_index(r::RootSpaceElem)
  i = findfirst(==(r), simple_roots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function Base.iszero(r::RootSpaceElem)
  return iszero(r.vec)
end

function reflect!(r::RootSpaceElem, s::Int)
  r.vec -=
    dot(view(cartan_matrix(root_system(r)), s, :), r.vec) *
    positive_root(root_system(r), s).vec
  return r
end

function root_system(r::RootSpaceElem)
  return r.root_system
end

###############################################################################
# DualRootSpaceElem

mutable struct DualRootSpaceElem
  root_system::RootSystem
  vec::QQMatrix # the coordinate (row) vector with respect to the simple coroots
end

function DualRootSpaceElem(root_system::RootSystem, vec::Vector{<:RationalUnion})
  return DualRootSpaceElem(root_system, matrix(QQ, 1, length(vec), vec))
end

function Base.:(*)(q::RationalUnion, r::DualRootSpaceElem)
  return DualRootSpaceElem(root_system(r), q * r.vec)
end

function Base.:(+)(r::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(r) === root_system(r2) "$r and $r2 must belong to the same root space"

  return DualRootSpaceElem(root_system(r), r.vec + r2.vec)
end

function Base.:(-)(r::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(r) === root_system(r2) "$r and $r2 must belong to the same root space"

  return DualRootSpaceElem(root_system(r), r.vec - r2.vec)
end

function Base.:(-)(r::DualRootSpaceElem)
  return DualRootSpaceElem(root_system(r), -r.vec)
end

function Base.:(==)(r::DualRootSpaceElem, r2::DualRootSpaceElem)
  return r.root_system === r2.root_system && r.vec == r2.vec
end

function Base.deepcopy_internal(r::DualRootSpaceElem, dict::IdDict)
  if haskey(dict, r)
    return dict[r]
  end

  w2 = DualRootSpaceElem(root_system(r), deepcopy_internal(r.vec, dict))
  dict[r] = w2
  return w2
end

@doc raw"""
    getindex(r::DualRootSpaceElem, i::Int) -> QQRingElem

Returns the coefficient of the `i`-th simple root in `r`.
"""
function Base.getindex(r::DualRootSpaceElem, i::Int)
  return coeff(r, i)
end

function Base.hash(r::DualRootSpaceElem, h::UInt)
  b = 0x721bec0418bdbe0f % UInt
  h = hash(r.root_system, h)
  h = hash(r.vec, h)
  return xor(b, h)
end

function coefficients(r::DualRootSpaceElem)
  return r.vec
end

function coeff(r::DualRootSpaceElem, i::Int)
  return r.vec[i]
end

@doc raw"""
    height(r::DualRootSpaceElem) -> QQFieldElem

For a coroot `r`, returns the height of `r`, i.e. the sum of the coefficients of the simple coroots.
If `r` is not a coroot, the return value is arbitrary.
"""
function height(r::DualRootSpaceElem)
  return sum(coefficients(r))
end

function is_coroot_with_index(r::DualRootSpaceElem)
  i = findfirst(==(r), coroots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_positive_coroot_with_index(r::DualRootSpaceElem)
  i = findfirst(==(r), positive_coroots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_negative_coroot_with_index(r::DualRootSpaceElem)
  i = findfirst(==(r), negative_coroots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_simple_coroot_with_index(r::DualRootSpaceElem)
  i = findfirst(==(r), simple_coroots(root_system(r)))
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function Base.iszero(r::DualRootSpaceElem)
  return iszero(r.vec)
end

function root_system(r::DualRootSpaceElem)
  return r.root_system
end

###############################################################################
# WeightLatticeElem

mutable struct WeightLatticeElem
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
  return WeightLatticeElem(root_system(w), n * w.vec)
end

function Base.:(+)(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req root_system(w) === root_system(w2) "$w and $w2 must belong to the same weight lattice"

  return RootSpaceElem(root_system(w), w.vec + w2.vec)
end

function Base.:(-)(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req root_system(w) === root_system(w2) "$w and $w2 must belong to the same weight lattice"

  return WeightLatticeElem(root_system(w), w.vec - w2.vec)
end

function Base.:(-)(w::WeightLatticeElem)
  return WeightLatticeElem(root_system(w), -w.vec)
end

function Base.:(==)(w::WeightLatticeElem, w2::WeightLatticeElem)
  return w.root_system === w2.root_system && w.vec == w2.vec
end

function Base.deepcopy_internal(w::WeightLatticeElem, dict::IdDict)
  if haskey(dict, w)
    return dict[w]
  end

  w2 = WeightLatticeElem(root_system(w), deepcopy_internal(w.vec, dict))
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
  while s <= rank(root_system(w))
    if conj.vec[s] < 0
      reflect!(conj, s)
      s = 1
    else
      s += 1
    end
  end

  return conj
end

@doc raw"""
    conjugate_dominant_weight_with_elem(w::WeightLatticeElem) -> Tuple{WeightLatticeElem, WeylGroupElem}

Returns the unique dominant weight `dom` conjugate to `w` and a Weyl group element `x`
such that `x*w == dom`.
"""
function conjugate_dominant_weight_with_elem(w::WeightLatticeElem)
  R = root_system(w)
  wt = deepcopy(w)

  # determine the Weyl group element taking w to the fundamental chamber
  word = sizehint!(UInt8[], count(<(0), coefficients(wt))^2)
  s = 1
  while s <= rank(R)
    if wt[s] < 0
      push!(word, UInt8(s))
      reflect!(wt, s)
      s = 1
    else
      s += 1
    end
  end

  # reversing word means it is in short revlex normal form
  # and it is the element taking w to wt
  return wt, weyl_group_elem(R, reverse!(word); normalize=false)
end

function expressify(w::WeightLatticeElem, s=:w; context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(w.vec)
    push!(sum.args, Expr(:call, :*, expressify(w.vec[i]; context), "$s$i"))
  end
  return sum
end
@enable_all_show_via_expressify WeightLatticeElem

function is_dominant(w::WeightLatticeElem)
  return all(>=(0), coefficients(w))
end

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
  addmul!(w.vec, view(cartan_matrix(root_system(w)), :, s:s), -w.vec[s])
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

  # sort roots by height
  perm = sortperm(roots; by=sum)
  invp = invperm(perm)

  table = zeros(UInt, rank, length(roots))
  for i in 1:length(roots), s in 1:rank
    table[s, i] = iszero(refl[s, perm[i]]) ? 0 : invp[refl[s, perm[i]]]
  end

  roots[perm], coroots[perm], table
end
