###############################################################################
#
#   Root Systems
#
###############################################################################

@doc raw"""
    root_system(cartan_matrix::ZZMatrix; check::Bool=true, detect_type::Bool=true) -> RootSystem
    root_system(cartan_matrix::Matrix{Int}; check::Bool=true, detect_type::Bool=true) -> RootSystem

Construct the root system defined by the Cartan matrix.
If `check` is `true`, checks that `cartan_matrix` is a generalized Cartan matrix.
Passing `detect_type=false` will skip the detection of the root system type.
"""
function root_system(cartan_matrix::ZZMatrix; check::Bool=true, detect_type::Bool=true)
  return RootSystem(cartan_matrix; check, detect_type)
end

function root_system(cartan_matrix::Matrix{<:Integer}; kwargs...)
  return root_system(matrix(ZZ, cartan_matrix); kwargs...)
end

@doc raw"""
    root_system(fam::Symbol, rk::Int) -> RootSystem

Construct the root system of the given type. See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.

# Examples
```jldoctest
julia> root_system(:A, 2)
Root system defined by Cartan matrix
  [ 2   -1]
  [-1    2]
```
"""
function root_system(fam::Symbol, rk::Int)
  cartan = cartan_matrix(fam, rk)
  R = root_system(cartan; check=false, detect_type=false)
  set_root_system_type!(R, [(fam, rk)])
  return R
end

function root_system(type::Vector{Tuple{Symbol,Int}})
  cartan = cartan_matrix(type)
  R = root_system(cartan; check=false, detect_type=false)
  set_root_system_type!(R, type)
  return R
end

function root_system(type::Tuple{Symbol,Int}...)
  return root_system(collect(type))
end

function Base.show(io::IO, mime::MIME"text/plain", R::RootSystem)
  @show_name(io, R)
  @show_special(io, mime, R)
  io = pretty(io)
  println(io, "Root system defined by Cartan matrix")
  print(io, Indent())
  show(io, mime, cartan_matrix(R))
  print(io, Dedent())
end

function Base.show(io::IO, R::RootSystem)
  @show_name(io, R)
  @show_special(io, R)
  if is_terse(io)
    print(io, "Root system")
  else
    print(io, "Root system defined by Cartan matrix $(cartan_matrix(R))")
  end
end

@attr ZZMatrix function bilinear_form(R::RootSystem)
  return cartan_bilinear_form(cartan_matrix(R); check=false)
end

@attr QQMatrix function _bilinear_form_QQ(R::RootSystem)
  return QQMatrix(bilinear_form(R))
end

@doc raw"""
    cartan_matrix(R::RootSystem) -> ZZMatrix

Return the Cartan matrix defining `R`.
"""
function cartan_matrix(R::RootSystem)
  return R.cartan_matrix
end

@attr QQMatrix function cartan_matrix_inv(R::RootSystem)
  return inv(matrix(QQ, cartan_matrix(R)))
end

@attr QQMatrix function cartan_matrix_inv_tr(R::RootSystem)
  return transpose(cartan_matrix_inv(R))
end

@attr ZZMatrix function cartan_matrix_tr(R::RootSystem)
  return transpose(cartan_matrix(R))
end

@attr Vector{ZZRingElem} function cartan_symmetrizer(R::RootSystem)
  return cartan_symmetrizer(cartan_matrix(R); check=false)
end

@attr ZZMatrix function _cartan_symmetrizer_mat(R::RootSystem)
  return diagonal_matrix(ZZ, cartan_symmetrizer(R))
end

@doc raw"""
    coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the `i`-th coroot of `R`, i.e. the `i`-th root of the dual root system of `R`.
This is a more efficient version for `coroots(R)[i]`.

Also see: `coroots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
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

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
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

Return the fundamental weights corresponding to the `simple_roots` of `R`.
"""
function fundamental_weights(R::RootSystem)
  return [fundamental_weight(R, i) for i in 1:rank(R)]
end

function Base.hash(R::RootSystem, h::UInt)
  # even though we don't have a == method for RootSystem, we add a hash method
  # to make hashing of RootSpaceElem and WeightLatticeElem more deterministic
  b = 0xeb5362118dea2a0e % UInt
  h = hash(cartan_matrix(R), h)
  return xor(b, h)
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

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function negative_root(R::RootSystem, i::Int)
  return -(R.positive_roots::Vector{RootSpaceElem})[i]
end

@doc raw"""
    negative_roots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the negative roots of `R`. The $i$-th element of the returned vector is the negative root corresponding to the $i$-th positive root.

Also see: `negative_root`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function negative_roots(R::RootSystem)
  return [-r for r in positive_roots(R)]
end

@doc raw"""
    negative_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroot corresponding to the `i`-th negative root of `R`
This is a more efficient version for `negative_coroots(R)[i]`.

Also see: `negative_coroots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function negative_coroot(R::RootSystem, i::Int)
  return -(R.positive_coroots::Vector{DualRootSpaceElem})[i]
end

@doc raw"""
    negative_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroots corresponding to the negative roots of `R`

Also see: `negative_coroots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
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

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function positive_root(R::RootSystem, i::Int)
  return (R.positive_roots::Vector{RootSpaceElem})[i]
end

@doc raw"""
    positive_roots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the positive roots of `R`, starting with the simple roots in the order of `simple_roots`,
and then increasing in height.

Also see: `positive_root`, `number_of_positive_roots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function positive_roots(R::RootSystem)
  return R.positive_roots::Vector{RootSpaceElem}
end

@doc raw"""
    positive_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroot corresponding to the `i`-th positive root of `R`
This is a more efficient version for `positive_coroots(R)[i]`.

Also see: `positive_coroots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function positive_coroot(R::RootSystem, i::Int)
  return (R.positive_coroots::Vector{DualRootSpaceElem})[i]
end

@doc raw"""
    positive_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroots corresponding to the positive roots of `R`

Also see: `positive_coroots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function positive_coroots(R::RootSystem)
  return R.positive_coroots::Vector{DualRootSpaceElem}
end

@doc raw"""
    rank(R::RootSystem) -> Int

Return the rank of `R`, i.e. the number of simple roots.
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

function set_root_system_type!(R::RootSystem, type::Vector{Tuple{Symbol,Int}})
  return set_root_system_type!(R, type, 1:sum(t[2] for t in type; init=0))
end

function set_root_system_type!(
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

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
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

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function roots(R::RootSystem)
  return [[r for r in positive_roots(R)]; [-r for r in positive_roots(R)]]
end

@doc raw"""
    simple_root(R::RootSystem, i::Int) -> RootSpaceElem

Returns the `i`-th simple root of `R`.
This is a more efficient version for `simple_roots(R)[i]`.

Also see: `simple_roots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function simple_root(R::RootSystem, i::Int)
  @req 1 <= i <= rank(R) "Invalid index"
  return positive_root(R, i)
end

@doc raw"""
    simple_roots(R::RootSystem) -> Vector{RootSpaceElem}

Returns the simple roots of `R`.

Also see: `simple_root`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function simple_roots(R::RootSystem)
  return positive_roots(R)[1:rank(R)]
end

@doc raw"""
    simple_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroot corresponding to the `i`-th simple root of `R`
This is a more efficient version for `simple_coroots(R)[i]`.

Also see: `simple_coroots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function simple_coroot(R::RootSystem, i::Int)
  @req 1 <= i <= rank(R) "Invalid index"
  return positive_coroot(R, i)
end

@doc raw"""
    simple_coroot(R::RootSystem, i::Int) -> RootSpaceElem

Returns the coroots corresponding to the simple roots of `R`

Also see: `simple_coroots`.

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function simple_coroots(R::RootSystem)
  return positive_coroots(R)[1:rank(R)]
end

@doc raw"""
    weyl_group(R::RootSystem) -> WeylGroup

Return the Weyl group of `R`.
"""
function weyl_group(R::RootSystem)
  return R.weyl_group::WeylGroup
end

@doc raw"""
    weyl_vector(R::RootSystem) -> WeightLatticeElem

Return the Weyl vector $\rho$ of `R`, which is the sum of all fundamental weights,
or half the sum of all positive roots.
"""
function weyl_vector(R::RootSystem)
  return WeightLatticeElem(R, matrix(ZZ, rank(R), 1, fill(1, rank(R))))
end

###############################################################################
#
#   Root space elements
#
###############################################################################

function RootSpaceElem(root_system::RootSystem, vec::Vector{<:RationalUnion})
  return RootSpaceElem(root_system, matrix(QQ, 1, length(vec), vec))
end

function RootSpaceElem(R::RootSystem, w::WeightLatticeElem)
  @req root_system(w) === R "Root system mismatch"
  coeffs = transpose!(cartan_matrix_inv(R) * coefficients(w))
  return RootSpaceElem(R, matrix(QQ, coeffs))
end

function RootSpaceElem(w::WeightLatticeElem)
  return RootSpaceElem(root_system(w), w)
end

function zero(::Type{RootSpaceElem}, R::RootSystem)
  return RootSpaceElem(R, zero_matrix(QQ, 1, rank(R)))
end

function zero(r::RootSpaceElem)
  return zero(RootSpaceElem, root_system(r))
end

function Base.:*(q::RationalUnion, r::RootSpaceElem)
  return RootSpaceElem(root_system(r), q * r.vec)
end

function Base.:+(r::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(r) === root_system(r2) "parent root system mismatch"

  return RootSpaceElem(root_system(r), r.vec + r2.vec)
end

function Base.:-(r::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(r) === root_system(r2) "parent root system mismatch"

  return RootSpaceElem(root_system(r), r.vec - r2.vec)
end

function Base.:-(r::RootSpaceElem)
  return RootSpaceElem(root_system(r), -r.vec)
end

function zero!(r::RootSpaceElem)
  r.vec = zero!(r.vec)
  return r
end

function add!(rr::RootSpaceElem, r1::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  rr.vec = add!(rr.vec, r1.vec, r2.vec)
  return rr
end

function neg!(rr::RootSpaceElem, r::RootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  rr.vec = neg!(rr.vec, r.vec)
  return rr
end

function sub!(rr::RootSpaceElem, r1::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  rr.vec = sub!(rr.vec, r1.vec, r2.vec)
  return rr
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

Return the coefficient of the `i`-th simple root in `r`.
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
  @req root_system(r1) === root_system(r2) "parent root system mismatch"

  # return dot(coefficients(r1) * _bilinear_form_QQ(root_system(r1)), coefficients(r2)) # currently the below is faster
  return only(
    coefficients(r1) * _bilinear_form_QQ(root_system(r1)) * transpose(coefficients(r2))
  )
end

@doc raw"""
    height(r::RootSpaceElem) -> QQFieldElem

For a root `r`, returns the height of `r`, i.e. the sum of the coefficients of the simple roots.
If `r` is not a root, the return value is arbitrary.
"""
function height(r::RootSpaceElem)
  return sum(coefficients(r))
end

function is_root(r::RootSpaceElem)
  return is_positive_root(r) || is_negative_root(r)
end

function is_root_with_index(r::RootSpaceElem)
  fl, i = is_positive_root_with_index(r)
  if fl
    return true, i
  end
  fl, j = is_negative_root_with_index(r)
  if fl
    return true, j + number_of_positive_roots(root_system(r))
  end
  return false, 0
end

function is_positive_root(r::RootSpaceElem)
  return haskey(root_system(r).positive_roots_map, coefficients(r))
end

function is_positive_root_with_index(r::RootSpaceElem)
  i = get(root_system(r).positive_roots_map, coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_negative_root(r::RootSpaceElem)
  return haskey(root_system(r).positive_roots_map, -coefficients(r))
end

function is_negative_root_with_index(r::RootSpaceElem)
  i = get(root_system(r).positive_roots_map, -coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_simple_root(r::RootSpaceElem)
  return is_simple_root_with_index(r)[1]
end

function is_simple_root_with_index(r::RootSpaceElem)
  i = get(root_system(r).positive_roots_map, coefficients(r), nothing)
  if isnothing(i) || i > number_of_simple_roots(root_system(r))
    return false, 0
  else
    return true, i
  end
end

function Base.iszero(r::RootSpaceElem)
  return iszero(r.vec)
end

function reflect(r::RootSpaceElem, s::Int)
  return reflect!(deepcopy(r), s)
end

function reflect!(r::RootSpaceElem, s::Int)
  r.vec -=
    dot(view(cartan_matrix(root_system(r)), s, :), r.vec) *
    simple_root(root_system(r), s).vec
  return r
end

function root_system(r::RootSpaceElem)
  return r.root_system
end

###############################################################################
#
#   Dual root space elements
#
###############################################################################

function DualRootSpaceElem(root_system::RootSystem, vec::Vector{<:RationalUnion})
  return DualRootSpaceElem(root_system, matrix(QQ, 1, length(vec), vec))
end

function zero(::Type{DualRootSpaceElem}, R::RootSystem)
  return DualRootSpaceElem(R, zero_matrix(QQ, 1, rank(R)))
end

function zero(r::DualRootSpaceElem)
  return zero(DualRootSpaceElem, root_system(r))
end

function Base.:*(q::RationalUnion, r::DualRootSpaceElem)
  return DualRootSpaceElem(root_system(r), q * r.vec)
end

function Base.:+(r::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(r) === root_system(r2) "parent root system mismatch"

  return DualRootSpaceElem(root_system(r), r.vec + r2.vec)
end

function Base.:-(r::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(r) === root_system(r2) "parent root system mismatch"

  return DualRootSpaceElem(root_system(r), r.vec - r2.vec)
end

function Base.:-(r::DualRootSpaceElem)
  return DualRootSpaceElem(root_system(r), -r.vec)
end

function zero!(r::DualRootSpaceElem)
  r.vec = zero!(r.vec)
  return r
end

function add!(rr::DualRootSpaceElem, r1::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  rr.vec = add!(rr.vec, r1.vec, r2.vec)
  return rr
end

function neg!(rr::DualRootSpaceElem, r::DualRootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  rr.vec = neg!(rr.vec, r.vec)
  return rr
end

function sub!(rr::DualRootSpaceElem, r1::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  rr.vec = sub!(rr.vec, r1.vec, r2.vec)
  return rr
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

function is_coroot(r::DualRootSpaceElem)
  return is_positive_coroot(r) || is_negative_coroot(r)
end

function is_coroot_with_index(r::DualRootSpaceElem)
  fl, i = is_positive_coroot_with_index(r)
  if fl
    return true, i
  end
  fl, j = is_negative_coroot_with_index(r)
  if fl
    return true, j + number_of_positive_roots(root_system(r))
  end
  return false, 0
end

function is_positive_coroot(r::DualRootSpaceElem)
  return haskey(root_system(r).positive_coroots_map, coefficients(r))
end

function is_positive_coroot_with_index(r::DualRootSpaceElem)
  i = get(root_system(r).positive_coroots_map, coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_negative_coroot(r::DualRootSpaceElem)
  return haskey(root_system(r).positive_coroots_map, -coefficients(r))
end

function is_negative_coroot_with_index(r::DualRootSpaceElem)
  i = get(root_system(r).positive_coroots_map, -coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

function is_simple_coroot(r::DualRootSpaceElem)
  return is_simple_coroot_with_index(r)[1]
end

function is_simple_coroot_with_index(r::DualRootSpaceElem)
  i = get(root_system(r).positive_roots_map, coefficients(r), nothing)
  if isnothing(i) || i > number_of_simple_roots(root_system(r))
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
#
#   Weight lattice elements
#
###############################################################################

@doc raw"""
    WeightLatticeElem(R::RootSystem, v::Vector{IntegerUnion}) -> WeightLatticeElem

Return the weight defined by the coefficients `v` of the fundamental weights with respect to the root system `R`.
"""
function WeightLatticeElem(R::RootSystem, v::Vector{<:IntegerUnion})
  return WeightLatticeElem(R, matrix(ZZ, rank(R), 1, v))
end

function WeightLatticeElem(R::RootSystem, r::RootSpaceElem)
  @req root_system(r) === R "Root system mismatch"
  coeffs = transpose!(coefficients(r) * cartan_matrix_tr(R))
  @req all(is_integer, coeffs) "RootSpaceElem does not correspond to a weight"
  return WeightLatticeElem(R, matrix(ZZ, coeffs))
end

function WeightLatticeElem(r::RootSpaceElem)
  return WeightLatticeElem(root_system(r), r)
end

function zero(::Type{WeightLatticeElem}, R::RootSystem)
  return WeightLatticeElem(R, zero_matrix(ZZ, rank(R), 1))
end

function zero(r::WeightLatticeElem)
  return zero(WeightLatticeElem, root_system(r))
end

function Base.:*(n::IntegerUnion, w::WeightLatticeElem)
  return WeightLatticeElem(root_system(w), n * w.vec)
end

function Base.:+(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req root_system(w) === root_system(w2) "parent root system mismatch"

  return WeightLatticeElem(root_system(w), w.vec + w2.vec)
end

function Base.:-(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req root_system(w) === root_system(w2) "parent root system mismatch"

  return WeightLatticeElem(root_system(w), w.vec - w2.vec)
end

function Base.:-(w::WeightLatticeElem)
  return WeightLatticeElem(root_system(w), -w.vec)
end

function zero!(w::WeightLatticeElem)
  w.vec = zero!(w.vec)
  return w
end

function add!(wr::WeightLatticeElem, w1::WeightLatticeElem, w2::WeightLatticeElem)
  @req root_system(wr) === root_system(w1) === root_system(w2) "parent root system mismatch"
  wr.vec = add!(wr.vec, w1.vec, w2.vec)
  return wr
end

function neg!(wr::WeightLatticeElem, w::WeightLatticeElem)
  @req root_system(wr) === root_system(w) "parent root system mismatch"
  wr.vec = neg!(wr.vec, w.vec)
  return wr
end

function sub!(wr::WeightLatticeElem, w1::WeightLatticeElem, w2::WeightLatticeElem)
  @req root_system(wr) === root_system(w1) === root_system(w2) "parent root system mismatch"
  wr.vec = sub!(wr.vec, w1.vec, w2.vec)
  return wr
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

Return the coefficient of the `i`-th fundamental weight in `w`.
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

Return whether `w` is zero.
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

Return the unique dominant weight conjugate to `w`.
"""
function conjugate_dominant_weight(w::WeightLatticeElem)
  return conjugate_dominant_weight!(deepcopy(w))
end

function conjugate_dominant_weight!(w::WeightLatticeElem)
  # w will be dominant once all fundamental weights have a positive coefficient,
  # so search for negative coefficients and make them positive by applying the corresponding reflection.
  s = 1
  while s <= rank(root_system(w))
    if w[s] < 0
      reflect!(w, s)
      s = 1
    else
      s += 1
    end
  end

  return w
end

@doc raw"""
    conjugate_dominant_weight_with_elem(w::WeightLatticeElem) -> Tuple{WeightLatticeElem, WeylGroupElem}

Returns the unique dominant weight `dom` conjugate to `w` and a Weyl group element `x`
such that `x*w == dom`.
"""
function conjugate_dominant_weight_with_elem(w::WeightLatticeElem)
  return conjugate_dominant_weight_with_elem!(deepcopy(w))
end

function conjugate_dominant_weight_with_elem!(w::WeightLatticeElem)
  R = root_system(w)

  # determine the Weyl group element taking w to the fundamental chamber
  word = UInt8[]
  #sizehint!(word, count(<(0), coefficients(w))^2)
  s = 1
  while s <= rank(R)
    if w[s] < 0
      push!(word, UInt8(s))
      reflect!(w, s)
      s = 1
    else
      s += 1
    end
  end

  # reversing word means it is in short revlex normal form
  # and it is the element taking original w to new w
  return w, weyl_group(R)(reverse!(word); normalize=false)
end

function dot(w1::WeightLatticeElem, w2::WeightLatticeElem)
  @req root_system(w1) === root_system(w2) "parent root system mismatch"
  R = root_system(w1)

  return dot(
    coefficients(w1),
    cartan_matrix_inv_tr(R) * (_cartan_symmetrizer_mat(R) * coefficients(w2)),
  )
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
    
Return the `w` reflected at the `s`-th simple root.
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
# more functions

function dot(r::RootSpaceElem, w::WeightLatticeElem)
  @req root_system(r) === root_system(w) "parent root system mismatch"
  R = root_system(r)

  return dot(coefficients(r), _cartan_symmetrizer_mat(R), coefficients(w))
end

function dot(w::WeightLatticeElem, r::RootSpaceElem)
  return dot(r, w)
end

# computes the maximum `p` such that `beta - p*alpha` is still a root
# beta is assumed to be a root
function _root_string_length_down(alpha::RootSpaceElem, beta::RootSpaceElem)
  p = 0
  beta_sub_p_alpha = beta - alpha
  while is_root(beta_sub_p_alpha)
    p += 1
    beta_sub_p_alpha = sub!(beta_sub_p_alpha, alpha)
  end
  return p
end

@doc raw"""
    dim_of_simple_module([T = Int], R::RootSystem, hw::WeightLatticeElem) -> T
    dim_of_simple_module([T = Int], R::RootSystem, hw::Vector{<:IntegerUnion}) -> T

Compute the dimension of the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw` using Weyl's dimension formula.
The return value is of type `T`.

# Example
```jldoctest
julia> R = root_system(:B, 2);

julia> dim_of_simple_module(R, [1, 0])
5
```
"""
function dim_of_simple_module(T::Type, R::RootSystem, hw::WeightLatticeElem)
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"
  rho = weyl_vector(R)
  hw_plus_rho = hw + rho
  num = one(ZZ)
  den = one(ZZ)
  for alpha in positive_roots(R)
    num *= ZZ(dot(hw_plus_rho, alpha))
    den *= ZZ(dot(rho, alpha))
  end
  return T(div(num, den))
end

function dim_of_simple_module(T::Type, R::RootSystem, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(T, R, WeightLatticeElem(R, hw))
end

function dim_of_simple_module(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(Int, R, hw)
end

function dim_of_simple_module(R::RootSystem, hw::WeightLatticeElem)
  return dim_of_simple_module(Int, R, hw)
end

@doc raw"""
    dominant_weights([T = WeightLatticeElem,] R::RootSystem, hw::WeightLatticeElem) -> Vector{T}
    dominant_weights([T = WeightLatticeElem,] R::RootSystem, hw::Vector{<:IntegerUnion}) -> Vector{T}

Computes the dominant weights occurring in the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw`,
sorted ascendingly by the total height of roots needed to reach them from `hw`.

When supplying `T = Vector{Int}`, the weights are returned as vectors of integers.

See [MP82](@cite) for details and the implemented algorithm.

# Example
```jldoctest
julia> R = root_system(:B, 3);

julia> dominant_weights(Vector{Int}, R, [3, 0, 1])
7-element Vector{Vector{Int64}}:
 [3, 0, 1]
 [1, 1, 1]
 [2, 0, 1]
 [0, 0, 3]
 [0, 1, 1]
 [1, 0, 1]
 [0, 0, 1]
```
"""
function dominant_weights(R::RootSystem, hw::WeightLatticeElem)
  return dominant_weights(WeightLatticeElem, R, hw)
end

function dominant_weights(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_weights(R, WeightLatticeElem(R, hw))
end

function dominant_weights(::Type{WeightLatticeElem}, R::RootSystem, hw::WeightLatticeElem)
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"

  pos_roots = positive_roots(R)
  pos_roots_w = WeightLatticeElem.(positive_roots(R))

  ws_with_level = Dict(hw => 0)
  todo = [hw]
  while !isempty(todo)
    new_todo = empty(todo)
    for w in todo
      for (alpha, alpha_w) in zip(pos_roots, pos_roots_w)
        w_sub_alpha = w - alpha_w
        if is_dominant(w_sub_alpha) && !haskey(ws_with_level, w_sub_alpha)
          push!(new_todo, w_sub_alpha)
          push!(ws_with_level, w_sub_alpha => ws_with_level[w] + Int(height(alpha)))
        end
      end
    end
    todo = new_todo
  end
  return first.(sort!(collect(ws_with_level); by=last)) # order by level needed for dominant_character
end

function dominant_weights(T::Type, R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_weights(T, R, WeightLatticeElem(R, hw))
end

function dominant_weights(
  T::Type{<:Vector{<:IntegerUnion}}, R::RootSystem, hw::WeightLatticeElem
)
  weights = dominant_weights(WeightLatticeElem, R, hw)
  return [T(_vec(coefficients(w))) for w in weights]
end

function _action_matrices_on_weights(W::WeylGroup)
  R = root_system(W)
  return map(1:rank(R)) do i
    x = gen(W, i)
    transpose!(
      matrix(
        ZZ, reduce(hcat, coefficients(x * fundamental_weight(R, j)) for j in 1:rank(R))
      ),
    )
  end
end

@doc raw"""
    dominant_character(R::RootSystem, hw::WeightLatticeElem) -> Dict{Vector{Int}, Int}
    dominant_character(R::RootSystem, hw::Vector{<:IntegerUnion}) -> Dict{Vector{Int}, Int}}

Computes the dominant weights occurring in the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw`, together with their multiplicities.

The return type may change in the future.

This function uses an optimized version of the Freudenthal formula, see [MP82](@cite) for details.

# Example
```jldoctest
julia> R = root_system(:B, 3);

julia> dominant_character(R, [2, 0, 1])
Dict{Vector{Int64}, Int64} with 4 entries:
  [1, 0, 1] => 3
  [0, 0, 1] => 6
  [2, 0, 1] => 1
  [0, 1, 1] => 1
```
"""
function dominant_character(R::RootSystem, hw::WeightLatticeElem)
  T = Int
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"
  W = weyl_group(R)
  rho = weyl_vector(R)
  hw_plus_rho = hw + rho
  dot_how_plus_rho = dot(hw_plus_rho, hw_plus_rho)

  pos_roots = positive_roots(R)
  pos_roots_w = WeightLatticeElem.(positive_roots(R))
  pos_roots_w_coeffs = transpose.(coefficients.(pos_roots_w))

  char = Dict(hw => T(1))

  todo = dominant_weights(R, hw)
  all_orbs = Dict{Vector{Int},Vector{Tuple{WeightLatticeElem,Int}}}()
  action_matrices_on_weights = _action_matrices_on_weights(W)

  for w in Iterators.drop(todo, 1)
    stab_inds = [i for (i, ci) in enumerate(coefficients(w)) if iszero(ci)]
    orbs = get!(all_orbs, stab_inds) do
      gens = action_matrices_on_weights[stab_inds]
      push!(gens, -identity_matrix(ZZ, rank(R)))
      G = matrix_group(gens)
      O = orbits(gset(G, *, [pos_roots_w_coeffs; -pos_roots_w_coeffs]))
      [
        (
          WeightLatticeElem(
            R,
            transpose(first(intersect(elements(o), pos_roots_w_coeffs))),
          ),
          length(o),
        ) for o in O
      ]
    end

    accum = sum(
      data -> begin
        rep, len = data
        accum2 = 0
        w_plus_i_rep = w + rep
        while true
          w_plus_i_rep_conj = conjugate_dominant_weight(w_plus_i_rep)
          haskey(char, w_plus_i_rep_conj) || break
          accum2 += char[w_plus_i_rep_conj] * dot(w_plus_i_rep, rep)
          add!(w_plus_i_rep, rep)
        end
        len * accum2
      end, orbs; init=zero(QQ))
    if !iszero(accum)
      w_plus_rho = w + rho
      denom = dot_how_plus_rho - dot(w_plus_rho, w_plus_rho)
      if !iszero(denom)
        char[w] = T(ZZ(div(accum, denom)))
      end
    end
  end
  # return char
  return Dict(Int.(_vec(coefficients(w))) => m for (w, m) in char)
end

function dominant_character(R::RootSystem, hw::Vector{<:IntegerUnion})
  return dominant_character(R, WeightLatticeElem(R, hw))
end

@doc raw"""
    character(R::RootSystem, hw::WeightLatticeElem) -> Vector{T}
    character(R::RootSystem, hw::Vector{<:IntegerUnion}) -> Vector{T}

Computes all weights occurring in the simple module of the Lie algebra defined by the root system `R`
with highest weight `hw`, together with their multiplicities.
This is achieved by acting with the Weyl group on the [`dominant_character`](@ref dominant_character(::RootSystem, ::WeightLatticeElem)).

The return type may change in the future.

# Example
```jldoctest
julia> R = root_system(:B, 3);

julia> character(R, [0, 0, 1])
Dict{Vector{Int64}, Int64} with 8 entries:
  [0, 1, -1]  => 1
  [-1, 1, -1] => 1
  [0, 0, 1]   => 1
  [1, -1, 1]  => 1
  [-1, 0, 1]  => 1
  [0, -1, 1]  => 1
  [0, 0, -1]  => 1
  [1, 0, -1]  => 1
```
"""
function character(R::RootSystem, hw::WeightLatticeElem)
  T = Int
  @req root_system(hw) === R "parent root system mismatch"
  @req is_dominant(hw) "not a dominant weight"
  dom_char = dominant_character(R, hw)
  char = Dict{WeightLatticeElem,T}()

  for (w_, m) in dom_char
    w = WeightLatticeElem(R, w_)
    for w_conj in weyl_orbit(w)
      push!(char, w_conj => m)
    end
  end

  # return char
  return Dict(Int.(_vec(coefficients(w))) => m for (w, m) in char)
end

function character(R::RootSystem, hw::Vector{<:IntegerUnion})
  return character(R, WeightLatticeElem(R, hw))
end

@doc raw"""
    tensor_product_decomposition(R::RootSystem, hw1::WeightLatticeElem, hw2::WeightLatticeElem) -> MSet{Vector{Int}}
    tensor_product_decomposition(R::RootSystem, hw1::Vector{<:IntegerUnion}, hw2::Vector{<:IntegerUnion}) -> MSet{Vector{Int}}

Computes the decomposition of the tensor product of the simple modules of the Lie algebra defined by the root system `R`
with highest weights `hw1` and `hw2` into simple modules with their multiplicities.
This function uses Klymik's formula.

The return type may change in the future.

# Example
```jldoctest
julia> R = root_system(:B, 2);

julia> tensor_product_decomposition(R, [1, 0], [0, 1])
MSet{Vector{Int64}} with 2 elements:
  [1, 1]
  [0, 1]

julia> tensor_product_decomposition(R, [1, 1], [1, 1])
MSet{Vector{Int64}} with 10 elements:
  [0, 0]
  [0, 4]
  [0, 2] : 2
  [2, 0]
  [1, 0]
  [2, 2]
  [3, 0]
  [1, 2] : 2
```
"""
function tensor_product_decomposition(
  R::RootSystem, hw1::WeightLatticeElem, hw2::WeightLatticeElem
)
  @req root_system(hw1) === R "parent root system mismatch"
  @req root_system(hw2) === R "parent root system mismatch"
  @req is_dominant(hw1) "not a dominant weight"
  @req is_dominant(hw2) "not a dominant weight"

  rho = weyl_vector(R)
  hw2_plus_rho = hw2 + rho

  mults = multiset(WeightLatticeElem)
  for (w_, m) in dominant_character(R, hw1)
    for w in weyl_orbit(WeightLatticeElem(R, w_))
      add!(w, hw2_plus_rho)
      w_dom, x = conjugate_dominant_weight_with_elem!(w)
      if all(!iszero, coefficients(w_dom))
        sub!(w_dom, rho)
        coeff = m * (-1)^length(x)
        push!(mults, w_dom, coeff)
      end
    end
  end

  # return mults
  return multiset(
    Dict(Int.(_vec(coefficients(w))) => multiplicity(mults, w) for w in unique(mults))
  )
end

function tensor_product_decomposition(
  R::RootSystem, hw1::Vector{<:IntegerUnion}, hw2::Vector{<:IntegerUnion}
)
  return tensor_product_decomposition(
    R, WeightLatticeElem(R, hw1), WeightLatticeElem(R, hw2)
  )
end

###############################################################################
# internal helpers

# cartan matrix in the format <a^v, b>
function positive_roots_and_reflections(cartan_matrix::ZZMatrix)
  rank, _ = size(cartan_matrix)

  roots = [[l == s ? one(ZZ) : zero(ZZ) for l in 1:rank] for s in 1:rank]
  coroots = [[l == s ? one(ZZ) : zero(ZZ) for l in 1:rank] for s in 1:rank]
  rootidx = Dict(roots[s] => s for s in 1:rank)
  refl = Dict((s, s) => 0 for s in 1:rank)

  i = 1
  while i <= length(roots)
    for s in 1:rank
      if haskey(refl, (s, i))
        continue
      end

      # pairing = dot(roots[i], view(cartan_matrix, s, :)) # currently the below is faster
      pairing = only(view(cartan_matrix, s:s, :) * roots[i])
      # copairing = dot(coroots[i], view(cartan_matrix, :, s)) # currently the below is faster
      copairing = only(coroots[i] * view(cartan_matrix, :, s:s))

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
