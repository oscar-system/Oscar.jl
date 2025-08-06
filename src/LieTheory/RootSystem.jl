###############################################################################
#
#   Root Systems
#
###############################################################################

@doc raw"""
    root_system(cartan_matrix::ZZMatrix; check::Bool=true, detect_type::Bool=true) -> RootSystem
    root_system(cartan_matrix::Matrix{<:IntegerUnion}; check::Bool=true, detect_type::Bool=true) -> RootSystem

Construct the root system defined by the given (generalized) Cartan matrix.

If `check=true` the function will verify that `cartan_matrix` is indeed a generalized Cartan matrix.
Passing `detect_type=false` will skip the detection of the root system type.

# Examples
```jldoctest
julia> root_system([2 -1; -1 2])
Root system of rank 2
  of type A2

julia> root_system(matrix(ZZ, 2, 2, [2, -1, -1, 2]); detect_type=false)
Root system of rank 2
  of unknown type

julia> root_system(matrix(ZZ, [2 -1 -2; -1 2 0; -1 0 2]))
Root system of rank 3
  of type C3 (with non-canonical ordering of simple roots)
```
"""
function root_system(cartan_matrix::ZZMatrix; check::Bool=true, detect_type::Bool=true)
  return RootSystem(cartan_matrix; check, detect_type)
end

function root_system(cartan_matrix::Matrix{<:IntegerUnion}; kwargs...)
  return root_system(matrix(ZZ, cartan_matrix); kwargs...)
end

@doc raw"""
    root_system(fam::Symbol, rk::Int) -> RootSystem

Construct the root system of the given type.

The input must be a valid Cartan type, see [`is_cartan_type(::Symbol, ::Int)`](@ref).

# Examples
```jldoctest
julia> root_system(:A, 2)
Root system of rank 2
  of type A2
```
"""
function root_system(fam::Symbol, rk::Int)
  cartan = cartan_matrix(fam, rk)
  R = root_system(cartan; check=false, detect_type=false)
  set_root_system_type!(R, [(fam, rk)])
  return R
end

@doc raw"""
    root_system(type::Vector{Tuple{Symbol,Int}}) -> RootSystem
    root_system(type::Tuple{Symbol,Int}...) -> RootSystem

Construct the root system of the given type.

Each element of `type` must be a valid Cartan type, see [`is_cartan_type(::Symbol, ::Int)`](@ref).
The vararg version needs at least one element.

# Examples
```jldoctest
julia> root_system([(:A, 2), (:F, 4)])
Root system of rank 6
  of type A2 x F4

julia> root_system(Tuple{Symbol,Int}[])
Root system of rank 0
  of type []
```
"""
function root_system(type::Vector{Tuple{Symbol,Int}})
  cartan = cartan_matrix(type)
  R = root_system(cartan; check=false, detect_type=false)
  set_root_system_type!(R, type)
  return R
end

function root_system(type1::Tuple{Symbol,Int}, type::Tuple{Symbol,Int}...)
  return root_system([type1; type...])
end

function Base.show(io::IO, mime::MIME"text/plain", R::RootSystem)
  @show_name(io, R)
  @show_special(io, mime, R)
  io = pretty(io)
  print(io, "Root system")
  print(io, " of rank ", rank(R))
  println(io, Indent())
  if has_root_system_type(R)
    type, ord = root_system_type_with_ordering(R)
    print(io, "of type ", _root_system_type_string(type))
    if !issorted(ord)
      print(io, " (with non-canonical ordering of simple roots)")
    end
  else
    print(io, "of unknown type")
  end
  print(io, Dedent())
end

function Base.show(io::IO, R::RootSystem)
  @show_name(io, R)
  @show_special(io, R)
  if is_terse(io)
    print(io, "Root system")
  else
    print(io, "Root system")
    if has_root_system_type(R) &&
      ((type, ord) = root_system_type_with_ordering(R); !isempty(type))
      print(io, " of type ", _root_system_type_string(type))
      if !issorted(ord)
        print(io, " (non-canonical ordering)")
      end
    else
      print(io, " of rank ", rank(R))
    end
  end
end

function _root_system_type_string(type::Vector{Tuple{Symbol,Int}})
  isempty(type) && return "[]"
  return join(
    [string(t[1]) * string(t[2]) for t in type], is_unicode_allowed() ? " × " : " x "
  )
end

# not exported atm
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

function Base.hash(R::RootSystem, h::UInt)
  # even though we don't have a == method for RootSystem, we add a hash method
  # to make hashing of RootSpaceElem and WeightLattice more deterministic
  b = 0xeb5362118dea2a0e % UInt
  h = hash(cartan_matrix(R), h)
  return xor(b, h)
end

@doc raw"""
    is_simple(R::RootSystem) -> Bool

Check if `R` is a simple root system.

!!! warning
    Currently only root systems of finite type are supported.
"""
function is_simple(R::RootSystem)
  if is_finite(weyl_group(R))
    return length(root_system_type(R)) == 1
  end
  error("Not implemented") # TODO: implement is_simple
end

@doc raw"""
    number_of_positive_roots(R::RootSystem) -> Int

Return the number of positive roots of `R`. This is the same as the number of negative roots.

See also: [`positive_roots(::RootSystem)`](@ref positive_roots), [`negative_roots(::RootSystem)`](@ref negative_roots).
"""
function number_of_positive_roots(R::RootSystem)
  return length(positive_roots(R))
end

@doc raw"""
    number_of_roots(R::RootSystem) -> Int

Return the number of roots of `R`.

See also: [`roots(::RootSystem)`](@ref roots).
"""
function number_of_roots(R::RootSystem)
  return 2 * number_of_positive_roots(R)
end

@doc raw"""
    number_of_simple_roots(R::RootSystem) -> Int

Return the number of simple roots of `R`.

See also: [`simple_roots(::RootSystem)`](@ref simple_roots).
"""
function number_of_simple_roots(R::RootSystem)
  return rank(R)
end

@doc raw"""
    rank(R::RootSystem) -> Int

Return the rank of `R`, i.e., the number of simple roots.

See also: [`number_of_simple_roots(::RootSystem)`](@ref number_of_simple_roots).
"""
function rank(R::RootSystem)
  return nrows(cartan_matrix(R))
end

@doc raw"""
    root_system_type(R::RootSystem) -> Vector{Tuple{Symbol,Int}}

Return the Cartan type of `R`.

If the type is already known, it is returned directly.
This can be checked with [`has_root_system_type(::RootSystem)`](@ref).

If the type is not known, it is determined and stored in `R`.

See also: [`root_system_type_with_ordering(::RootSystem)`](@ref).

!!! warning
    This function will error if the type is not known yet and the Weyl group is infinite.
"""
function root_system_type(R::RootSystem)
  assure_root_system_type(R)
  return R.type
end

@doc raw"""
    root_system_type_with_ordering(R::RootSystem) -> Vector{Tuple{Symbol,Int}}, Vector{Int}

Return the Cartan type of `R`, together with a vector `ordering` such that [`simple_roots`](@ref)`(R)[ordering]` is a canonical ordering of the simple roots.

If the type is already known, it is returned directly.
This can be checked with [`has_root_system_type(::RootSystem)`](@ref).

If the type is not known, it is determined and stored in `R`.

See also: [`root_system_type(::RootSystem)`](@ref), [`cartan_type_with_ordering(::ZZMatrix)`](@ref).

!!! warning
    This function will error if the type is not known yet and the Weyl group is infinite.
"""
function root_system_type_with_ordering(R::RootSystem)
  assure_root_system_type(R)
  return R.type, R.type_ordering
end

@doc raw"""
    has_root_system_type(R::RootSystem) -> Bool

Check if the root system `R` already knows its Cartan type.

The type can then be queried with [`root_system_type(::RootSystem)`](@ref)
and [`root_system_type_with_ordering(::RootSystem)`](@ref).
"""
function has_root_system_type(R::RootSystem)
  return isdefined(R, :type) && isdefined(R, :type_ordering)
end

function assure_root_system_type(R::RootSystem)
  has_root_system_type(R) && return nothing
  @req is_finite(weyl_group(R)) "Root system type cannot be determined for infinite Weyl groups"
  set_root_system_type!(R, cartan_type_with_ordering(cartan_matrix(R); check=false)...)
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

@doc raw"""
    weight_lattice(R::RootSystem) -> WeightLattice

Return the weight lattice of `R`, i.e. the lattice spanned by the fundamental weights.

This is the parent of all weights of `R`.

# Examples
```jldoctest
julia> weight_lattice(root_system([2 -1; -1 2]))
Weight lattice
  of root system of rank 2
    of type A2

julia> weight_lattice(root_system(matrix(ZZ, 2, 2, [2, -1, -1, 2]); detect_type=false))
Weight lattice
  of root system of rank 2
    of unknown type

julia> weight_lattice(root_system(matrix(ZZ, [2 -1 -2; -1 2 0; -1 0 2])))
Weight lattice
  of root system of rank 3
    of type C3 (with non-canonical ordering of simple roots)
```
"""
function weight_lattice(R::RootSystem)
  return R.weight_lattice::WeightLattice
end

@doc raw"""
    weyl_group(R::RootSystem) -> WeylGroup

Return the Weyl group of `R`.

# Examples
```jldoctest
julia> weyl_group(root_system([2 -1; -1 2]))
Weyl group
  of root system of rank 2
    of type A2

julia> weyl_group(root_system(matrix(ZZ, 2, 2, [2, -1, -1, 2]); detect_type=false))
Weyl group
  of root system of rank 2
    of unknown type

julia> weyl_group(root_system(matrix(ZZ, [2 -1 -2; -1 2 0; -1 0 2])))
Weyl group
  of root system of rank 3
    of type C3 (with non-canonical ordering of simple roots)
```
"""
function weyl_group(R::RootSystem)
  return R.weyl_group::WeylGroup
end

###############################################################################
# root constructors

@doc raw"""
    root(R::RootSystem, i::Int) -> RootSpaceElem

Return the `i`-th root of `R`.

This is a more efficient version for `roots(R)[i]`.

See also: [`roots(::RootSystem)`](@ref).

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

Return the roots of `R`, starting with the positive roots and then the negative roots,
in the order of [`positive_roots(::RootSystem)`](@ref positive_roots) and [`negative_roots(::RootSystem)`](@ref negative_roots).

See also: [`root(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function roots(R::RootSystem)
  return [[r for r in positive_roots(R)]; [-r for r in positive_roots(R)]]
end

@doc raw"""
    positive_root(R::RootSystem, i::Int) -> RootSpaceElem

Return the `i`-th positive root of `R`.

This is a more efficient version for `positive_roots(R)[i]`.

See also: [`positive_roots(::RootSystem)`](@ref).

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

Return the positive roots of `R`, starting with the simple roots in the order
of [`simple_roots(::RootSystem)`](@ref simple_roots), and then increasing in height.

See also: [`positive_root(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.

# Examples
```jldoctest
julia> positive_roots(root_system(:A, 2))
3-element Vector{RootSpaceElem}:
 a_1
 a_2
 a_1 + a_2
```
"""
function positive_roots(R::RootSystem)
  return R.positive_roots::Vector{RootSpaceElem}
end

@doc raw"""
    negative_root(R::RootSystem, i::Int) -> RootSpaceElem

Return the `i`-th negative root of `R`, i.e. the negative of the `i`-th positive root.

This is a more efficient version for `negative_roots(R)[i]`.

See also: [`negative_roots(::RootSystem)`](@ref).

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

Return the negative roots of `R`.

The $i$-th element of the returned vector is the negative root corresponding to the $i$-th positive root.

See also: [`positive_root(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.

# Examples
```jldoctest
julia> negative_roots(root_system(:A, 2))
3-element Vector{RootSpaceElem}:
 -a_1
 -a_2
 -a_1 - a_2
```
"""
function negative_roots(R::RootSystem)
  return [-r for r in positive_roots(R)]
end

@doc raw"""
    simple_root(R::RootSystem, i::Int) -> RootSpaceElem

Return the `i`-th simple root of `R`.

This is a more efficient version for `simple_roots(R)[i]`.

See also: [`simple_roots(::RootSystem)`](@ref).

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

Return the simple roots of `R`.

See also: [`simple_root(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function simple_roots(R::RootSystem)
  return positive_roots(R)[1:rank(R)]
end

@doc raw"""
    highest_root(R::RootSystem) -> RootSpaceElem

Return the highest root of `R`, where `R` is a simple root system of finite type.

This is equivalent to `positive_root(R, number_of_positive_roots(R))`.

See also: [`positive_root(::RootSystem, ::Int)`](@ref positive_root).

# Examples
```jldoctest
julia> highest_root(root_system(:F, 4))
2*a_1 + 3*a_2 + 4*a_3 + 2*a_4
```
"""
function highest_root(R::RootSystem)
  @req is_finite(weyl_group(R)) "Root system is not of finite type"
  @req is_simple(R) "Root system is not simple"
  return positive_root(R, number_of_positive_roots(R))
end

###############################################################################
# coroot constructors

@doc raw"""
    coroot(R::RootSystem, i::Int) -> DualRootSpaceElem

Return the coroot of the `i`-th root of `R`.

This is a more efficient version for `coroots(R)[i]`.

See also: [`coroots(::RootSystem)`](@ref).

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
    coroots(R::RootSystem) -> Vector{DualRootSpaceElem}

Return the coroots of `R` in the order of [`roots(::RootSystem)`](@ref roots),
i.e. starting with the coroots of positive roots and then those of negative roots,
each in the order of [`positive_coroots(::RootSystem)`](@ref positive_coroots) and [`negative_coroots(::RootSystem)`](@ref negative_coroots), respectively.

See also: [`coroot(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function coroots(R::RootSystem)
  return [[r for r in positive_coroots(R)]; [-r for r in positive_coroots(R)]]
end

@doc raw"""
    positive_coroot(R::RootSystem, i::Int) -> DualRootSpaceElem

Return the coroot of the `i`-th positive root of `R`.

This is a more efficient version for `positive_coroots(R)[i]`.

See also: [`positive_coroots(::RootSystem)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function positive_coroot(R::RootSystem, i::Int)
  return (R.positive_coroots::Vector{DualRootSpaceElem})[i]
end

@doc raw"""
    positive_coroots(R::RootSystem) -> DualRootSpaceElem

Return the coroots corresponding to the positive roots of `R`,
in the order of [`positive_roots(::RootSystem)`](@ref positive_roots).

See also: [`positive_coroot(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.

# Examples
```jldoctest
julia> positive_coroots(root_system(:A, 2))
3-element Vector{DualRootSpaceElem}:
 a_1^v
 a_2^v
 a_1^v + a_2^v
```
"""
function positive_coroots(R::RootSystem)
  return R.positive_coroots::Vector{DualRootSpaceElem}
end

@doc raw"""
    negative_coroot(R::RootSystem, i::Int) -> DualRootSpaceElem

Return the coroot of the `i`-th negative root of `R`.

This is a more efficient version for `negative_coroots(R)[i]`.

See also: [`negative_coroots(::RootSystem)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function negative_coroot(R::RootSystem, i::Int)
  return -(R.positive_coroots::Vector{DualRootSpaceElem})[i]
end

@doc raw"""
    negative_coroots(R::RootSystem) -> DualRootSpaceElem

Return the coroots corresponding to the negative roots of `R`,
in the order of [`negative_roots(::RootSystem)`](@ref negative_roots).

See also: [`negative_coroot(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.

# Examples
```jldoctest
julia> negative_coroots(root_system(:A, 2))
3-element Vector{DualRootSpaceElem}:
 -a_1^v
 -a_2^v
 -a_1^v - a_2^v
```
"""
function negative_coroots(R::RootSystem)
  return [-r for r in positive_coroots(R)]
end

@doc raw"""
    simple_coroot(R::RootSystem, i::Int) -> DualRootSpaceElem

Return the coroot of the `i`-th simple root of `R`.

This is a more efficient version for `simple_coroots(R)[i]`.

See also: [`simple_coroots(::RootSystem)`](@ref).

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
    simple_coroots(R::RootSystem) -> DualRootSpaceElem

Return the coroots corresponding to the simple roots of `R`,
in the order of [`simple_roots(::RootSystem)`](@ref simple_roots).

See also: [`simple_coroot(::RootSystem, ::Int)`](@ref).

!!! note
    This function does not return a copy of the asked for object,
    but the internal field of the root system.
    Mutating the returned object will lead to undefined behavior.
"""
function simple_coroots(R::RootSystem)
  return positive_coroots(R)[1:rank(R)]
end

###############################################################################
# weight constructors

@doc raw"""
    fundamental_weight(R::RootSystem, i::Int) -> WeightLatticeElem

Return the `i`-th fundamental weight of `R`.

This is a more efficient version for `fundamental_weights(R)[i]`.

See also: [`fundamental_weight(::RootSystem)`](@ref).
"""
function fundamental_weight(R::RootSystem, i::Int)
  return gen(weight_lattice(R), i)
end

@doc raw"""
    fundamental_weights(R::RootSystem) -> Vector{WeightLatticeElem}

Return the fundamental weights corresponding to the simple roots of `R`,
in the order of [`simple_roots(::RootSystem)`](@ref simple_roots).

See also: [`fundamental_weight(::RootSystem, ::Int)`](@ref).

# Examples
```jldoctest
julia> fundamental_weights(root_system(:A, 2))
2-element Vector{WeightLatticeElem}:
 w_1
 w_2
```
"""
function fundamental_weights(R::RootSystem)
  return gens(weight_lattice(R))
end

@doc raw"""
    weyl_vector(R::RootSystem) -> WeightLatticeElem

Return the Weyl vector $\rho$ of `R`,
that is the sum of all fundamental weights or,
equivalently, half the sum of all positive roots.
"""
function weyl_vector(R::RootSystem)
  return WeightLatticeElem(R, matrix(ZZ, 1, rank(R), fill(1, rank(R))))
end

###############################################################################
#
#   Root space elements
#
###############################################################################

@doc raw"""
    RootSpaceElem(R::RootSystem, vec::Vector{<:RationalUnion}) -> RootSpaceElem

Construct a root space element in the root system `R` with the given coefficients w.r.t. the simple roots of `R`.
"""
function RootSpaceElem(root_system::RootSystem, vec::Vector{<:RationalUnion})
  return RootSpaceElem(root_system, matrix(QQ, 1, length(vec), vec))
end

@doc raw"""
    RootSpaceElem(w::WeightLatticeElem) -> RootSpaceElem

Construct a root space element from the weight lattice element `w`.
"""
function RootSpaceElem(w::WeightLatticeElem)
  R = root_system(w)
  coeffs = coefficients(w) * cartan_matrix_inv_tr(R)
  return RootSpaceElem(R, matrix(QQ, coeffs))
end

@doc raw"""
    zero(::Type{RootSpaceElem}, R::RootSystem) -> RootSpaceElem

Return the neutral additive element in the root space of `R`.
"""
function zero(::Type{RootSpaceElem}, R::RootSystem)
  return RootSpaceElem(R, zero_matrix(QQ, 1, rank(R)))
end

function zero(r::RootSpaceElem)
  return zero(RootSpaceElem, root_system(r))
end

function Base.:*(q::RationalUnion, r::RootSpaceElem)
  return RootSpaceElem(root_system(r), q * r.vec)
end

function Base.:*(r::RootSpaceElem, q::RationalUnion)
  return RootSpaceElem(root_system(r), r.vec * q)
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
  return RootSpaceElem(root_system(r), zero!(r.vec))
end

function neg!(rr::RootSpaceElem, r::RootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), neg!(rr.vec, r.vec))
end

function add!(rr::RootSpaceElem, r1::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), add!(rr.vec, r1.vec, r2.vec))
end

function sub!(rr::RootSpaceElem, r1::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), sub!(rr.vec, r1.vec, r2.vec))
end

function mul!(rr::RootSpaceElem, r::RootSpaceElem, q::RationalUnionOrPtr)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), mul!(rr.vec, r.vec, q))
end

function mul!(rr::RootSpaceElem, q::RationalUnionOrPtr, r::RootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), mul!(rr.vec, q, r.vec))
end

function addmul!(rr::RootSpaceElem, r::RootSpaceElem, q::RationalUnionOrPtr)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), addmul!(rr.vec, r.vec, q))
end

function addmul!(rr::RootSpaceElem, q::RationalUnionOrPtr, r::RootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), addmul!(rr.vec, q, r.vec))
end

# ignore temp storage
addmul!(rr::RootSpaceElem, r::RootSpaceElem, q::RationalUnionOrPtr, t) = addmul!(rr, r, q)
addmul!(rr::RootSpaceElem, q::RationalUnionOrPtr, r::RootSpaceElem, t) = addmul!(rr, q, r)

function submul!(rr::RootSpaceElem, r::RootSpaceElem, q::RationalUnionOrPtr)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), submul!(rr.vec, r.vec, q))
end

function submul!(rr::RootSpaceElem, q::RationalUnionOrPtr, r::RootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return RootSpaceElem(root_system(rr), submul!(rr.vec, q, r.vec))
end

# ignore temp storage
submul!(rr::RootSpaceElem, r::RootSpaceElem, q::RationalUnionOrPtr, t) = submul!(rr, r, q)
submul!(rr::RootSpaceElem, q::RationalUnionOrPtr, r::RootSpaceElem, t) = submul!(rr, q, r)

function mat_entry_ptr(r::RootSpaceElem, i::Int)
  return mat_entry_ptr(r.vec, 1, i)
end

function is_zero_entry(r::RootSpaceElem, i::Int)
  return is_zero_entry(r.vec, 1, i)
end

function is_positive_entry(r::RootSpaceElem, i::Int)
  return is_positive_entry(r.vec, 1, i)
end

function is_negative_entry(r::RootSpaceElem, i::Int)
  return is_negative_entry(r.vec, 1, i)
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

function Base.hash(r::RootSpaceElem, h::UInt)
  b = 0xbe7603eb38c985ad % UInt
  h = hash(root_system(r), h)
  h = hash(r.vec, h)
  return xor(b, h)
end

@doc raw"""
    coefficients(r::RootSpaceElem) -> QQMatrix

Return the coefficients of the root space element `r`
w.r.t. the simple roots as a row vector.

!!! note
    The return type may not be relied on;
    we only guarantee that it is a one-dimensional iterable with `eltype` `QQFieldElem`
    that can be indexed with integers.
"""
function coefficients(r::RootSpaceElem)
  return r.vec
end

@doc raw"""
    coeff(r::RootSpaceElem, i::Int) -> QQFieldElem

Return the coefficient of the `i`-th simple root in `r`.

This can be also accessed via `r[i]`.
"""
function coeff(r::RootSpaceElem, i::Int)
  return r.vec[i]
end

function Base.getindex(r::RootSpaceElem, i::Int)
  return coeff(r, i)
end

function dot(r1::RootSpaceElem, r2::RootSpaceElem)
  @req root_system(r1) === root_system(r2) "parent root system mismatch"

  # return dot(coefficients(r1) * _bilinear_form_QQ(root_system(r1)), coefficients(r2)) # currently the below is faster
  return only(
    coefficients(r1) * _bilinear_form_QQ(root_system(r1)) * transpose(coefficients(r2))
  )
end

function expressify(r::RootSpaceElem; context=nothing)
  if is_unicode_allowed()
    return expressify(r, :α; context)
  else
    return expressify(r, :a; context)
  end
end

function expressify(r::RootSpaceElem, s; context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(r.vec)
    push!(sum.args, Expr(:call, :*, expressify(r.vec[i]; context), "$(s)_$(i)"))
  end
  return sum
end
@enable_all_show_via_expressify RootSpaceElem

@doc raw"""
    height(r::RootSpaceElem) -> QQFieldElem

For a root `r`, returns the height of `r`, i.e. the sum of the coefficients of the simple roots.
If `r` is not a root, the return value is arbitrary.
"""
function height(r::RootSpaceElem)
  return sum(coefficients(r))
end

@doc raw"""
    is_root(r::RootSpaceElem) -> Bool

Check if `r` is a root of its root system.

See also: [`is_root_with_index(::RootSpaceElem)`](@ref).
"""
function is_root(r::RootSpaceElem)
  return is_positive_root(r) || is_negative_root(r)
end

@doc raw"""
    is_root_with_index(r::RootSpaceElem) -> Bool, Int

Check if `r` is a root of its root system
and return this together with the index of the root in [`roots(::RootSystem)`](@ref).

If `r` is not a root, the second return value is arbitrary.

See also: [`is_root(::RootSpaceElem)`](@ref).
"""
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

@doc raw"""
    is_positive_root(r::RootSpaceElem) -> Bool

Check if `r` is a positive root of its root system.

See also: [`is_positive_root_with_index(::RootSpaceElem)`](@ref).
"""
function is_positive_root(r::RootSpaceElem)
  return haskey(root_system(r).positive_roots_map, coefficients(r))
end

@doc raw"""
    is_positive_root_with_index(r::RootSpaceElem) -> Bool, Int

Check if `r` is a positive root of its root system
and return this together with the index of the root in [`positive_roots(::RootSystem)`](@ref).

If `r` is not a positive root, the second return value is arbitrary.

See also: [`is_positive_root(::RootSpaceElem)`](@ref).
"""
function is_positive_root_with_index(r::RootSpaceElem)
  i = get(root_system(r).positive_roots_map, coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

@doc raw"""
    is_negative_root(r::RootSpaceElem) -> Bool

Check if `r` is a negative root of its root system.

See also: [`is_negative_root_with_index(::RootSpaceElem)`](@ref).
"""
function is_negative_root(r::RootSpaceElem)
  return haskey(root_system(r).positive_roots_map, -coefficients(r))
end

@doc raw"""
    is_negative_root_with_index(r::RootSpaceElem) -> Bool, Int

Check if `r` is a negative root of its root system
and return this together with the index of the root in [`negative_roots(::RootSystem)`](@ref).

If `r` is not a negative root, the second return value is arbitrary.

See also: [`is_negative_root(::RootSpaceElem)`](@ref).
"""
function is_negative_root_with_index(r::RootSpaceElem)
  i = get(root_system(r).positive_roots_map, -coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

@doc raw"""
    is_simple_root(r::RootSpaceElem) -> Bool

Check if `r` is a simple root of its root system.

See also: [`is_simple_root_with_index(::RootSpaceElem)`](@ref).
"""
function is_simple_root(r::RootSpaceElem)
  return is_simple_root_with_index(r)[1]
end

@doc raw"""
    is_simple_root_with_index(r::RootSpaceElem) -> Bool, Int

Check if `r` is a simple root of its root system
and return this together with the index of the root in [`simple_roots(::RootSystem)`](@ref).

If `r` is not a simple root, the second return value is arbitrary.

See also: [`is_simple_root(::RootSpaceElem)`](@ref).
"""
function is_simple_root_with_index(r::RootSpaceElem)
  i = get(root_system(r).positive_roots_map, coefficients(r), nothing)
  if isnothing(i) || i > number_of_simple_roots(root_system(r))
    return false, 0
  else
    return true, i
  end
end

@doc raw"""
    iszero(r::RootSpaceElem) -> Bool

Check if `r` is the neutral additive element in the root space of its root system.
"""
function Base.iszero(r::RootSpaceElem)
  return iszero(r.vec)
end

@doc raw"""
    reflect(r::RootSpaceElem, s::Int) -> RootSpaceElem
  
Return the reflection of `r` in the hyperplane orthogonal to the `s`-th simple root.

See also: [`reflect!(::RootSpaceElem, ::Int)`](@ref).
"""
function reflect(r::RootSpaceElem, s::Int)
  return reflect!(deepcopy(r), s)
end

@doc raw"""
    reflect!(r::RootSpaceElem, s::Int) -> RootSpaceElem

Reflect `r` in the hyperplane orthogonal to the `s`-th simple root, and return it.

This is a mutating version of [`reflect(::RootSpaceElem, ::Int)`](@ref).
"""
function reflect!(r::RootSpaceElem, s::Int)
  GC.@preserve r begin
    sub!(mat_entry_ptr(r, s), dot(view(cartan_matrix(root_system(r)), s, :), r.vec))
  end
  return r
end

@doc raw"""
    reflect(r::RootSpaceElem, beta::RootSpaceElem) -> RootSpaceElem
  
Return the reflection of `r` in the hyperplane orthogonal to root `beta`.

See also: [`reflect!(::RootSpaceElem, ::RootSpaceElem)`](@ref).
"""
function reflect(r::RootSpaceElem, beta::RootSpaceElem)
  return reflect!(deepcopy(r), beta)
end

@doc raw"""
    reflect!(r::RootSpaceElem, beta::RootSpaceElem) -> RootSpaceElem

Reflect `r` in the hyperplane orthogonal to the root `beta`, and return it.

This is a mutating version of [`reflect(::RootSpaceElem, ::RootSpaceElem)`](@ref).
"""
function reflect!(r::RootSpaceElem, beta::RootSpaceElem)
  @req root_system(r) === root_system(beta) "Incompatible root systems"
  for s in word(reflection(beta))
    reflect!(r, Int(s))
  end
  return r
end

@doc raw"""
    root_system(r::RootSpaceElem) -> RootSystem

Return the root system `r` belongs to.
"""
function root_system(r::RootSpaceElem)
  return r.root_system
end

ConformanceTests.equality(a::RootSpaceElem, b::RootSpaceElem) = a == b

###############################################################################
#
#   Dual root space elements
#
###############################################################################

@doc raw"""
    DualRootSpaceElem(R::RootSystem, vec::Vector{<:RationalUnion}) -> DualRootSpaceElem

Construct a dual root space element in the root system `R` with the given coefficients w.r.t. the simple coroots of `R`.
"""
function DualRootSpaceElem(root_system::RootSystem, vec::Vector{<:RationalUnion})
  return DualRootSpaceElem(root_system, matrix(QQ, 1, length(vec), vec))
end

@doc raw"""
    zero(::Type{DualRootSpaceElem}, R::RootSystem) -> DualRootSpaceElem

Return the neutral additive element in the dual root space of `R`.
"""
function zero(::Type{DualRootSpaceElem}, R::RootSystem)
  return DualRootSpaceElem(R, zero_matrix(QQ, 1, rank(R)))
end

function zero(r::DualRootSpaceElem)
  return zero(DualRootSpaceElem, root_system(r))
end

function Base.:*(q::RationalUnion, r::DualRootSpaceElem)
  return DualRootSpaceElem(root_system(r), q * r.vec)
end

function Base.:*(r::DualRootSpaceElem, q::RationalUnion)
  return DualRootSpaceElem(root_system(r), r.vec * q)
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
  return DualRootSpaceElem(root_system(r), zero!(r.vec))
end

function add!(rr::DualRootSpaceElem, r1::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), add!(rr.vec, r1.vec, r2.vec))
end

function neg!(rr::DualRootSpaceElem, r::DualRootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), neg!(rr.vec, r.vec))
end

function sub!(rr::DualRootSpaceElem, r1::DualRootSpaceElem, r2::DualRootSpaceElem)
  @req root_system(rr) === root_system(r1) === root_system(r2) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), sub!(rr.vec, r1.vec, r2.vec))
end

function mul!(rr::DualRootSpaceElem, r::DualRootSpaceElem, q::RationalUnionOrPtr)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), mul!(rr.vec, r.vec, q))
end

function mul!(rr::DualRootSpaceElem, q::RationalUnionOrPtr, r::DualRootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), mul!(rr.vec, q, r.vec))
end

function addmul!(rr::DualRootSpaceElem, r::DualRootSpaceElem, q::RationalUnionOrPtr)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), addmul!(rr.vec, r.vec, q))
end

function addmul!(rr::DualRootSpaceElem, q::RationalUnionOrPtr, r::DualRootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), addmul!(rr.vec, q, r.vec))
end

# ignore temp storage
addmul!(rr::DualRootSpaceElem, r::DualRootSpaceElem, q::RationalUnionOrPtr, t) =
  addmul!(rr, r, q)
addmul!(rr::DualRootSpaceElem, q::RationalUnionOrPtr, r::DualRootSpaceElem, t) =
  addmul!(rr, q, r)

function submul!(rr::DualRootSpaceElem, r::DualRootSpaceElem, q::RationalUnionOrPtr)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), submul!(rr.vec, r.vec, q))
end

function submul!(rr::DualRootSpaceElem, q::RationalUnionOrPtr, r::DualRootSpaceElem)
  @req root_system(rr) === root_system(r) "parent root system mismatch"
  return DualRootSpaceElem(root_system(rr), submul!(rr.vec, q, r.vec))
end

# ignore temp storage
submul!(rr::DualRootSpaceElem, r::DualRootSpaceElem, q::RationalUnionOrPtr, t) =
  submul!(rr, r, q)
submul!(rr::DualRootSpaceElem, q::RationalUnionOrPtr, r::DualRootSpaceElem, t) =
  submul!(rr, q, r)

function mat_entry_ptr(r::DualRootSpaceElem, i::Int)
  return mat_entry_ptr(r.vec, 1, i)
end

function is_zero_entry(r::DualRootSpaceElem, i::Int)
  return is_zero_entry(r.vec, 1, i)
end

function is_positive_entry(r::DualRootSpaceElem, i::Int)
  return is_positive_entry(r.vec, 1, i)
end

function is_negative_entry(r::DualRootSpaceElem, i::Int)
  return is_negative_entry(r.vec, 1, i)
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

function Base.hash(r::DualRootSpaceElem, h::UInt)
  b = 0x721bec0418bdbe0f % UInt
  h = hash(root_system(r), h)
  h = hash(r.vec, h)
  return xor(b, h)
end

@doc raw"""
    coefficients(r::DualRootSpaceElem) -> QQMatrix

Return the coefficients of the dual root space element `r`
w.r.t. the simple coroots as a row vector.

!!! note
    The return type may not be relied on;
    we only guarantee that it is a one-dimensional iterable with `eltype` `QQFieldElem`
    that can be indexed with integers.
"""
function coefficients(r::DualRootSpaceElem)
  return r.vec
end

@doc raw"""
    coeff(r::DualRootSpaceElem, i::Int) -> QQFieldElem

Return the coefficient of the `i`-th simple coroot in `r`.

This can be also accessed via `r[i]`.
"""
function coeff(r::DualRootSpaceElem, i::Int)
  return r.vec[i]
end

function Base.getindex(r::DualRootSpaceElem, i::Int)
  return coeff(r, i)
end

function expressify(r::DualRootSpaceElem; context=nothing)
  if is_unicode_allowed()
    return expressify(r, :α̌; context)
  else
    return expressify(r, :a, Symbol("^v"); context)
  end
end

function expressify(r::DualRootSpaceElem, s, s_suffix=Symbol(""); context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(r.vec)
    push!(sum.args, Expr(:call, :*, expressify(r.vec[i]; context), "$(s)_$(i)$(s_suffix)"))
  end
  return sum
end
@enable_all_show_via_expressify DualRootSpaceElem

@doc raw"""
    height(r::DualRootSpaceElem) -> QQFieldElem

For a coroot `r`, returns the height of `r`, i.e. the sum of the coefficients of the simple coroots.
If `r` is not a coroot, the return value is arbitrary.
"""
function height(r::DualRootSpaceElem)
  return sum(coefficients(r))
end

@doc raw"""
    is_coroot(r::DualRootSpaceElem) -> Bool

Check if `r` is a coroot of its root system.

See also: [`is_coroot_with_index(::DualRootSpaceElem)`](@ref).
"""
function is_coroot(r::DualRootSpaceElem)
  return is_positive_coroot(r) || is_negative_coroot(r)
end

@doc raw"""
    is_coroot_with_index(r::DualRootSpaceElem) -> Bool, Int

Check if `r` is a coroot of its root system
and return this together with the index of the coroot in [`coroots(::RootSystem)`](@ref).

If `r` is not a coroot, the second return value is arbitrary.

See also: [`is_coroot(::DualRootSpaceElem)`](@ref).
"""
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

@doc raw"""
    is_positive_coroot(r::DualRootSpaceElem) -> Bool

Check if `r` is a positive coroot of its root system.

See also: [`is_positive_coroot_with_index(::DualRootSpaceElem)`](@ref).
"""
function is_positive_coroot(r::DualRootSpaceElem)
  return haskey(root_system(r).positive_coroots_map, coefficients(r))
end

@doc raw"""
    is_positive_coroot_with_index(r::DualRootSpaceElem) -> Bool, Int

Check if `r` is a positive coroot of its root system
and return this together with the index of the coroot in [`positive_coroots(::RootSystem)`](@ref).

If `r` is not a positive coroot, the second return value is arbitrary.

See also: [`is_positive_coroot(::DualRootSpaceElem)`](@ref).
"""
function is_positive_coroot_with_index(r::DualRootSpaceElem)
  i = get(root_system(r).positive_coroots_map, coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

@doc raw"""
    is_negative_coroot(r::DualRootSpaceElem) -> Bool

Check if `r` is a negative coroot of its root system.

See also: [`is_negative_coroot_with_index(::DualRootSpaceElem)`](@ref).
"""
function is_negative_coroot(r::DualRootSpaceElem)
  return haskey(root_system(r).positive_coroots_map, -coefficients(r))
end

@doc raw"""
    is_negative_coroot_with_index(r::DualRootSpaceElem) -> Bool, Int

Check if `r` is a negative coroot of its root system
and return this together with the index of the coroot in [`negative_coroots(::RootSystem)`](@ref).

If `r` is not a negative coroot, the second return value is arbitrary.

See also: [`is_negative_coroot(::DualRootSpaceElem)`](@ref).
"""
function is_negative_coroot_with_index(r::DualRootSpaceElem)
  i = get(root_system(r).positive_coroots_map, -coefficients(r), nothing)
  if isnothing(i)
    return false, 0
  else
    return true, i
  end
end

@doc raw"""
    is_simple_coroot(r::DualRootSpaceElem) -> Bool

Check if `r` is a simple coroot of its root system.

See also: [`is_simple_coroot_with_index(::DualRootSpaceElem)`](@ref).
"""
function is_simple_coroot(r::DualRootSpaceElem)
  return is_simple_coroot_with_index(r)[1]
end

@doc raw"""
    is_simple_coroot_with_index(r::DualRootSpaceElem) -> Bool, Int

Check if `r` is a simple coroot of its root system
and return this together with the index of the coroot in [`simple_coroots(::RootSystem)`](@ref).

If `r` is not a simple coroot, the second return value is arbitrary.

See also: [`is_simple_coroot(::DualRootSpaceElem)`](@ref).
"""
function is_simple_coroot_with_index(r::DualRootSpaceElem)
  i = get(root_system(r).positive_roots_map, coefficients(r), nothing)
  if isnothing(i) || i > number_of_simple_roots(root_system(r))
    return false, 0
  else
    return true, i
  end
end

@doc raw"""
    iszero(r::DualRootSpaceElem) -> Bool

Check if `r` is the neutral additive element in the dual root space of its root system.
"""
function Base.iszero(r::DualRootSpaceElem)
  return iszero(r.vec)
end

@doc raw"""
    root_system(r::DualRootSpaceElem) -> RootSystem

Return the root system `r` belongs to.
"""
function root_system(r::DualRootSpaceElem)
  return r.root_system
end

ConformanceTests.equality(a::DualRootSpaceElem, b::DualRootSpaceElem) = a == b

###############################################################################
# more functions

function dot(r::RootSpaceElem, w::WeightLatticeElem)
  @req root_system(r) === root_system(w) "parent root system mismatch"
  R = root_system(r)

  return dot(coefficients(r) * _cartan_symmetrizer_mat(R), coefficients(w))
end

function dot(w::WeightLatticeElem, r::RootSpaceElem)
  return dot(r, w)
end

###############################################################################
# internal helpers

@doc raw"""
    _positive_roots_and_reflections(cartan_matrix::ZZMatrix) -> Vector{Vector{ZZRingElem}}, Vector{Vector{ZZRingElem}}, Matrix{UInt64}

Compute the positive roots, the positive coroots, and a matrix `refl` of size $m \times n$,
where $m$ is the rank of the root system and $n$ is the number of minimal roots.

The minimal roots and coroots are given as coefficient vectors w.r.t. the simple roots and simple coroots, respectively.
The minimal roots are indexed by `1:n`, with the first `m` of them corresponding
to the simple roots, and the other roots sorted by height.

If `beta = alpha_j * s_i` is a minimal root, then `refl_table[i, j]` stores the index of beta, and otherwise `0`.
Note that `refl_table[i, i] = 0` for every simple root `alpha_i`.
"""
function _positive_roots_and_reflections(cartan_matrix::ZZMatrix)
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

  return roots[perm], coroots[perm], table
end
