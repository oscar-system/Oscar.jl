###############################################################################
#
#   Weight lattices
#
###############################################################################

function elem_type(::Type{WeightLattice})
  return WeightLatticeElem
end

@doc raw"""
    rank(P::WeightLattice) -> Int

Return the rank of the weight lattice `P`.
"""
function rank(P::WeightLattice)
  return rank(root_system(P))
end

@doc raw"""
    root_system(P::WeightLattice) -> RootSystem

Return the underlying root system of `P`.
"""
function root_system(P::WeightLattice)
  return P.root_system
end

@doc raw"""
    zero(P::WeightLattice) -> WeightLatticeElem

Return the neutral additive element in the weight lattice `P`.
"""
function zero(P::WeightLattice)
  return WeightLatticeElem(P, zero_matrix(ZZ, 1, rank(P)))
end

function Base.show(io::IO, mime::MIME"text/plain", P::WeightLattice)
  @show_name(io, P)
  @show_special(io, mime, P)
  io = pretty(io)
  println(io, "Weight lattice")
  print(io, Indent(), "of ", Lowercase())
  show(io, mime, root_system(P))
  print(io, Dedent())
end

function Base.show(io::IO, P::WeightLattice)
  @show_name(io, P)
  @show_special(io, P)
  io = pretty(io)
  if is_terse(io)
    print(io, "Weight lattice")
  else
    print(io, "Weight lattice of ", Lowercase(), root_system(P))
  end
end

function number_of_generators(P::WeightLattice)
  return rank(P)
end

@doc raw"""
    gen(P::WeightLattice, i::Int) -> WeightLatticeElem

Return the `i`-th generator of the weight lattice `P`,
i.e. the `i`-th fundamental weight of the root system of `P`.

This is a more efficient version for `gens(P)[i]`.

See also: [`fundamental_weight(::RootSystem, ::Int)`](@ref).
"""
function gen(P::WeightLattice, i::Int)
  @req 1 <= i <= rank(P) "invalid index"
  return WeightLatticeElem(P, matrix(ZZ, 1, rank(P), i .== 1:rank(P)))
end

@doc raw"""
    gens(P::WeightLattice) -> Vector{WeightLatticeElem}

Return the generators of the weight lattice `P`,
i.e. the fundamental weights of the root system of `P`.

See also: [`gen(::WeightLattice, ::Int)`](@ref), [`fundamental_weights(::RootSystem)`](@ref).
"""
function gens(P::WeightLattice)
  return [gen(P, i) for i in 1:rank(P)]
end

function is_abelian(P::WeightLattice)
  return true
end

@doc raw"""
    is_finite(P::WeightLattice) -> Bool

Check if the weight lattice `P` is finite, i.e. if it has rank 0.
"""
function is_finite(P::WeightLattice)
  return iszero(rank(P))
end

function Base.hash(P::WeightLattice, h::UInt)
  # even though we don't have a == method for WeightLattice, we add a hash method
  # to make hashing of WeightLattice and WeightLatticeElem more deterministic
  b = 0x770fdd486cbdeea2 % UInt
  h = hash(root_system(P), h)
  return xor(b, h)
end

###############################################################################
#
#   Weight lattice elements
#
###############################################################################

@doc raw"""
    WeightLatticeElem(R::RootSystem, vec::ZZMatrix) -> WeightLatticeElem

Construct a weight lattice element in the root system `R` with the given coefficient vector w.r.t. the fundamental weights of `R`.

`vec` must be a row vector of the same length as the rank of `R`.
"""
function WeightLatticeElem(R::RootSystem, vec::ZZMatrix)
  return WeightLatticeElem(weight_lattice(R), vec)
end

@doc raw"""
    WeightLatticeElem(R::RootSystem, vec::Vector{<:IntegerUnion}) -> WeightLatticeElem

Construct a weight lattice element in the root system `R` with the given coefficients w.r.t. the fundamental weights of `R`.
"""
function WeightLatticeElem(R::RootSystem, v::Vector{<:IntegerUnion})
  return WeightLatticeElem(R, matrix(ZZ, 1, rank(R), v))
end

@doc raw"""
    WeightLatticeElem(P::WeightLattice, vec::Vector{<:IntegerUnion}) -> WeightLatticeElem

Construct a weight lattice element in `P` with the given coefficients w.r.t. the fundamental weights of corresponding root system.
"""
function WeightLatticeElem(P::WeightLattice, v::Vector{<:IntegerUnion})
  return WeightLatticeElem(P, matrix(ZZ, 1, rank(P), v))
end

@doc raw"""
    WeightLatticeElem(r::RootSpaceElem) -> WeightLatticeElem

Construct a weight lattice element from the root space element `r`.
"""
function WeightLatticeElem(r::RootSpaceElem)
  R = root_system(r)
  coeffs = coefficients(r) * cartan_matrix_tr(R)
  @req all(is_integer, coeffs) "RootSpaceElem does not correspond to a weight"
  return WeightLatticeElem(R, matrix(ZZ, coeffs))
end

function parent_type(::Type{WeightLatticeElem})
  return WeightLattice
end

function parent(w::WeightLatticeElem)
  return w.parent_lat
end

function root_system(w::WeightLatticeElem)
  return root_system(parent(w))
end

function zero(w::WeightLatticeElem)
  return zero(parent(w))
end

function Base.:*(n::IntegerUnion, w::WeightLatticeElem)
  return WeightLatticeElem(parent(w), n * w.vec)
end

function Base.:*(w::WeightLatticeElem, n::IntegerUnion)
  return WeightLatticeElem(parent(w), w.vec * n)
end

function Base.:+(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req parent(w) === parent(w2) "parent mismatch"

  return WeightLatticeElem(parent(w), w.vec + w2.vec)
end

function Base.:-(w::WeightLatticeElem, w2::WeightLatticeElem)
  @req parent(w) === parent(w2) "parent mismatch"

  return WeightLatticeElem(parent(w), w.vec - w2.vec)
end

function Base.:-(w::WeightLatticeElem)
  return WeightLatticeElem(parent(w), -w.vec)
end

function zero!(w::WeightLatticeElem)
  w.vec = zero!(w.vec)
  return w
end

function add!(wr::WeightLatticeElem, w1::WeightLatticeElem, w2::WeightLatticeElem)
  @req parent(wr) === parent(w1) === parent(w2) "parent mismatch"
  wr.vec = add!(wr.vec, w1.vec, w2.vec)
  return wr
end

function neg!(wr::WeightLatticeElem, w::WeightLatticeElem)
  @req parent(wr) === parent(w) "parent mismatch"
  wr.vec = neg!(wr.vec, w.vec)
  return wr
end

function sub!(wr::WeightLatticeElem, w1::WeightLatticeElem, w2::WeightLatticeElem)
  @req parent(wr) === parent(w1) === parent(w2) "parent mismatch"
  wr.vec = sub!(wr.vec, w1.vec, w2.vec)
  return wr
end

function mul!(wr::WeightLatticeElem, w::WeightLatticeElem, n::IntegerUnion)
  @req parent(wr) === parent(w) "parent mismatch"
  wr.vec = mul!(wr.vec, w.vec, n)
  return wr
end

function mul!(wr::WeightLatticeElem, n::IntegerUnion, w::WeightLatticeElem)
  @req parent(wr) === parent(w) "parent mismatch"
  wr.vec = mul!(wr.vec, n, w.vec)
  return wr
end

function addmul!(wr::WeightLatticeElem, w::WeightLatticeElem, n::IntegerUnion)
  @req parent(wr) === parent(w) "parent mismatch"
  wr.vec = addmul!(wr.vec, w.vec, n)
  return wr
end

function addmul!(wr::WeightLatticeElem, n::IntegerUnion, w::WeightLatticeElem)
  @req parent(wr) === parent(w) "parent mismatch"
  wr.vec = addmul!(wr.vec, n, w.vec)
  return wr
end

# ignore temp storage
addmul!(wr::WeightLatticeElem, w::WeightLatticeElem, n::IntegerUnion, t) = addmul!(wr, w, n)
addmul!(wr::WeightLatticeElem, n::IntegerUnion, w::WeightLatticeElem, t) = addmul!(wr, n, w)

function Base.:(==)(w1::WeightLatticeElem, w2::WeightLatticeElem)
  return parent(w1) === parent(w2) && w1.vec == w2.vec
end

function Base.deepcopy_internal(w::WeightLatticeElem, dict::IdDict)
  if haskey(dict, w)
    return dict[w]
  end

  w2 = WeightLatticeElem(parent(w), deepcopy_internal(w.vec, dict))
  dict[w] = w2
  return w2
end

function Base.hash(w::WeightLatticeElem, h::UInt)
  b = 0x7b2fefadacf46f4e % UInt
  h = hash(parent(w), h)
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

@doc raw"""
    coefficients(w::WeightLatticeElem) -> ZZMatrix

Return the coefficients of the weight lattice element `w`
w.r.t. the fundamental weights as a row vector.

!!! note
    The return type may not be relied on;
    we only guarantee that it is a one-dimensional iterable with `eltype` `ZZRingElem`
    that can be indexed with integers.
"""
function coefficients(w::WeightLatticeElem)
  return w.vec
end

@doc raw"""
    coeff(w::WeightLatticeElem, i::Int) -> ZZRingElem

Return the coefficient of the `i`-th fundamental weight in `w`.

This can be also accessed via `w[i]`.
"""
function coeff(w::WeightLatticeElem, i::Int)
  return w.vec[i]
end

function Base.getindex(w::WeightLatticeElem, i::Int)
  return coeff(w, i)
end

@doc raw"""
    conjugate_dominant_weight(w::WeightLatticeElem) -> WeightLatticeElem

Return the unique dominant weight conjugate to `w`.

See also: [`conjugate_dominant_weight_with_elem(::WeightLatticeElem)`](@ref).
"""
function conjugate_dominant_weight(w::WeightLatticeElem)
  return conjugate_dominant_weight!(deepcopy(w))
end

function conjugate_dominant_weight!(w::WeightLatticeElem)
  # w will be dominant once all fundamental weights have a positive coefficient,
  # so search for negative coefficients and make them positive by applying the corresponding reflection.
  s = 1
  while s <= rank(parent(w))
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
such that `w * x == dom`.
"""
function conjugate_dominant_weight_with_elem(w::WeightLatticeElem)
  return conjugate_dominant_weight_with_elem!(deepcopy(w))
end

function conjugate_dominant_weight_with_elem!(w::WeightLatticeElem)
  # determine the Weyl group element taking w to the fundamental chamber
  word = UInt8[]
  #sizehint!(word, count(<(0), coefficients(w))^2)
  s = 1
  while s <= rank(parent(w))
    if w[s] < 0
      push!(word, UInt8(s))
      reflect!(w, s)
      s = 1
    else
      s += 1
    end
  end

  # word is already in short lex normal form
  # and it is the element taking original w to new w
  return w, weyl_group(root_system(w))(word; normalize=false)
end

function dot(w1::WeightLatticeElem, w2::WeightLatticeElem)
  @req parent(w1) === parent(w2) "parent mismatch"
  R = root_system(w1)

  return dot(
    coefficients(w1),
    (coefficients(w2) * _cartan_symmetrizer_mat(R)) * cartan_matrix_inv(R),
  )
end

function expressify(w::WeightLatticeElem; context=nothing)
  if is_unicode_allowed()
    return expressify(w, :ω; context)
  else
    return expressify(w, :w; context)
  end
end

function expressify(w::WeightLatticeElem, s; context=nothing)
  sum = Expr(:call, :+)
  for i in 1:length(w.vec)
    push!(sum.args, Expr(:call, :*, expressify(w.vec[i]; context), "$(s)_$(i)"))
  end
  return sum
end
@enable_all_show_via_expressify WeightLatticeElem

@doc raw"""
    is_dominant(w::WeightLatticeElem) -> Bool

Check if `w` is a dominant weight, i.e. if all coefficients are non-negative.
"""
function is_dominant(w::WeightLatticeElem)
  return all(>=(0), coefficients(w))
end

@doc raw"""
    is_fundamental_weight(w::WeightLatticeElem) -> Bool

Check if `w` is a fundamental weight, i.e. exactly one coefficient is equal to 1 and all others are zero.

See also: [`is_fundamental_weight_with_index(::WeightLatticeElem)`](@ref).
"""
function is_fundamental_weight(w::WeightLatticeElem)
  fl, _ = is_fundamental_weight_with_index(w)
  return fl
end

function is_gen(w::WeightLatticeElem)
  return is_fundamental_weight(w)
end

@doc raw"""
    is_fundamental_weight_with_index(w::WeightLatticeElem) -> Bool, Int

Check if `w` is a fundamental weight and return this together with the index of the fundamental weight in [`fundamental_weights(::RootSystem)`](@ref fundamental_weights(root_system(w))).

If `w` is not a fundamental weight, the second return value is arbitrary.

See also: [`is_fundamental_weight(::WeightLatticeElem)`](@ref).
"""
function is_fundamental_weight_with_index(w::WeightLatticeElem)
  ind = 0
  coeffs = coefficients(w)
  for i in 1:size(coeffs, 2)
    if is_zero_entry(coeffs, 1, i)
      continue
    elseif is_one(coeffs[1, i])
      ind != 0 && return false, 0
      ind = i
    else
      return false, 0
    end
  end
  return ind != 0, ind
end

@doc raw"""
    reflect(w::WeightLatticeElem, s::Int) -> WeightLatticeElem

Return the reflection of `w` in the hyperplane orthogonal to the `s`-th simple root.

See also: [`reflect!(::WeightLatticeElem, ::Int)`](@ref).
"""
function reflect(w::WeightLatticeElem, s::Int)
  return reflect!(deepcopy(w), s)
end

@doc raw"""
    reflect!(w::WeightLatticeElem, s::Int) -> WeightLatticeElem
  
Reflect `w` in the hyperplane orthogonal to the `s`-th simple root, and return it.

This is a mutating version of [`reflect(::WeightLatticeElem, ::Int)`](@ref).
"""
function reflect!(w::WeightLatticeElem, s::Int)
  w.vec = addmul!(w.vec, view(cartan_matrix_tr(root_system(w)), s:s, :), -w.vec[s]) # change to submul! once available
  return w
end
