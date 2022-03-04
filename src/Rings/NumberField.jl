################################################################################
#
#  Number field as multivariate quotients K[x]/I
#
################################################################################

# We represent number fields as quotients of the form K[x]/I where I is a zero-
# dimensional irreducible ideal of a multivaraite polynomial ring over the
# rationals or any number field.
#
# Elements are just represented using elements of K[x].
#
# IMPORTANT
# We reduce the elements only
#  - before printing
#  - before multiplying
#  - before comparisons
#  - before computing the denominator (in the absolute case)

# We are tagging T to record the element type of the base field. This is needed
# in Hecke to distinguish between relative and absolute fields.
@attributes mutable struct NfNSGen{T, S} <: Hecke.NonSimpleNumField{T}
  I::MPolyIdeal{S}
  S::Vector{Symbol}
  degree::Int
  basis
  auxilliary_data::Vector{Any}

  function NfNSGen{U, V}(I, S::Vector{Symbol}) where {U, V}
    Kx = base_ring(I)::parent_type(V)
    K = base_ring(Kx)
    r = new{U, V}()
    # Determine the basis and degree
    r.I = I
    r.S = S
    B = basis(r)
    r.degree = length(B)
    return r
  end
end

mutable struct NfNSGenElem{T, S} <: Hecke.NumFieldElem{T}
  f::S
  parent::NfNSGen{T, S}
end

const NfAbsNSGen = NfNSGen{fmpq, fmpq_mpoly}

const NfAbsNSGenElem = NfNSGenElem{fmpq, fmpq_mpoly}

parent_type(::NfNSGenElem{T, S}) where {T, S} = NfNSGen{T, S}

parent_type(::Type{NfNSGenElem{T, S}}) where {T, S} = NfNSGen{T, S}

elem_type(::NfNSGen{T, S}) where {T, S} = NfNSGenElem{T, S}

elem_type(::Type{NfNSGen{T, S}}) where {T, S} = NfNSGenElem{T, S}

################################################################################
#
#  Constructors
#
################################################################################

function number_field(I::MPolyIdeal{fmpq_mpoly}, var::Vector{Symbol})
  n = length(var)
  R = base_ring(I)
  nn = nvars(R)
  @req length(var) == nvars(R) """
      Number of symbols $(n) must be the number of variables $(nn)
      """
  K = NfNSGen{fmpq, fmpq_mpoly}(I, var)
  return K, gens(K)
end

number_field(I::MPolyIdeal{fmpq_mpoly}, var::Vector{String}) =
    number_field(I, map(Symbol, var))

function number_field(I::MPolyIdeal{fmpq_mpoly}, var::String)
  R = base_ring(I)
  return number_field(I, String["$var$i" for i in 1:nvars(R)])
end

function number_field(I::MPolyIdeal{fmpq_mpoly}, var::Symbol)
  R = base_ring(I)
  return number_field(I, Symbol[Symbol(var, i) for i in 1:nvars(R)])
end

number_field(I::MPolyIdeal{fmpq_mpoly}) = number_field(I, "_\$")

# Now for number fields

function number_field(I::MPolyIdeal{Generic.MPoly{nf_elem}}, var::Vector{Symbol})
  n = length(var)
  R = base_ring(I)
  nn = nvars(R)
  @req length(var) == nvars(R) """
      Number of symbols $(n) must be the number of variables $(nn)
      """
  K = NfNSGen{nf_elem, Generic.MPoly{nf_elem}}(I, var)
  return K, gens(K)
end

number_field(I::MPolyIdeal{Generic.MPoly{nf_elem}}, var::Vector{String}) =
    number_field(I, map(Symbol, var))

function number_field(I::MPolyIdeal{Generic.MPoly{nf_elem}}, var::String)
  R = base_ring(I)
  return number_field(I, String["$var$i" for i in 1:nvars(R)])
end

function number_field(I::MPolyIdeal{Generic.MPoly{nf_elem}}, var::Symbol)
  R = base_ring(I)
  return number_field(I, Symbol[Symbol(var, i) for i in 1:nvars(R)])
end

number_field(I::MPolyIdeal{Generic.MPoly{nf_elem}}) = number_field(I, "_\$")

################################################################################
#
#  Basic field access
#
################################################################################

degree(K::NfNSGen) = K.degree

defining_ideal(K::NfNSGen) = K.I

polynomial_ring(K::NfNSGen{T, S}) where {T, S} = base_ring(defining_ideal(K))::parent_type(S)

gens(K::NfNSGen) = [K(x) for x = gens(polynomial_ring(K))]

gen(K::NfNSGen, i::Int) = K(gen(polynomial_ring(K), i))

ngens(K::NfNSGen) = ngens(polynomial_ring(K))

symbols(K::NfNSGen) = K.S

parent(a::NfNSGenElem) = a.parent

data(a::NfNSGenElem) = a.f

base_field(K::NfNSGen{fmpq, fmpq_mpoly}) = FlintQQ

base_field(K::NfNSGen) = base_ring(polynomial_ring(K))

################################################################################
#
#  String I/O
#
################################################################################

function show(io::IO, K::NfNSGen)
  Hecke.@show_name(io, K)
  Hecke.@show_special(io, K)
  io = IOContext(io, :compact => true)
  print(io, "Number field defined by\n")
  G = gens(defining_ideal(K))
  print(io, "[")
  for i in 1:length(G)
    print(io, G[i])
    if i < length(G)
      print(io, ", ")
    end
  end
  print(io, "]")
end

function AbstractAlgebra.expressify(a::NfNSGenElem; context = nothing)
  return expressify(a.f, symbols(parent(a)), context = context)
end

function show(io::IO, a::NfNSGenElem)
  reduce!(a)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

################################################################################
#
#  Reduction
#
################################################################################

function assert_has_gb(K::NfNSGen)
  I = defining_ideal(K)
  if isdefined(I, :gb) && isdefined(I.gb, :S) && I.gb.S.isGB
    return nothing
  end
  I = defining_ideal(K)
  groebner_assure(I, degrevlex(gens(base_ring(I))))
  GI = first(values(I.gb))
  singular_assure(GI)
  GI.S.isGB = true
  return nothing
end

function reduce!(a::NfNSGenElem)
  K = parent(a)
  assert_has_gb(K)
  I = defining_ideal(K)
  GI = collect(values(I.gb))[1]
  Sx = base_ring(GI.S)
  f = a.f
  a.f = I.gens.Ox(reduce(Sx(f), GI.S))
  return a
end

function check_parent(a::NfNSGenElem{S, T}, b::NfNSGenElem{S, T}) where {S, T}
  if parent(a) !== parent(b)
    throw(ArgumentError("Parents of elements must be equal"))
  end
end

################################################################################
#
#  Ring interface functions
#
################################################################################

canonical_unit(a::NfNSGenElem) = a

zero(K::NfNSGen) = K(0)

one(K::NfNSGen) = K(1)

# Our ideals are non-trivial, so we don't have to reduce
isone(a::NfNSGenElem) = isone(data(reduce!(a)))

iszero(a::NfNSGenElem) = iszero(data(reduce!(a)))

###############################################################################3
#
#  Binary operations
#
################################################################################

function +(a::NfNSGenElem, b::NfNSGenElem)
  check_parent(a, b)
  return NfNSGenElem(data(a) + data(b), parent(a))
end

function -(a::NfNSGenElem, b::NfNSGenElem)
  check_parent(a, b)
  NfNSGenElem(data(a) - data(b), parent(a))
end

function divexact(a::NfNSGenElem, b::NfNSGenElem; check::Bool = true)
  @req !iszero(b) "Division by zero"
  check_parent(a, b)
  return a * inv(b)
end

Base.://(a::NfNSGenElem, b::NfNSGenElem) = divexact(a, b)

function Base.:(*)(a::NfNSGenElem, b::NfNSGenElem)
  check_parent(a, b)
  reduce!(a)
  reduce!(b)
  return NfNSGenElem(data(a) * data(b), parent(a))
end

################################################################################
#
#  Unary operations
#
################################################################################

-(a::NfNSGenElem) = NfNSGenElem(-data(a), parent(a))

################################################################################
#
#  Power by squaring
#
################################################################################

function Base.:(^)(a::NfNSGenElem, i::Int)
  if iszero(i) # x^0 = 1 for all x
    return one(parent(a))
  elseif iszero(a)
    return zero(parent(a))
  elseif isone(i)
    return deepcopy(a)
  else
    bit = ~((~UInt(0)) >> 1)
    while (UInt(bit) & i) == 0
      bit >>= 1
    end
    z = deepcopy(a)
    bit >>= 1
    while bit != 0
      mul!(z, z, z)
      if (UInt(bit) & i) != 0
        mul!(z, z, a)
      end
      bit >>= 1
    end
    return z
  end
end

################################################################################
#
#  Inverse
#
################################################################################

function inv(a::NfNSGenElem)
  @req !iszero(a) "Element must be non-zero"
  f = minpoly(a)
  z = parent(a)(coeff(f, degree(f)))
  for i=degree(f)-1:-1:1
    z = z*a + coeff(f, i)
  end
  return -z*inv(coeff(f, 0))
end

################################################################################
#
#  Ad hoc binary operations
#
################################################################################

# with base type

+(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::nf_elem) =
    NfNSGenElem(data(a) + b, parent(a))

*(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::nf_elem) =
    NfNSGenElem(data(a) * b, parent(a))

Base.://(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::nf_elem) =
    NfNSGenElem(divexact(data(a), b), parent(a))

divexact(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::nf_elem) = a//b

+(a::nf_elem, b::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}) = b + a

*(a::nf_elem, b::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}) = b * a

Base.://(a::nf_elem, b::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}) = a * inv(b)

divexact(a::nf_elem, b::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}) = a//b

# with fmpq
+(a::NfNSGenElem{T, S}, b::fmpq) where {S, T} =
    NfNSGenElem(data(a) + b, parent(a))

*(a::NfNSGenElem{T, S}, b::fmpq) where {S, T} =
    NfNSGenElem(data(a) * b, parent(a))

Base.://(a::NfNSGenElem{T, S}, b::fmpq) where {S, T} =
    NfNSGenElem(divexact(data(a), b), parent(a))

divexact(a::NfNSGenElem{T, S}, b::fmpq) where {S, T} = a//b

+(a::fmpq, b::NfNSGenElem{T, S}) where {S, T} = b + a

*(a::fmpq, b::NfNSGenElem{T, S}) where {S, T} = b * a

Base.://(a::fmpq, b::NfNSGenElem{T, S}) where {S, T} = a * inv(b)

divexact(a::fmpq, b::NfNSGenElem{T, S}) where {S, T} = a//b

# with fmpz
+(a::NfNSGenElem{T, S}, b::fmpz) where {S, T} =
    NfNSGenElem(data(a) + b, parent(a))

*(a::NfNSGenElem{T, S}, b::fmpz) where {S, T} =
    NfNSGenElem(data(a) * b, parent(a))

Base.://(a::NfNSGenElem{T, S}, b::fmpz) where {S, T} =
    NfNSGenElem(divexact(data(a), b), parent(a))

divexact(a::NfNSGenElem{T, S}, b::fmpz) where {S, T} = a//b

+(a::fmpz, b::NfNSGenElem{T, S}) where {S, T} = b + a

*(a::fmpz, b::NfNSGenElem{T, S}) where {S, T} = b * a

Base.://(a::fmpz, b::NfNSGenElem{T, S}) where {S, T} = a * inv(b)

divexact(a::fmpz, b::NfNSGenElem{T, S}) where {S, T} = a//b

# with Integer
+(a::NfNSGenElem{T, S}, b::Base.Integer) where {S, T} =
    NfNSGenElem(data(a) + b, parent(a))

*(a::NfNSGenElem{T, S}, b::Base.Integer) where {S, T} =
    NfNSGenElem(data(a) * b, parent(a))

Base.://(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::Base.Integer) =
    NfNSGenElem(divexact(data(a), b), parent(a))

divexact(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::Base.Integer) = a//b

+(a::Base.Integer, b::NfNSGenElem{T, S}) where {S, T} = b + a

*(a::Base.Integer, b::NfNSGenElem{T, S}) where {S, T} = b * a

Base.://(a::Base.Integer, b::NfNSGenElem{T, S}) where {S, T} = a * inv(b)

divexact(a::Base.Integer, b::NfNSGenElem{T, S}) where {S, T} = a//b

# with Rational
+(a::NfNSGenElem{T, S}, b::Base.Rational{<:Base.Integer}) where {T, S} =
    NfNSGenElem(data(a) + b, parent(a))

*(a::NfNSGenElem{T, S}, b::Base.Rational{<:Base.Integer}) where {S, T} =
    NfNSGenElem(data(a) * b, parent(a))

Base.://(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::Base.Rational{<:Base.Integer}) =
    NfNSGenElem(divexact(data(a), b), parent(a))

divexact(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::Base.Rational{<:Base.Integer}) =
    a//b

+(a::Base.Rational{<:Base.Integer}, b::NfNSGenElem{T, S}) where {S, T} = b + a

*(a::Base.Rational{<:Base.Integer}, b::NfNSGenElem{T, S}) where {S, T} = b * a

Base.://(a::Base.Rational{<:Base.Integer}, b::NfNSGenElem{T, S}) where {S, T} = a * inv(b)

divexact(a::Base.Rational{<:Base.Integer}, b::NfNSGenElem{T, S}) where {S, T} =
    a//b

################################################################################
#
#  Deepcopy
#
################################################################################

function Base.deepcopy_internal(a::NfNSGenElem, dict::IdDict)
  return NfNSGenElem(Base.deepcopy_internal(a.f, dict), parent(a))
end

################################################################################
#
#  In place operations
#
################################################################################

function Oscar.mul!(a::NfNSGenElem, b::NfNSGenElem, c::NfNSGenElem)
  a.f = b.f*c.f
  return a
end

function Oscar.add!(a::NfNSGenElem, b::NfNSGenElem, c::NfNSGenElem)
  a.f = b.f + c.f
  return a
end

function Oscar.addeq!(a::NfNSGenElem, b::NfNSGenElem)
  a.f += b.f
  return a
end

################################################################################
#
#  Comparison
#
################################################################################

function ==(a::NfNSGenElem{T, S}, b::NfNSGenElem{T, S}) where {T, S}
  reduce!(a)
  reduce!(b)
  return a.f == b.f
end

################################################################################
#
#  Ad hoc comparisons
#
################################################################################

# Our ideals are non-trivial, so we don't have to reduce before comparing
# with scalars
for t in [Base.Integer, Base.Rational{<:Integer}, fmpz, fmpq]
  @eval begin
    ==(a::NfNSGenElem{T, S}, b::$t) where {T, S} = data(a) == b

    ==(a::$t, b::NfNSGenElem{T, S}) where {T, S} = data(b) == a
  end
end

==(a::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}, b::nf_elem) = data(a) == b

==(a::nf_elem, b::NfNSGenElem{nf_elem, Generic.MPoly{nf_elem}}) where {T, S} = a == data(b)

==(a::NfNSGenElem{fmpq, fmpq_mpoly}, b::fmpq) = data(a) == b

==(a::fmpq, b::NfNSGenElem{fmpq, fmpq_mpoly}) where {T, S} = data(b) == a

################################################################################
#
#  Parent call overloading
#
################################################################################

(K::NfNSGen)() = zero(K)

function (K::NfNSGen{T, S})(a::NfNSGenElem{T, S}) where {T, S}
  parent(a) !== K && error("Parents do not match")
  return a
end

function (K::NfNSGen{T, S})(a::S) where {T, S}
  parent(a) !== polynomial_ring(K) &&
      error("Parents do not match")

  return NfNSGenElem(a, K)
end

function (K::NfNSGen{T, S})(a::T) where {T, S}
  parent(a) !== base_field(K) && 
      error("Parents do not match")

  return K(polynomial_ring(K)(a))
end

for t in [Base.Integer, Base.Rational{<:Base.Integer}, fmpz, fmpq]
  @eval begin
    function (K::NfNSGen{T, S})(a::$t) where {T, S}
      return K(polynomial_ring(K)(a))
    end
  end
end

################################################################################
#
#  Denominator (for absolute fields)
#
################################################################################

Hecke.denominator(a::NfNSGenElem{fmpq, fmpq_mpoly}) =
    denominator(data(reduce!(a)))

################################################################################
#
#  Basis
#
################################################################################

function basis(K::NfNSGen; copy::Bool = true)
  if isdefined(K, :basis)
    B = K.basis::Vector{elem_type(K)}
    return copy ? Base.deepcopy(B) : B
  else
    I = defining_ideal(K)
    assert_has_gb(K)
    GI = first(values(I.gb))
    s = Singular.kbase(GI.S)
    if iszero(s)
      error("ideal was not zero-dimensional")
    end
    B = elem_type(K)[K(base_ring(defining_ideal(K))(x)) for x = gens(s)]
    if !isone(B[1])
      i = findfirst(isone, B)
      B[1], B[i] = B[i], B[1]
    end
    K.basis = B
    return copy ? Base.deepcopy(B) : B
  end
end

################################################################################
#
#  Coordinates
#
################################################################################

function coordinates(a::NfNSGenElem{T, S}) where {T, S}
  K = parent(a)
  B = basis(K, copy = false)
  v = Vector{T}(undef, degree(K))
  reduce!(a)
  g = a.f
  for j in 1:length(B)
    v[j] = coeff(g, data(B[j]))
  end
  @assert dot(v, B) == a
  return v
end

absolute_coordinates(a::NfNSGenElem{fmpq, fmpq_mpoly}) = coordinates(a)

function absolute_coordinates(a::NfNSGenElem)
  K = parent(a)
  v = Vector{fmpq}(undef, absolute_degree(K))
  va = coordinates(a)
  ind = 1
  for i = 1:length(va)
    vi = absolute_coordinates(va[i])
    for j = 1:length(vi)
      v[ind] = vi[j]
      ind += 1
    end
  end
  return v
end

function elem_to_mat_row!(M::MatElem{T}, i::Int, a::NfNSGenElem{T, S}) where {T, S}
  v = coordinates(a)
  for j in 1:length(v)
    M[i, j] = v[j]
  end
  return nothing
end

function elem_to_mat_row!(M::fmpz_mat, i::Int, d::fmpz, a::NfNSGenElem{fmpq, fmpq_mpoly})
  dd = denominator(a)
  v = coordinates(dd * a)
  for j in 1:length(v)
    @assert isone(denominator(v[j]))
    M[i, j] = numerator(v[j])
  end
  ccall((:fmpz_set, Nemo.libflint), Nothing, (Ref{fmpz}, Ref{fmpz}), d, dd)
  return nothing
end

function basis_matrix(v::Vector{NfNSGenElem{fmpq, fmpq_mpoly}},
                      ::Type{FakeFmpqMat})
  d = degree(parent(v[1]))
  z = zero_matrix(FlintQQ, length(v), d)
  for i in 1:length(v)
    elem_to_mat_row!(z, i, v[i])
  end
  return FakeFmpqMat(z)
end

################################################################################
#
#  Representation matrix
#
################################################################################

function representation_matrix(a::NfNSGenElem)
  K = parent(a)
  b = basis(K, copy = false)
  M = zero_matrix(base_field(K), degree(K), degree(K))
  for i=1:degree(K)
    elem_to_mat_row!(M, i, a*b[i])
  end
  return M
end

function representation_matrix_q(a::NfNSGenElem{fmpq, fmpq_mpoly})
  M = representation_matrix(a)
  return Hecke._fmpq_mat_to_fmpz_mat_den(M)
end

function Hecke.elem_from_mat_row(K::NfNSGen{fmpq, fmpq_mpoly},
                                 M::fmpz_mat, i::Int, d::fmpz)
  z = zero(K)
  B = basis(K, copy = false)
  for j in 1:ncols(M)
    z = z + M[i, j] * B[j]
  end
  return K(fmpq(1, d)) * z
end

################################################################################
#
#  Minimal polynomial
#
################################################################################

function minpoly(a::NfNSGenElem)
  K = parent(a)
  n = degree(K)
  k = base_field(K)
  M = zero_matrix(k, degree(K)+1, degree(K))
  z = a^0
  elem_to_mat_row!(M, 1, z)
  z *= a
  elem_to_mat_row!(M, 2, z)
  i = 2
  Qt, _ = PolynomialRing(k, "t", cached = false)
  while true
    if n % (i-1) == 0 && rank(M) < i
      N = nullspace(transpose(sub(M, 1:i, 1:ncols(M))))
      @assert N[1] == 1
      v = Vector{elem_type(k)}(undef, i)
      for j in 1:i
        v[j] = N[2][j, 1]
      end
      f = Qt(v)
      return f*inv(leading_coefficient(f))
    end
    z *= a
    elem_to_mat_row!(M, i+1, z)
    i += 1
  end
end

################################################################################
#
#  Characteristic polynomial
#
################################################################################

function Hecke.charpoly(a::NfNSGenElem)
  f = minpoly(a)
  return f^div(degree(parent(a)), degree(f))
end

################################################################################
#
#  Trace
#
################################################################################

function Hecke.tr(a::NfNSGenElem)
  return -coeff(charpoly(a), degree(parent(a)) - 1)
end

################################################################################
#
#  Norm
#
################################################################################

function Hecke.norm(a::NfNSGenElem)
  f = charpoly(a)
  n = degree(f)
  return isodd(n) ? -coeff(f, 0) : coeff(f, 0)
end


################################################################################
#
#  Computation of maximal order
#
################################################################################

function Hecke.any_order(K::NfNSGen)
  B = basis(K, copy = false)
  for b in B
    @assert isone(denominator(b.f))
  end

  O = Order(K, B)

  return O
end

function Hecke.MaximalOrder(K::NfNSGen)
  E = any_order(K)
  return Hecke._maximal_order_round_four(E)
end

Hecke.isdefining_polynomial_nice(::NfNSGen) = false

################################################################################
#
#  Random
#
################################################################################

RandomExtensions.maketype(K::NfNSGen, r) = elem_type(K)

function rand(rng::AbstractRNG, sp::RandomExtensions.SamplerTrivial{<:RandomExtensions.Make2{<:NfNSGenElem,<:NfNSGen,<:UnitRange}})
  K, r = sp[][1:end]
  # TODO: This is super slow
  b = basis(K, copy = false)
  z::Random.gentype(sp) = K() # type-assert to help inference on Julia 1.0 and 1.1
  for i in 1:degree(K)
    z += rand(rng, r) * b[i]
  end
  return z
end

rand(K::NfNSGen, r::UnitRange) = rand(Random.GLOBAL_RNG, K, r)
rand(rng::AbstractRNG, K::NfNSGen, r::UnitRange) = rand(rng, make(K, r))

################################################################################
#
#  Primitive element
#
################################################################################

function primitive_element(K::NfNSGen)
  r = Random.MersenneTwister(1)
  d = degree(K)
  b = 1
  cnt = 0
  while true
    cnt += 1
    z = rand(r, K, -b:b)
    if degree(minpoly(z)) == d
      return z
    end
    if cnt >= 10
      b = 2 * b
    end
  end
end

################################################################################
#
#  Simple extension
#
################################################################################

function Hecke.simple_extension(K::NfNSGen; simplify = false)
  a = primitive_element(K)
  if simplify
    L = NumberField(minpoly(a), "\$", cached = false)[1]
    h = hom(L, K, a)
    LL, mLL = Hecke.simplify(L, cached = false)
    return LL, mLL * h
  else
    L = NumberField(minpoly(a), "\$", cached = false)[1]
    return L, hom(L, K, a)
  end
end

################################################################################
#
#  Hook into Hecke number field maps
#
################################################################################

# From NfAbsNSGen into something
mutable struct MapDataFromNfAbsNSGen{T}
  images::T
  isid::Bool

  function MapDataFromNfAbsNSGen{T}(x::T) where {T}
    z = new{T}(x, false)
    return z
  end
  
  function MapDataFromNfAbsNSGen{T}(x::Bool) where {T}
    @assert x
    z = new{T}()
    z.isid = true
    return z
  end
end

function Hecke._isequal(K, L, u::MapDataFromNfAbsNSGen{T}, v::MapDataFromNfAbsNSGen{T}) where {T}
  # If one is the identity, this means that K === L
  if (u.isid && !v.isid)
    return v.images == gens(K)
  elseif (v.isid && !u.isid)
    return u.images == gens(K)
  elseif u.isid && v.isid
    return true
  end

  return v.images == u.images 
end

function image(f::MapDataFromNfAbsNSGen, L, y)
  f.isid && return L(y)
  return Hecke.msubst(y.f, f.images)
end

Hecke.map_data_type(K::NfAbsNSGen, L) = MapDataFromNfAbsNSGen{Vector{elem_type(L)}}

Hecke.map_data_type(T::Type{NfAbsNSGen}, L::Type) = MapDataFromNfAbsNSGen{Vector{elem_type(L)}}

Hecke.map_data(K::NfAbsNSGen, L, ::Bool) = MapDataFromNfAbsNSGen{Vector{elem_type(L)}}(true)

function Hecke.map_data(K::NfAbsNSGen, L, x::Vector; check = true)
  if length(x) != ngens(K)
    error("Data does not define a morphism")
  end

  local xx::Vector{elem_type(L)}

  if x isa Vector{elem_type(L)}
    if parent(x[1]) !== L
      error("Data does not define a morphism")
    end
    xx = x
  else
    xx = map(L, x)::Vector{elem_type(L)}
  end

  if check
    for f in gens(defining_ideal(K))
      if !iszero(evaluate(f, xx))
        error("")
      end
    end
  end 

  @assert typeof(xx) == Vector{elem_type(L)}

  return MapDataFromNfAbsNSGen{typeof(xx)}(xx)
end 

function Hecke.image_generators(f::Hecke.NumFieldMor{<:NfAbsNSGen})
  return f.image_data.images
end

# into NfAbsNSGen
function Hecke._compute_inverse_data(f#= image data =#, K, L::NfAbsNSGen)
  return Hecke._compute_inverse_data(f, K, L, L)
end

function Hecke._compute_inverse_data(f#= image data =#, K, LL, L::NfAbsNSGen)
  preimg_gens = elem_type(K)[]
  for g in gens(L)
    fl, preimg = haspreimage(f, LL(g))
    @assert fl
    push!(preimg_gens, preimg)
  end
  return MapDataFromNfAbsNSGen{typeof(preimg_gens)}(preimg_gens)
end

# compose
function Hecke._compose(f::MapDataFromNfAbsNSGen, g#= map data =#, K, L, M)
  return Hecke.map_data_type(K, M)(elem_type(M)[image(g, M, image(f, L, gg)) for gg in gens(K)])
end

# From NfNSGen into something
mutable struct MapDataFromNfNSGen{T, S}
  images::T
  base_field_map_data::S
  isid::Bool

  function MapDataFromNfNSGen{T, S}(x::T, y::S) where {T, S}
    z = new{T, S}(x, y, false)
    return z
  end
  
  function MapDataFromNfNSGen{T, S}(x::Bool) where {T, S}
    @assert x
    z = new{T, S}()
    z.isid = true
    return z
  end
end

function Hecke._isequal(K, L, u::MapDataFromNfNSGen{T, S}, v::MapDataFromNfNSGen{T, S}) where {T, S}
   if u.isid && v.isid
    return true
  end

  return all(g -> image(u, L, g) == image(v, L, g), gens(K)) && Hecke._isequal(base_field(K), base_field(L), u.base_field_map_data, v.base_field_map_data)
end

function image(f::MapDataFromNfNSGen, L, y)
  f.isid && return L(y)
  z = map_coefficients(w -> image(f.base_field_map_data, L, w), data(y), cached = false)
  return evaluate(z, f.images)
end

function Hecke.map_data_type(T::Type{<: NfNSGen}, L::Type{<:Any})
  MapDataFromNfNSGen{Vector{elem_type(L)}, map_data_type(base_field_type(T), L)}
end

Hecke.map_data_type(K::NfNSGen, L) = MapDataFromNfNSGen{Vector{elem_type(L)}, Hecke.map_data_type(base_field(K), L)}

function Hecke.map_data(K::NfNSGen, L, ::Bool)
  z = MapDataFromNfNSGen{Vector{elem_type(L)}, Hecke.map_data_type(base_field(K), L)}(true)
  z.base_field_map_data = Hecke.map_data(base_field(K), L, true)
  return z
end

function Hecke.map_data(K::NfNSGen, L, x...; check = true)
  z = Hecke.map_data(base_field(K), L, Base.front(x)...; check = check)

  local yy::Vector{elem_type(L)}

  if Base.last(x) isa Hecke.NumFieldMor
    domain(Base.last(x)) !== K && error("")
    _y = image_generators(Base.last(x))
    if parent(_y[1]) === L
      yy = _y
    else
      yy = map(L, _y)::Vector{elem_type(L)}
    end
  else
    y = Base.last(x)
    if !(y isa Vector)
      error("")
    end
    if parent(y[1]) === L
      yy = y
    else
      yy = map(L, y)::Vector{elem_type(L)}
    end
  end

  if check
    for i in 1:ngens(K)
      w = evaluate(map_coefficients(w -> image(z, L, w), gens(defining_ideal(K))[i], cached = false), yy)
      !iszero(w) && error("Data does not define a morphism")
    end
  end
  
  @assert typeof(yy) == Vector{elem_type(L)}
  @assert typeof(z) == Hecke.map_data_type(base_field(K), L)

  return MapDataFromNfNSGen{typeof(yy), typeof(z)}(yy, z)
end

# into NfNSGen
function Hecke._compute_inverse_data(f#= image data =#, K, L::NfNSGen)
  return Hecke._compute_inverse_data(f, K, L, L)
end

function Hecke._compute_inverse_data(f, K, LL, L::NfNSGen)
  preimg_gens = elem_type(K)[]
  for g in gens(L)
    fl, preimg = haspreimage(f, LL(g))
    push!(preimg_gens, preimg)
  end
  inverse_data_base_field = Hecke._compute_inverse_data(f, K, LL, base_field(L))
  return MapDataFromNfRelNS{typeof(preimg_gens), typeof(inverse_data_base_field)}(preimg_gens, inverse_data_base_field)
end

# compose
function Hecke._compose(f::MapDataFromNfNSGen, g#= map data =#, K, L, M)
  return Hecke.map_data_type(K, M)(elem_type(M)[image(g, M, image(f, L, u)) for u in gens(K)],
                             Hecke._compose(f.base_field_map_data, g, base_field(K), L, M))
end

# Image of genertors

function image_generators(f::Hecke.NumFieldMor{<:NfNSGen})
  return f.image_data.images
end
