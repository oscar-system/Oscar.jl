module StandardFiniteFields

using Oscar
export standard_finite_field


const IntegerUnion = Oscar.IntegerUnion
const PrimeField = Union{fpField,FpField}
const PrimeFieldElem = Union{fpFieldElem,FpFieldElem}
const PrimeFieldMatrix = Union{FpMatrix,fpMatrix}

# NOTE: These give missing features to OSCAR/Nemo that will likely be added in the near future.

function pop_largest_factor!(f::Fac{ZZRingElem})
  D = f.fac
  m = maximum(f)
  if isone(m[2])
    Base.delete!(D, m[1])
  else
    D[m[1]] -= 1
  end
  return m
end

# TODO : Should be fixed in Nemo
function (k::FpField)(a::Vector)
  @assert isone(length(a))
  return k(a[1])
end

function (k::FqPolyRepField)(a::Vector)
  return k(polynomial(Native.GF(ZZ(characteristic(k))), a))
end

# TODO : Should be fixed in Hecke for prime fields
function coords(x::FinFieldElem)
  return absolute_coordinates(x)
end
function coords(x::PrimeFieldElem)
  return [x]
end

function largest_factor(n::IntegerUnion)
  nfactorization = factor(ZZ(n))
  return maximum(nfactorization)
end

# _attributes = [
# :is_standard_finite_field        - bool
# :is_standard_prime_field         - bool
# :primitive_powers_in_tower_basis - Matrix{FinFieldElem}
# :tower_basis                     - Matrix{FFieldElem}
#                                       = inv(:primitive_powers_in_tower_basis)
# :steinitz_prime_degree           - Dict{Int, Dict{Int, ZZRingElem}}
# :standard_extensions             - Dict{Int, FinField}
# ]
function set_standard_prime_field!(F::PrimeField)
  get_attribute!(F, :is_standard_prime_field) do
    set_attribute!(F, :standard_extensions, Dict{ZZRingElem,FinField}(1 => F))
    set_attribute!(F, :primitive_powers_in_tower_basis, identity_matrix(F, 1))
    true
  end
end
function set_standard_finite_field!(F::FinField)
  set_attribute!(F, :is_standard_finite_field, true)
end
function set_primitive_powers_in_tower_basis!(F::FinField, M::PrimeFieldMatrix)
  set_attribute!(F, :primitive_powers_in_tower_basis, M)
  set_attribute!(F, :tower_basis, inv(M))
end
function set_standard_extension!(F::PrimeField, n::IntegerUnion, K::FinField)
  ext = get_attribute!(F, :standard_extensions, Dict{Int,FinField}())
  ext[n] = K
end


function is_standard_prime_field(F::PrimeField)
  get_attribute(F, :is_standard_prime_field, false)
end
function is_standard_finite_field(F::FinField)
  get_attribute(F, :is_standard_finite_field, false)
end
function primitive_powers_in_tower_basis(F::FinField)
  get_attribute(F, :primitive_powers_in_tower_basis, nothing)
end
function tower_basis(F::FinField)
  get_attribute(F, :tower_basis, nothing)
end
function get_steinitz_prime_degree(F::FinField, r::IntegerUnion, k::IntegerUnion)
  spd = get_attribute(F, :steinitz_prime_degree, nothing)
  spdr = get(spd, r, nothing)
  get(spdr, k, nothing)
end
function get_steinitz_prime_degree!(
  f::Function,
  F::FinField,
  r::IntegerUnion,
  k::IntegerUnion,
)
  spd = get_attribute!(F, :steinitz_prime_degree, Dict{Int,Dict{Int,ZZRingElem}}())
  spdr = get!(spd, r, Dict{Int,ZZRingElem}())
  get!(spdr, k, f())
end
function get_standard_extensions(F::PrimeField)
  get_attribute(F, :standard_extensions, nothing)
end
function get_standard_extensions!(F::PrimeField)
  get_attribute!(F, :standard_extensions) do
    Dict{Int,FinField}()
  end
end
function get_standard_extension(F::PrimeField, k::IntegerUnion)
  ext = get_standard_extensions(F)
  get(ext, k, nothing)
end
function get_standard_extension!(F::PrimeField, k::IntegerUnion, L::FinField)
  ext = get_standard_extensions!(F)
  get!(ext, k, L)
end
function get_standard_extension!(f::Function, F::PrimeField, k::IntegerUnion)
  ext = get_standard_extensions!(F)
  get!(f, ext, k)
end


function standard_affine_shift_data(q::IntegerUnion)
  m = div(4 * q, 5)
  while gcd(m, q) != 1
    m -= 1
  end
  a = div(2 * q, 3)
  return m, a
end

function standard_affine_shift(q::IntegerUnion, i::IntegerUnion)
  m, a = standard_affine_shift_data(q)
  return mod((m * i + a), q)
end

# Given a standard finite field F and Steinitz number n,
# give the corresponding field element.
function element_from_steinitz_number(F::PrimeField, n::IntegerUnion)
  @req 0 <= n <= order(F) "We need to have 0 <= n <= q"
  return F(n)
end
function element_from_steinitz_number(F::FinField, n::IntegerUnion)
  @req is_standard_finite_field(F) "First input must be a standard finite field"
  p = characteristic(F)
  @req 0 <= n <= order(F) "We need to have 0 <= n <= q"
  iszero(n) && return zero(F)

  # this forms a linear combo of F.towervasis rows using vectorrep as coefficients,
  # and then convert this vector to an element of F.
  vectorrep = digits(n, base = Int(p))
  return F(vectorrep * @view tower_basis(F)[1:length(vectorrep), :])
end

# Returns an element a in F that is NOT an rth root of unity
# we REQUIRE that F is a standard finite field
function non_rth_root(F::FinField, r::IntegerUnion)
  @assert is_standard_finite_field(F) || is_standard_prime_field(F)
  q = order(F)
  if iszero(mod(q - 1, r))
    i = 0
    a = zero(F)
    k = divexact(q - 1, r)
    while iszero(a) || isone(a^k)
      i += 1
      a = element_from_steinitz_number(F, standard_affine_shift(q, i))
    end
    return a
  else
    return nothing
  end
end

function standard_irreducible_coefficient_list(
  F::FinField,
  r::IntegerUnion,
  a::FinFieldElem,
)
  q = order(F)
  l = zeros(F, Int(r) + 1)
  l[Int(r)+1] = one(F)
  l[1] = a
  l[2] = one(F)
  # inc is the expected number of nonzero coefficients
  inc = 1
  let t = q
    while t < 2 * r
      t *= q
      inc += 1
    end
  end
  # allowing non-zero coeffs up to position d
  # after every r attempts allow inc more non-zero coeffs
  d = 0
  count = 0
  qq = 0
  while !is_irreducible(polynomial(F, l))
    if iszero(mod(count, r)) && d < r - 1
      d += inc
      if d >= r
        d = r - 1
      end
      qq = q^(d - 1)
    end
    # q can be very very large so Int is not big enough...
    st = digits(standard_affine_shift(qq, count), base = BigInt(q), pad = d - 1)
    # TODO: we can remove this while loop when fix for padding digits live
    while length(st) < d - 1
      push!(st, 0)
    end
    for k in 2:d
      l[k] = element_from_steinitz_number(F, st[k-1])
    end
    count += 1
  end
  return l
end

# returns the Steinitz number corresponding to the polynomial g(X),
# where f = X^r + g(X) is the standard irreducible polynomial over FF(p, r^(k-1))
function steinitz_number_for_prime_degree(p::IntegerUnion, r::IntegerUnion, k::IntegerUnion)
  Fp = standard_finite_field(p, 1)

  get_steinitz_prime_degree!(Fp, r, k) do
    # now we need to create the polynomial depending on the prime r
    if r == p
      # Artin-Schreier case
      # k = 1 we get [(Xr[1])^p - (Xr[1]) -1]
      # k > 1 we get (Xr[k])^p - (Xr[k]) - (prod(Xr[j] : j in [1..k-1]))^(p-1))
      q = ZZ(p)^(p^(k - 1))
      return (p - 1) * (q + divexact(q, p))
    elseif r == 2 && mod(p, 4) == 3
      if k == 1
        # (Xr[1])^2 +1
        return 1
      elseif k == 2
        a = non_rth_root(standard_finite_field(p, 2), r)
        # Xr[2]^2 -a
        return steinitz_number(-a)
      else
        # Xr[i]^2 - Xr[i-1]
        return (p - 1) * ZZ(p)^(r^(k - 2))
      end
    elseif r == 2
      if k == 1
        # Xr[1]^2 -a
        a = non_rth_root(standard_finite_field(p, 1), r)
        return steinitz_number(-a)
      else
        # Xr[j]^r - Xr[j-1]
        return (p - 1) * ZZ(p)^(r^(k - 2))
      end
    else
      # Here we use pseudo-random polynomials...
      F = standard_finite_field(p, r^(k - 1))
      if k == 1
        a = -one(F)
      else
        a = -gen(F)
      end
      l = standard_irreducible_coefficient_list(F, r, a)
      pop!(l)
      while is_zero(l[end])
        pop!(l)
      end
      return evaluate(polynomial(ZZ, steinitz_number.(l)), order(F))
    end
  end
end

# x will be represented internally as a polynomial over Fp in the generator of F.
# We need to first convert this to an Fp-vector, then to the corresponding vector
# with respect to the Tower Basis.
# Then we think of this vector as a polynomial (over ZZ) in a temporary indeterminate z,
# and evaluate at z = char(F) to get the Steinitz number.
# NOTE for whatever reason, evaluate(polynomial(), ) is faster than evalpoly()
function steinitz_number(F::PrimeField, x::PrimeFieldElem)
  @assert parent(x) === F
  return ZZ(lift(x))
end
function steinitz_number(F::FinField, x::FinFieldElem)
  @assert parent(x) === F
  v = lift.(absolute_coordinates(x) * primitive_powers_in_tower_basis(F))
  return evaluate(polynomial(ZZ, v), characteristic(F))
end
function steinitz_number(x::FinFieldElem)
  return steinitz_number(parent(x), x)
end

# describes monomials in tower basis plus their degrees
function standard_monomial(n::IntegerUnion)
  error("not implemented")
end
# just returns degrees of monomials in tower basis
# TODO : pass in factorization? Do we need this with memoization?
function standard_monomial_degrees(n::IntegerUnion)
  if n == 1
    return [1]
  end
  # need the largest prime factor a of n
  a, k = largest_factor(n)
  res = standard_monomial_degrees(divexact(n, a))
  m = a^k
  new = map(x -> lcm(x, m), res)
  for i in 1:a-1
    append!(res, new)
  end
  return res::Vector{Int}
end
# map of monomials for degree n -> monomials of degree m by positions
function standard_monomial_map(n::IntegerUnion, m::IntegerUnion)
  d = standard_monomial_degrees(m)
  return [i for i in 1:length(d) if mod(n, d[i]) == 0]
end

# Embed an element x of Fp^n into Fp^m by Steinitz numbers
# where nr = steinitz_number(Fp^n, x)
# I hate hate hate these variable names copied (mostly) from Lübeck
function embed_steinitz(p::IntegerUnion, n::IntegerUnion, m::IntegerUnion, nr::IntegerUnion)
  if n == m || iszero(nr)
    return nr
  end
  l = digits(nr, base = Int(p))
  M = @view standard_monomial_map(n, m)[1:length(l)]
  c = zeros(ZZRingElem, M[end])
  c[M] = l
  return evaluate(polynomial(ZZ, c), p)
end


# Given a field K, we construct an extension L with [L:K] = deg
# We use the irreducible polynomial f = X^deg  + g(X)
#    where lcoeffs are the coefficients of g(X).
# We assume b is a generator for K, and so bX will be a generator for L
function _extension_with_tower_basis(
  K::PrimeField,
  deg::IntegerUnion,
  lcoeffs::Vector,
  b::PrimeFieldElem,
)
  @assert parent(b) === K

  while length(lcoeffs) < deg
    push!(lcoeffs, zero(K))
  end
  push!(lcoeffs, one(K))
  pmat = identity_matrix(K, Int(deg))
  L, X = Native.finite_field(polynomial(K, lcoeffs), Symbol(:x, deg))
  set_standard_finite_field!(L)
  set_primitive_powers_in_tower_basis!(L, pmat)

  return L
end
function _extension_with_tower_basis(
  K::T,
  deg::IntegerUnion,
  lcoeffs::Vector,
  b::FinFieldElem,
) where {T<:FinField}
  @assert parent(b) === K

  dK = absolute_degree(K)
  if dK == 1
    println("_extension_with_tower_basis is running unoptimized...")
  end
  d = Int(dK * deg)
  F = prime_field(K)
  while length(lcoeffs) < deg
    push!(lcoeffs, zero(K))
  end
  push!(lcoeffs, one(K))
  pK = primitive_powers_in_tower_basis(K)

  # The idea is to collect (bX)^i mod f for 1 in 0..d*dK-1
  # and use this to compute the minimal polynomial of bX over F.
  # Should we just form the polynomial and compute "mod"???
  vec = zeros(F, d)
  vec[1] = one(F)
  v = zeros(K, Int(deg))
  v[1] = one(K)

  vecs = Vector{Vector{eltype(F)}}(undef, d)
  pols = Vector{Vector{eltype(F)}}(undef, d)
  pmat = zero_matrix(F, d, d)
  poly = Vector{eltype(F)}[]

  for i in 1:d+1
    # println("i: ", i, " vec: ", vec, " v: ", v)
    if i <= d
      pmat[i, :] = vec
    end

    poly = zeros(F, i)
    poly[end] = one(F)

    w = copy(vec)
    piv = findfirst(!iszero, w)
    # TODO : figure out the purpose of this loop and FIX it
    while piv !== nothing && isassigned(vecs, piv)
      x = -w[piv]
      if isassigned(pols, piv)
        # println("p: ", p, " piv ", piv, " pols[piv]: ", pols[piv])
        poly[1:length(pols[piv])] += x .* pols[piv]
      end
      w[piv:d] += x .* @view vecs[piv][piv:d]
      piv = findnext(!iszero, w, piv + 1)
    end
    # NOTE : exits the while loop when either piv == nothing, or vecs[piv] is undefined.
    #        what happens if piv == nothing???
    # println("Exiting loop with piv = ", piv)
    if i <= d
      # println(order(K), " ", i, " ", piv, " ", v, " ", w)
      x = inv(w[piv])
      poly .= x .* poly
      w .= x .* w
      pols[piv] = copy(poly)
      vecs[piv] = copy(w)


      # Multiply by vX and find the next vec in the Tower Basis
      v = collect(
        coefficients(
          mod(polynomial(K, pushfirst!(b .* v, zero(K))), polynomial(K, lcoeffs)),
        ),
      )

      while length(v) < deg
        push!(v, zero(K))
      end
      # println("new v:", v)
      vec = vcat(map(a -> coords(a) * pK, v)...)

    end
  end
  # Now p is the minimal polynomial over F
  # pmat gives the primitive powers in the tower basis for the new extension

  L, X = Native.finite_field(polynomial(F, poly), Symbol(:x, d))
  set_standard_finite_field!(L)
  set_primitive_powers_in_tower_basis!(L, pmat)

  return L

end


@doc raw"""
    standard_finite_field(p::Union{ZZRingElem, Integer}, n::Union{ZZRingElem, Integer}) -> FinField

Return a finite field of order $p^n$.

# Examples
```jldoctest
julia> standard_finite_field(3, 24)
Finite field of degree 24 over GF(3)
```
"""
function standard_finite_field(p::IntegerUnion, n::IntegerUnion)
  @req is_prime(p) "first argument must be a prime"
  F = Native.GF(p)
  set_standard_prime_field!(F)

  function _sff(N::Fac{ZZRingElem})
    # local m::ZZRingElem, k::IntegerUnion, nK::ZZRingElem, K::FinField, stn::ZZRingElem,
    #         n1::ZZRingElem, q1::ZZRingElem, l::Vector{ZZRingElem}, c::Vector{ZZRingElem}, b::FinFieldElem
    m, k = pop_largest_factor!(N)
    nK = evaluate(N)

    K = get_standard_extension!(F, nK) do
      _sff(N)
    end
    stn = steinitz_number_for_prime_degree(p, m, k)
    n1 = ZZ(m)^(k - 1)
    q1 = ZZ(p)^n1

    l = digits(stn, base = BigInt(q1))
    c = map(y -> element_from_steinitz_number(K, embed_steinitz(p, n1, nK, y)), l)

    d = divexact(nK, n1)
    b = element_from_steinitz_number(
      K,
      p^(findfirst(==(d), standard_monomial_degrees(nK)) - 1),
    )

    return _extension_with_tower_basis(K, m, c, b)
  end

  return get_standard_extension!(F, n) do
    N = factor(ZZ(n))
    return _sff(N)
  end
end

end

using .StandardFiniteFields
export standard_finite_field
