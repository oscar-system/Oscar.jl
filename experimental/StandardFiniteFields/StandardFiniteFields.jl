using Oscar

IntegerUnion   = Union{Integer, ZZRingElem}
PrimeField     = Union{Nemo.fpField, Nemo.FpField}
PrimeFieldElem = Union{fpFieldElem, FpFieldElem}
PrimeFieldMatrix = Union{FpMatrix, fpMatrix}

# NOTE: These give missing features to OSCAR/Nemo that will likely be added in the near future.

# TODO : Should be fixed in Nemo
function (k::Nemo.FpField)(a::Vector)
  @assert length(a) == 1
  return k(a[1])
end
function (k::FqPolyRepField)(a::Vector)
  return k(polynomial(GF(ZZ(characteristic(k))), a))
end


# TODO : Should be fixed in Hecke for prime fields
function coords(x::FinFieldElem)
    return absolute_coordinates(x)
end
function coords(x::PrimeFieldElem)
    return [x]
end

# TODO : this should be pushed to Nemo.jl/src/flint/gfp_fmpz_mat.jl
# BUG the fact that this is missing is a major performance problem
import Base: inv
using FLINT_jll
const libflint = FLINT_jll.libflint
function inv(a::FpMatrix)
  !is_square(a) && error("Matrix must be a square matrix")
  z = similar(a)
  r = ccall((:fmpz_mod_mat_inv, libflint), Int,
          (Ref{FpMatrix}, Ref{FpMatrix}), z, a)
  !Bool(r) && error("Matrix not invertible")
  return z
end

# TODO : this just makes sense for writing general code
import Base: big
big(a::ZZRingElem) = a


# _attributes = [
    # :is_standard_finite_field      - bool
    # :is_standard_prime_field       - bool
    # :primitive_powers_in_tower_basis - Matrix{FinFieldElem}
    # :tower_basis                  - Matrix{FFieldElem} = inv(:primitive_powers_in_tower_basis)
    # NOTE I guess we used ZZ ring elements for these instead of Int; does it matter?
    # :steinitz_prime_degree         - Dict{Int, Dict{Int, Int}}
    # :standard_extensions          - Dict{Int, FinField} ]
function set_standard_prime_field!(F::PrimeField)
    get_attribute!(F, :is_standard_prime_field) do
        set_attribute!(F, :standard_extensions, Dict{ZZRingElem, FinField}(1 => F))
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
function set_steinitz_prime_degree!(F::FinField, r::IntegerUnion, k::IntegerUnion, nr::IntegerUnion)
    spd = get_attribute!(F, :steinitz_prime_degree, Dict{ZZRingElem, Dict{ZZRingElem, ZZRingElem}}())
    spdr = get!(spd, r, Dict{ZZRingElem, ZZRingElem})
    spdr[k] = nr
end
function set_steinitz_prime_degree!(f::Function, F::FinField, r::IntegerUnion, k::IntegerUnion)
    spd = get_attribute!(F, :steinitz_prime_degree, Dict{ZZRingElem, Dict{ZZRingElem, ZZRingElem}}())
    spdr = get!(spd, r, Dict{ZZRingElem, ZZRingElem})
    spd[r][k] = f()
end
function set_standard_extension!(F::PrimeField, n::IntegerUnion, K::FinField)
    ext = get_attribute!(F, :standard_extensions, Dict{ZZRingElem, FinField}())
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
function get_steinitz_prime_degree!(f::Function, F::FinField, r::IntegerUnion, k::IntegerUnion)
    spd = get_attribute!(F, :steinitz_prime_degree, Dict{ZZRingElem, Dict{ZZRingElem, ZZRingElem}}() )
    spdr = get!(spd, r, Dict{ZZRingElem, ZZRingElem}())
    get!(spdr, k, f())
end
function get_standard_extensions(F::PrimeField)
  get_attribute(F, :standard_extensions, nothing)
end
function get_standard_extensions!(F::PrimeField)
  get_attribute(F, :standard_extensions, Dict{ZZRingElem, FinField}())
end
function get_standard_extension(F::PrimeField, k::IntegerUnion)
  ext = get_attribute(F, :standard_extensions, nothing)
  get(ext, k, nothing)
end
function get_standard_extension!(F::PrimeField, k::IntegerUnion, L::FinField)
  ext = get_attribute!(F, :standard_extensions, Dict{ZZRingElem, FinField}())
  get!(ext, k, L)
end
function get_standard_extension!(f::Function, F::PrimeField, k::IntegerUnion)
  ext = get_attribute!(F, :standard_extensions, Dict{ZZRingElem, FinField}())
  get!(f, ext, k)
end


# TODO: Lübeck speeds this up by caching triples [q,m,a] resulting from this
function standard_affine_shift(q::IntegerUnion, i::IntegerUnion)
    m = div(4*q, 5)
    while gcd(m,q) != 1
        m -= 1
    end
    a = div(2*q, 3)
    return mod((m*i + a), q)
end

# Given a field F and Steinitz number n, give the corresponding field element.
# we REQUIRE that F is a standard finite field TODO: using @assert?
function element_from_steinitz_number(F::PrimeField, n::IntegerUnion)
    return F(n)
end
function element_from_steinitz_number(F::FinField, n::IntegerUnion)
    p = characteristic(F)
    q = order(F)
    if n < 0 || n > q
        error("We need to have 0 <= n <= q")
    end
    if n == 0
        return zero(F)
    else
        # this forms a linear combo of F.towervasis rows using vectorrep as coefficients,
        # and then convert this vector to an element of F.
        vectorrep = digits(n, base = Int(p))
        return F(vectorrep * @view tower_basis(F)[1:length(vectorrep), :])
    end
end

# Returns an element a in F that is NOT an rth root of unity
# we REQUIRE that F is a standard finite field TODO: using @assert?
function non_rth_root(F::FinField, r::IntegerUnion)
    @assert is_standard_finite_field(F) || is_standard_prime_field(F)
    q = order(F)
    if mod(q-1, r) == 0
        i = 0
        a = zero(F)
        k = div(q-1,r)
        while iszero(a) || isone(a^k)
            i += 1
            a = element_from_steinitz_number(F, standard_affine_shift(q, i))
        end
        return a
    else
        return nothing
    end
end

function standard_irreducible_coefficient_list(F::FinField, r::IntegerUnion, a::FinFieldElem)
    q = order(F)
    l = zeros(F, Int(r)+1)
    l[Int(r)+1] = one(F)
    l[1] = a
    l[2] = one(F)
    # inc is the expected number of nonzero coefficients
    inc = 1
    while q^inc < 2*r
        inc += 1
    end
    # allowing non-zero coeffs up to position d
    # after every r attempts allow inc more non-zero coeffs
    d = 0
    count = 0
    qq = 0
    while !is_irreducible(polynomial(F, l))
        if mod(count, r) == 0 && d < r-1
            d += inc
            if d >= r
                d = r-1
            end
            qq = q^(d-1)
        end
        # q can be very very large so Int is not big enough...
        st = digits(standard_affine_shift(qq,count), base = BigInt(q), pad = d-1)
        # TODO: we can remove this while loop when digits bug for n = 0 is fixed
        while length(st) < d-1
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
# TODO Maybe want to ensure this always returns a BigInt?
function steinitz_number_for_prime_degree(p::IntegerUnion, r::IntegerUnion, k::IntegerUnion)
    Fp = standard_finite_field(p,1)

    get_steinitz_prime_degree!(Fp, r, k) do
        # now we need to create the polynomial depending on the prime r
        if r == p
            # Artin-Schreier case
            # k = 1 we get [(Xr[1])^p - (Xr[1]) -1]
            # k > 1 we get (Xr[k])^p - (Xr[k]) - (prod(Xr[j] : j in [1..k-1]))^(p-1))
            q = big(p)^(p^(k-1))
            return (p-1)*(q + div(q,p))
        elseif r == 2 && mod(p,4) == 3
            if k == 1
                # (Xr[1])^2 +1
                return 1
            elseif k == 2
                a = non_rth_root(standard_finite_field(p,2), r)
                # Xr[2]^2 -a
                return steinitz_number(-a)
            else
                # Xr[i]^2 - Xr[i-1]
                return (p-1)*big(p)^(r^(k-2))
            end
        elseif r == 2
            if k == 1
                # Xr[1]^2 -a
                a = non_rth_root(standard_finite_field(p, 1), r)
                return steinitz_number(-a)
            else
                # Xr[j]^r - Xr[j-1]
                return (p-1)*big(p)^(r^(k-2))
            end
        else
            # Here we use pseudo-random polynomials...
            F = standard_finite_field(p, r^(k-1))
            if k == 1
                a = -one(F)
            else
                a = -gen(F)
            end
            l = standard_irreducible_coefficient_list(F,r,a)
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
  return lift(x)
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
# just return degrees
function standard_monomial_degrees(n::IntegerUnion)
    if n == 1
        return [1]
    end
    # need the largest prime factor a of n
    nfactorization = factor(ZZ(n))
    nfactors = sort([r for (r,e) in nfactorization])
    a = Int(nfactors[end])
    res = standard_monomial_degrees(div(n,a))
    k = a^nfactorization[a]
    new = map( x -> lcm(x, k), res)
    for i = 1:a-1
        append!(res, new)
    end
    return res
end
# map of monomials for degree n -> monomials of degree m by positions
function standard_monomial_map(n::IntegerUnion, m::IntegerUnion)
    d = standard_monomial_degrees(m)
    return [i for i = 1:length(d) if mod(n, d[i]) == 0]
end

# Embed an element x of Fp^n into Fp^m by Steinitz numbers
# where nr = steinitz_number(Fp^n, x)
# I hate hate hate these variable names copied (mostly) from Lübeck
function embed_steinitz(p::IntegerUnion, n::IntegerUnion, m::IntegerUnion, nr::IntegerUnion)
    if n == m || iszero(nr)
        return nr
    end
    l = digits(nr, base = Int(p))
    m = @view standard_monomial_map(n,m)[1:length(l)]
    c = zeros(ZZRingElem, m[end])
    c[m] = l
    return evaluate(polynomial(ZZ, c), p)
end


# Given a field K, we construct an extension L with [L:K] = deg
# We use the irreducible polynomial f = X^deg  + g(X)
#    where lcoeffs are the coefficients of g(X).
# We assume b is a generator for K, and so bX will be a generator for L
function _extension_with_tower_basis(K::PrimeField, deg::IntegerUnion, lcoeffs::Vector, b::PrimeFieldElem)
    @assert parent(b) === K

    while length(lcoeffs) < deg
        push!(lcoeffs, zero(K))
    end
    push!(lcoeffs, one(K))
    pmat = identity_matrix(K, Int(deg))
    vname = "x" * string(deg)
    L, X = FiniteField(polynomial(K, lcoeffs), vname)
    set_standard_finite_field!(L)
    set_primitive_powers_in_tower_basis!(L, pmat)

    return L
end
function _extension_with_tower_basis(K::T, deg::IntegerUnion, lcoeffs::Vector, b::FinFieldElem) where T<:FinField
    @assert parent(b) === K

    dK = absolute_degree(K)
    if dK == 1 then
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
            pmat[i, : ] = vec
        end

        poly = zeros(F, i)
        poly[end] = one(F)

        w = copy(vec)
        piv = findfirst(!iszero, w)
        # TODO : figure out the purpose of this loop and FIX it
        while piv != nothing && isassigned(vecs, piv)
            x = -w[piv]
            if isassigned(pols, piv)
                # println("p: ", p, " piv ", piv, " pols[piv]: ", pols[piv])
                poly[1:length(pols[piv])] += x .* pols[piv]
            end
            w[piv:d] += x .* @view vecs[piv][piv:d]
            piv = findnext(!iszero, w, piv+1)
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
            v = collect(coefficients(mod(polynomial(K, pushfirst!(b .* v, zero(K))),
                                         polynomial(K, lcoeffs))))

            while length(v) < deg
                push!(v, zero(K))
            end
            # println("new v:", v)
            vec = vcat(map( a -> coords(a) * pK, v)...)

        end
    end
    # Now p is the minimal polynomial over F
    # pmat gives the primitive powers in the tower basis for the new extension

    vname = "x" * string(d)
    L, X = FiniteField(polynomial(F, poly), vname)
    set_standard_finite_field!(L)
    set_primitive_powers_in_tower_basis!(L, pmat)

    return L

end


# TODO: this should work also if p is an integer
function standard_finite_field(p::T, n::IntegerUnion) where T<:IntegerUnion
    if !isprime(p)
        error()
    end
    F = GF(p)
    set_standard_prime_field!(F)
    get_standard_extension!(F, n) do
      nfactorization = factor(ZZ(n));
      nfactors = sort([r for (r,e) in nfactorization]);
      lastfactor = nfactors[end]
      nK = div(n,lastfactor)
      K = standard_finite_field(p, nK)

      stn = steinitz_number_for_prime_degree(p, Int(lastfactor), nfactorization[lastfactor])
      n1 = big(lastfactor)^(nfactorization[lastfactor]-1)
      # BUG this overflows if p is an Int64
      q1 = big(p)^T(n1)

      # for each element y in this list, we want to
      # 1. call EmbedSteinitz(p, n1, nK, y)
      # 2. this should give a number, we want to use ElementSteinitzNumber to get an element of K.
      l = digits(stn, base = BigInt(q1))
      c = map(y -> element_from_steinitz_number(K, embed_steinitz(p, n1, nK, y)), l)
      b = element_from_steinitz_number(K, p^( findfirst(x -> x == div(nK, n1), standard_monomial_degrees(nK))-1))

      _extension_with_tower_basis(K, lastfactor, c, b)
    end
end
