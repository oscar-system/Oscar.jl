
underlying_cycle(D::AbsWeilDivisor) = underlying_cycle(underlying_divisor(D))

### forwarding of all essential functionality
underlying_cycle(D::WeilDivisor) = D.C

### type getters 
scheme_type(D::WeilDivisor{S, U, V}) where{S, U, V} = S
scheme_type(::Type{WeilDivisor{S, U, V}}) where{S, U, V} = S
coefficient_ring_type(D::WeilDivisor{S, U, V}) where{S, U, V} = U
coefficient_ring_type(::Type{WeilDivisor{S, U, V}}) where{S, U, V} = U
coefficient_type(D::WeilDivisor{S, U, V}) where{S, U, V} = V
coefficient_type(::Type{WeilDivisor{S, U, V}}) where{S, U, V} = V

@doc raw"""
    WeilDivisor(X::CoveredScheme, R::Ring)

Return the zero `WeilDivisor` over `X` with coefficients 
in `R`.
"""
function WeilDivisor(X::AbsCoveredScheme, R::Ring)
  D = IdDict{AbsIdealSheaf, elem_type(R)}()
  return WeilDivisor(X, R, D, check=false)
end

function zero(W::WeilDivisor)
  return WeilDivisor(ambient_scheme(W), coefficient_ring(W))
end

# provide non-camelcase methods
@doc raw"""
    weil_divisor(X::AbsCoveredScheme, R::Ring) -> WeilDivisor

Return the zero weil divisor on `X` with coefficients in the ring `R`.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P, I);

julia> Ycov = covered_scheme(Y);

julia> weil_divisor(Ycov, QQ)
Zero weil divisor
  on scheme over QQ covered with 3 patches
with coefficients in rational field
```
"""
weil_divisor(X::AbsCoveredScheme, R::Ring) = WeilDivisor(X, R)

@doc raw"""
    WeilDivisor(I::AbsIdealSheaf)

Return the `WeilDivisor` ``D = 1 ⋅ I`` with coefficients
in ``ℤ`` for a sheaf of ideals ``I`` of pure codimension `1`.
"""
function WeilDivisor(I::AbsIdealSheaf; check::Bool=true)
  WeilDivisor(I, ZZ, check=check)
end

@doc raw"""
    weil_divisor(I::AbsIdealSheaf) -> WeilDivisor

Given an ideal sheaf `I` of pure codimension ``1``, return the weil divisor $D = 1 ⋅ I$ with
coefficients in the integer ring.

# Example
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> Y = proj(P);

julia> I = ideal([(x^3-y^2*z)]);

julia> II = IdealSheaf(Y, I);

julia> weil_divisor(II)
Effective weil divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  1 * sheaf of ideals
  
julia> JJ = II^2;

julia> D = weil_divisor(JJ)
Effective weil divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  1 * product of 2 ideal sheaves

julia> irreducible_decomposition(D)
Effective weil divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
given as the formal sum of
  2 * sheaf of prime ideals

```
"""
weil_divisor(I::AbsIdealSheaf; check::Bool=true) = WeilDivisor(I, check=check)

function WeilDivisor(I::AbsIdealSheaf, R::Ring; check::Bool=true)
  D = WeilDivisor(space(I), R)
  @check is_equidimensional(I) "ideal sheaf must be equidimensional"
  @check dim(space(I)) - dim(I) == 1 "components of a divisor must be of codimension one"
  setindex!(D, one(R), I; check)
  return D
end

@doc raw"""
    weil_divisor(I::AbsIdealSheaf, R::Ring; check::Bool=true)
    
Given an ideal sheaf `I` of pure codimension ``1`` and a ring `R`, return the weil divisor $D = 1 ⋅ I$ with
coefficients in `R`.
"""
weil_divisor(I::AbsIdealSheaf, R::Ring; check::Bool=true) = WeilDivisor(I, R, check=check)

### copy constructor
function copy(D::AbsWeilDivisor) 
  new_dict = IdDict{AbsIdealSheaf, elem_type(coefficient_ring_type(D))}()
  for I in keys(coefficient_dict(D))
    new_dict[I] = D[I]
  end
  return WeilDivisor(ambient_scheme(D), coefficient_ring(D), new_dict, check=false)
end

function irreducible_decomposition(D::AbsWeilDivisor)
  E = irreducible_decomposition(underlying_cycle(D))
  return WeilDivisor(E; check=false)
end

# If we know something about the Weil divisor, we write it! Always good to have
# relevant information for free
function Base.show(io::IO, D::AbsWeilDivisor)
  io = pretty(io)
  X = ambient_scheme(D)
  if get(io, :show_semi_compact, false)
    cov = Oscar._covering_for_printing(io, X)
    _show_semi_compact(io, D, cov)
  else
    C = underlying_cycle(D)
    eff = all(i >= 0 for i in collect(values(coefficient_dict(C))))
    prim = eff && get_attribute(D, :is_prime, false)
    if has_name(D)
      print(io, name(D))
    elseif is_terse(io)
      print(io, "Weil divisor")
    # if the divisor is prime and the ideal sheaf has a name print that
    elseif length(components(D)) == 1 && has_attribute(first(components(D)), :name)
      I = first(components(D))
      I_name = get_attribute(I, :name)
      print(io, Lowercase(), I_name)
    elseif length(components(D)) == 0
      print(io, "Zero weil divisor on ", Lowercase(),  X)
    elseif eff
      if prim
        print(io, "Prime Weil divisor on ", Lowercase(), X)
      else
        print(io, "Effective Weil divisor on ", Lowercase(), X)
      end
    else
      print(io, "Weil divisor on ", Lowercase(), X)
    end
  end
end

# Used in nested printing, where we assume that the associated scheme is already
# printed in the nest - we keep track of the good covering `cov` to describe
# everything consistently.
function _show_semi_compact(io::IO, D::AbsWeilDivisor, cov::Covering)
  io = pretty(io)
  X = ambient_scheme(D)
  C = underlying_cycle(D)
  eff = all(i >= 0 for i in collect(values(coefficient_dict(C))))
  prim = eff && get_attribute(D, :is_prime, false)
  if has_name(D)
    print(io, name(D))
  elseif length(components(D)) == 0
    print(io, "Zero weil divisor on ", Lowercase())
  elseif eff
    if prim
      print(io, "Prime weil divisor on ", Lowercase())
    else
      print(io, "Effective Weil divisor on ", Lowercase())
    end
  else
    print(io, "Weil divisor on ", Lowercase())
  end
  show(IOContext(io, :show_semi_compact => true, :covering => cov), X)
end

# Take care of some offsets to make sure that the coefficients are all aligned
# on the right.
function Base.show(io::IO, ::MIME"text/plain", D::AbsWeilDivisor)
  io = pretty(io)
  X = ambient_scheme(D)
  cov = Oscar._covering_for_printing(io, X)
  C = underlying_cycle(D)
  eff = all(i >= 0 for i in collect(values(coefficient_dict(C))))
  prim = eff && get_attribute(D, :is_prime, false)
  if length(components(C)) == 0
    print(io, "Zero weil divisor")
  else
    if eff
      if prim
        print(io, "Prime weil divisor")
      else
        print(io, "Effective weil divisor")
      end
    else
      print(io, "Weil divisor")
    end
  end
  if has_name(D)
    print(io, " ", name(D))
  end
  println(io)
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :covering => cov), X)
  println(io, Dedent())
  print(io, "with coefficients in ", Lowercase(), coefficient_ring(C))
  if length(components(C)) != 0
    println(io)
    print(io, Dedent(), "given as the formal sum of")
    print(io, Indent())
    co_str = String["$(C[I])" for I in components(C)]
    k = max(length.(co_str)...)
    for i in 1:length(components(C))
      println(io)
      I = components(C)[i]
      kI = length(co_str[i])
      print(io, " "^(k-kI)*"$(C[I]) * ")
      print(io, Indent(), Lowercase())
      show(IOContext(io, :show_scheme => false), I)
      print(io, Dedent())
    end
  end
end

function Base.:+(D::AbsWeilDivisor, E::AbsWeilDivisor) 
  return underlying_divisor(D) + underlying_divisor(E)
end

function Base.:+(D::WeilDivisor, E::AbsWeilDivisor) 
  return D + underlying_divisor(E)
end

function Base.:+(D::AbsWeilDivisor, E::WeilDivisor) 
  return underlying_divisor(D) + E
end

function Base.:+(D::WeilDivisor, E::WeilDivisor) 
  return WeilDivisor(underlying_cycle(D) + underlying_cycle(E), check=false)
end



function Base.:-(D::AbsWeilDivisor)
  return -underlying_divisor(D)
end

function Base.:-(D::WeilDivisor)
  return WeilDivisor(-underlying_cycle(D), check=false)
end


-(D::AbsWeilDivisor, E::AbsWeilDivisor) = D + (-E)

function Base.:*(a::RingElem, E::AbsWeilDivisor)
  return a*underlying_divisor(E)
end

# Method ambiguity requires the following two methods:
function Base.:*(a::ZZRingElem, E::AbsWeilDivisor)
  return a*underlying_divisor(E)
end

function Base.:*(a::ZZRingElem, E::WeilDivisor)
  return WeilDivisor(a*underlying_cycle(E), check=false)
end

function Base.:*(a::RingElem, E::WeilDivisor)
  return WeilDivisor(a*underlying_cycle(E), check=false)
end

Base.:*(a::T, E::AbsWeilDivisor) where {T<:IntegerUnion} = coefficient_ring(E)(a)*E
# method ambiguity requires us to also implement the following:
Base.:*(a::Int, E::AbsWeilDivisor) = coefficient_ring(E)(a)*E

Base.:+(D::AbsWeilDivisor, I::AbsIdealSheaf) = D + WeilDivisor(I)

function ==(D::AbsWeilDivisor, E::AbsWeilDivisor) 
  return underlying_divisor(D) == underlying_divisor(E)
end

function ==(D::WeilDivisor, E::AbsWeilDivisor) 
  return D == underlying_divisor(E)
end

function ==(D::AbsWeilDivisor, E::WeilDivisor) 
  return underlying_divisor(D) == E
end

function ==(D::WeilDivisor, E::WeilDivisor) 
  return underlying_cycle(D) == underlying_cycle(E)
end

@doc raw"""
    intersect(D::AbsWeilDivisor, E::AbsWeilDivisor; covering::Covering=default_covering(ambient_scheme(D)))

Return the intersection number of the the Weil divisors `D` and `E`
on a complete smooth surface as defined in [Har77](@cite).

# Input 
The optional keyword argument `covering` specifies the covering to be used for the computation.
"""
function intersect(D::AbsWeilDivisor, E::AbsWeilDivisor;
    covering::Covering=default_covering(ambient_scheme(D))
  )
  X = ambient_scheme(D)
  @req dim(X) == 2 "intersection of Weil divisors is only implemented for surfaces."
  X === ambient_scheme(E) || error("divisors do not live on the same scheme")
  R = coefficient_ring(D)
  R === coefficient_ring(E) || error("divisors do not have the same coefficient ring")
  result = zero(R)
  for c1 in components(D)
    a1 = D[c1]
    for c2 in components(E)
      a2 = E[c2]
      if c1 === c2
        result = result + a1*a2*_self_intersection(c1)
      else
        I = c1 + c2
        if !has_dimension_leq_zero(I) # potentially faster for localized ideals
          if c1 == c2
            result = result + a1*a2*_self_intersection(c1)
          else
            error("self intersection unknown")
          end
        else
          result = result + a1 * a2 * colength(I, covering=covering)
        end
      end
    end
  end
  return result
end


function pushforward(inc::CoveredClosedEmbedding, W::WeilDivisor)
  X = domain(inc)
  Y = codomain(inc)
  X === ambient_scheme(W) || error("divisor not defined on the domain")
  kk = coefficient_ring(W)
  ideal_dict = IdDict{AbsIdealSheaf, elem_type(kk)}()
  for I in components(W)
    pfI = pushforward(inc)(I)
    ideal_dict[pfI] = W[I]
  end
  return WeilDivisor(Y, kk, ideal_dict, check=false)
end



"""
    _self_intersection(I::AbsIdealSheaf) -> Integer

For ``I`` a sheaf of pure codimension ``1`` on a surface,
return the self-intersection of ``I`` viewed as a Weil-Divisor.
"""
function _self_intersection(I::AbsIdealSheaf)
  has_attribute(I, :_self_intersection) || error("self intersection unknown")
  return get_attribute(I, :_self_intersection)::Int
end

function _is_known_to_be_one(I::AbsIdealSheaf, U::AbsAffineScheme;
    dec_inf::Vector = elem_type(OO(U))[]
  )
  @vprintln :Divisors 5 "checking triviality of $(I) on $(U)"
  K = ideal(OO(U), dec_inf)
  if U in keys(object_cache(I))
    i = I(U)
    has_attribute(i, :is_one) && get_attribute(i, :is_one) === true && return true
    one(OO(U)) in gens(i) && return true
    is_one(K + i) && return true
  end
  j = cheap_sub_ideal(I, U)
  is_one(j + K) && return true
  return false
end

function _is_known_to_be_one(I::SumIdealSheaf, U::AbsAffineScheme;
    dec_inf::Vector = elem_type(OO(U))[]
  )
  @vprintln :Divisors 5 "triviality on $(U)"
  K = ideal(OO(U), dec_inf)
  for J in summands(I)
    if U in keys(object_cache(J))
      j = J(U)
      has_attribute(j, :is_one) && get_attribute(j, :is_one) === true && return true
      one(OO(U)) in gens(j) && return true
      is_one(K + j) && return true
    end
  end
  if U in keys(object_cache(I))
    i = I(U)
    has_attribute(i, :is_one) && get_attribute(i, :is_one) === true && return true
    one(OO(U)) in gens(i) && return true
    is_one(K + i) && return true
  end
  j = cheap_sub_ideal(I, U)
  is_one(j + K) && return true
  return false
end


@doc raw"""
    is_in_linear_system(f::VarietyFunctionFieldElem, D::WeilDivisor; regular_on_complement::Bool=true, check::Bool=true) -> Bool

Return whether the rational function `f` is in the linear system ``|D|``, i.e. if $(f) + D \geq 0$.

# Input 
- `regular_on_complement` -- set to `true` if `f` is regular on the complement of the support of `D`. 
"""
function is_in_linear_system(f::VarietyFunctionFieldElem, D::AbsWeilDivisor; regular_on_complement::Bool=false, check::Bool=true)
  X = ambient_scheme(D) 
  X === variety(parent(f)) || error("schemes not compatible")
  C = simplified_covering(X)
  for I in components(D)
    @check is_prime(I) "components of the divisor must be prime"
    order_of_vanishing(f, I, check=false) >= -D[I] || return false
  end
  regular_on_complement && return true
  for U in patches(C)
    # we have to check that f[U] has no poles outside the support of D[U]
    J = intersect([J(U) for J in components(D)])
    incH = ClosedEmbedding(U, J)
    W = complement(incH) # This is a AffineSchemeOpenSubscheme
    is_regular(f, W) || return false
  end
  return true
end


function Base.show(io::IO, L::LinearSystem)
  if is_terse(io)
    print(io, "Linear system")
  else
    io = pretty(io)
    print(io, "Linear system of ", Lowercase(), weil_divisor(L))
  end
end

function Base.show(io::IO, ::MIME"text/plain", L::LinearSystem)
  io = pretty(io)
  X = scheme(L)
  cov = default_covering(X)
  println(io, "Linear system")
  print(io, Indent(), "of ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => cov), weil_divisor(L))
  gg = gens(L)
  if length(gg) > 0
    println(io)
    print(io, Dedent(), "generated by")
    print(io, Indent())
    ll = String["$(representative(f))" for f in gg]
    k = max(length.(ll)...)
    offset = Int[k-length(s) for s in ll]
    for i in 1:length(gg)
      f = gg[i]
      println(io)
      print(io, Lowercase())
      show(IOContext(io, :show_semi_compact => true, :covering => cov, :offset => offset[i]), f)
    end
    print(io, Dedent())
  end
end

@doc raw"""
    linear_system(f::Vector, D::AbsWeilDivisor; check::Bool=true)

Return the linear system ``L`` generated by rational functions ``f₁,…,fᵣ ∈ K(X)``
with $L \subseteq |D|$ for `D` on a variety `X`. If `check` is set,
confirm that $L \subseteq |D|$.
"""
linear_system(f::Vector, D::AbsWeilDivisor; check::Bool=true) = LinearSystem(f, D, check=check)

### essential getters 
@doc raw"""
    weil_divisor(L::LinearSystem)

Return the divisor $D$ of the linear system $L = |D|$.
"""
function weil_divisor(L::LinearSystem) 
  return L.D
end

gens(L::LinearSystem) = L.f
number_of_generators(L::LinearSystem) = length(L.f)
gen(L::LinearSystem, i::Int) = L.f[i]

@doc raw"""
    variety(L::LinearSystem)

Return the variety on which `L` is defined.
"""
variety(L::LinearSystem) = ambient_scheme(weil_divisor(L))
# an alias for the user's convenience 
scheme(L::LinearSystem) = variety(L)

@doc raw"""
    subsystem(L::LinearSystem, D::AbsWeilDivisor) -> LinearSystem, MatElem

Given a linear system $L = |E|$ and a divisor $D \leq E$ compute $|D|$
and the matrix representing the inclusion $|D| \hookrightarrow |E|$
with respect to the given bases of both systems.
"""
function subsystem(L::LinearSystem, D::AbsWeilDivisor)
  E = weil_divisor(L)
  @req D <= E "input does not define a subsystem"
  Lnew = L
  T = identity_matrix(base_ring(scheme(L)), ngens(L))
  for P in components(E)
    if coeff(D,P) == coeff(E,P)
      continue
    end
    Lnew, Tnew = _subsystem(Lnew, P, -coeff(D,P))
    T = Tnew*T
  end
  return Lnew, T
end


@doc raw"""
    _subsystem(L::LinearSystem, P::AbsIdealSheaf, n) -> LinearSystem

Given a linear system ``L = |D|``, a sheaf of prime ideals `P` 
and an integer `n`, return a pair ``(K, A)`` consisting
of the subsystem of elements in ``|D|`` that vanish to order at least n at ``P``.
The matrix ``A`` for its inclusion into ``L`` on the given set
of generators.
"""
function _subsystem(L::LinearSystem, P::AbsIdealSheaf, n; covering::Covering=simplified_covering(scheme(P)))
  # find one chart in which P is supported
  if coeff(weil_divisor(L),P) == -n
    return L, identity_matrix(ZZ, length(gens(L)))
  end
  X = variety(L)
  X === space(P) || error("input incompatible")
  C = covering

  U = _find_good_representative_chart(P; covering)

  # Now U it is.
  # Assemble the local representatives
  R = ambient_coordinate_ring(U)
  loc_rep = [g[U] for g in gens(L)]
  if length(loc_rep) == 0
    common_denominator = R(1)
  else
    common_denominator = lcm([denominator(g) for g in loc_rep])
  end
  numerators = [numerator(g)*divexact(common_denominator, denominator(g)) for g in loc_rep]

  # compute a symbolic power
  RP, _ = localization(OO(U), complement_of_prime_ideal(saturated_ideal(P(U))))
  PP = RP(prime_ideal(inverted_set(RP)))
  K = function_field(X)

  denom_mult = order_of_vanishing(K(common_denominator), P, check=false)
  #denom_mult = (_minimal_power_such_that(PP, I -> !(RP(common_denominator) in I))[1])-1
  w = n + denom_mult # Adjust!
  if w < 0
    pPw = ideal(R, one(R))
  else
    pPw = saturated_ideal(PP^w) # the symbolic power
  end
  # reduce the numerators modulo P^(w)
  images = elem_type(R)[]
  for a in numerators
    push!(images, normal_form(a, pPw))
  end
  # collect a monomial basis in which to represent the results
  all_mons = elem_type(R)[]
  for b in images
    all_mons = vcat(all_mons, [m for m in monomials(b) if !(b in all_mons)])
  end

  kk = base_ring(X)
  A = zero_matrix(kk, ngens(L), length(all_mons))
  for i in 1:ngens(L)
    for (c, m) in zip(coefficients(images[i]), monomials(images[i]))
      k = findfirst(==(m), all_mons)
      A[i, k] = c
    end
  end
  K = kernel(A; side = :left)
  r = nrows(K)
  new_gens = [sum([K[i,j]*gen(L, j) for j in 1:ncols(K)]) for i in 1:r]
  W = weil_divisor(L)
  PW = WeilDivisor(P, check=false)
  k = coeff(W,P)
  D = W + (min(-n,k)-k)*PW
  return LinearSystem(new_gens, D, check=false), K[1:r,:]
end

function subsystem(L::LinearSystem, P::AbsWeilDivisor, n::Int; check::Bool=true)
  @check is_prime(P) "P must be a prime divisor"
  I = components(P)[1]
  return subsystem(L, I, n)
end


@doc raw"""
    order_of_vanishing(f::VarietyFunctionFieldElem, D::AbsWeilDivisor; check::Bool=true)

Return the order of vanishing of the rational function `f` on the prime divisor `D`.
"""
function order_of_vanishing(
    f::VarietyFunctionFieldElem,
    D::AbsWeilDivisor;
    check::Bool=true
  )
  @check is_prime(D) || error("divisor must be prime")
  I = components(D)[1]
  return order_of_vanishing(f, I, check=false)
end
