function (C::EffectiveCartierDivisor)(U::AbsAffineScheme)
  return gens(C.I(U))
end

iszero(C::EffectiveCartierDivisor) = isone(ideal_sheaf(C))

@doc raw"""
    ideal_sheaf(C::EffectiveCartierDivisor)
    
Return the sheaf of ideals $\mathcal{I}_C \subseteq \mathcal{O}_X$ representing `C`.
"""
ideal_sheaf(C::EffectiveCartierDivisor) = C.I

@doc raw"""
    ambient_scheme(C::EffectiveCartierDivisor)
    
Return the ambient scheme containing `C`. 
"""
ambient_scheme(C::EffectiveCartierDivisor) = C.X

@doc raw"""
    trivializing_covering(C::EffectiveCartierDivisor)
    
Return the trivializing covering of the effective Cartier divisor `C`.

A covering $(U_i)_{i \in I}$ is called trivializing for $C$ if 
$C(U_i)$ is principal for all $i \in I$.
"""
trivializing_covering(C::EffectiveCartierDivisor) = C.C

function EffectiveCartierDivisor(I::AbsIdealSheaf; 
    trivializing_covering::Covering=default_covering(scheme(I)),
    check::Bool=true
  )
  X = scheme(I)
  eq_dict = IdDict{AbsAffineScheme, RingElem}()
  for U in patches(trivializing_covering)
    isone(ngens(I(U))) || error("ideal sheaf is not principal on the given covering")
    eq_dict[U] = first(gens(I(U)))
  end
  return EffectiveCartierDivisor(X, eq_dict, trivializing_covering=trivializing_covering, check=check)
end

@doc raw"""
    ambient_scheme(C::CartierDivisor)
    
Return the ambient scheme containing `C`. 
"""
ambient_scheme(C::CartierDivisor) = C.X

@doc raw"""
    coefficient_ring(C::CartierDivisor)
    
Return the ring of coefficients of `C`.
"""
coefficient_ring(C::CartierDivisor) = C.R

coefficient_dict(C::CartierDivisor) = C.coeff_dict
getindex(C::CartierDivisor, k::EffectiveCartierDivisor) = coefficient_dict(C)[k]

@doc raw"""
    components(C::CartierDivisor)
    
Return a list of effective Cartier divisors $C_i$ such that $C$ is a linear combination of the $C_i$.
"""
components(C::CartierDivisor) = collect(keys(coefficient_dict(C)))

function +(C::CartierDivisor, D::CartierDivisor) 
  ambient_scheme(C) === ambient_scheme(D) || error("divisors must be defined over the same scheme")
  coefficient_ring(C) === coefficient_ring(D) || error("divisors must have the same coefficient rings")
  R = coefficient_ring(C)
  coeff_dict = IdDict{EffectiveCartierDivisor, elem_type(R)}()
  for k in keys(coefficient_dict(C))
    coeff_dict[k] = C[k]
  end
  for k in keys(coefficient_dict(D))
    if haskey(coeff_dict, k)
      c = coeff_dict[k] + D[k]
      if iszero(c)
        delete!(coeff_dict, k)
      else
        coeff_dict[k] = c
      end
    else
      coeff_dict[k] = D[k]
    end
  end
  return CartierDivisor(ambient_scheme(C), coefficient_ring(C), coeff_dict)
end

function -(C::CartierDivisor, E::EffectiveCartierDivisor)
  return C - 1*E
end

function +(C::CartierDivisor, D::EffectiveCartierDivisor) 
  return C + CartierDivisor(D)
end

function +(C::EffectiveCartierDivisor, D::EffectiveCartierDivisor) 
  return CartierDivisor(C) + CartierDivisor(D)
end

function +(C::EffectiveCartierDivisor, D::CartierDivisor) 
  return CartierDivisor(C) + D
end

zero(D::CartierDivisor) = CartierDivisor(ambient_scheme(D),coefficient_ring(D))

function *(a::RingElem, C::CartierDivisor)
  parent(a) === coefficient_ring(C) || return coefficient_ring(C)(a)*C
  coeff_dict = IdDict{EffectiveCartierDivisor, typeof(a)}()
  for k in keys(coefficient_dict(C))
    c = a*C[k]
    if iszero(c)
      # do nothing; forget about the generator
    else
      coeff_dict[k] = c
    end
  end
  return CartierDivisor(ambient_scheme(C), coefficient_ring(C), coeff_dict)
end

function *(a::Integer, C::CartierDivisor)
  return coefficient_ring(C)(a)*C
end

function -(C::CartierDivisor, D::CartierDivisor) 
  return C + (-one(coefficient_ring(D)))*D
end

function iszero(C::CartierDivisor)
  iszero(length(keys(coefficient_dict(C)))) && return true
  all(iszero, values(coefficient_dict(C))) && return true
  all(iszero, keys(coefficient_dict(C))) && return true
  
  # write C = P - M with P and M effective.
  # TODO: Multiplying everything together is quick and dirty
  P = sum(ai*Ci for (Ci,ai) in coefficient_dict(C) if ai>0; init=zero(C))
  M = sum(-ai*Ci for (Ci,ai) in coefficient_dict(C) if ai<0; init=zero(C))
  if iszero(length(coefficient_dict(P))) || iszero(length(coefficient_dict(M)))
      # we know that there is at least one non-zero summand and no cancellation
      return false
  end
  IP = prod(ideal_sheaf(Ci)^ai for (Ci,ai) in coefficient_dict(P))
  IM = prod(ideal_sheaf(Ci)^ai for (Ci,ai) in coefficient_dict(M))
  return IP==IM
end

@doc raw"""
    cartier_divisor(E::EffectiveCartierDivisor) -> CartierDivisor

Convert an `EffectiveCartierDivisor` into a `CartierDivisor` with
coefficient $1$ in the ring of integers.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P);

julia> II = IdealSheaf(Y, I);

julia> E = effective_cartier_divisor(II)
Effective cartier divisor
  on scheme over QQ covered with 3 patches
    1: [(y//x), (z//x)]   affine 2-space
    2: [(x//y), (z//y)]   affine 2-space
    3: [(x//z), (y//z)]   affine 2-space
defined by
  sheaf of ideals with restrictions
    1: Ideal (-(y//x)^2*(z//x) + 1)
    2: Ideal ((x//y)^3 - (z//y))
    3: Ideal ((x//z)^3 - (y//z)^2)

julia> cartier_divisor(E)
Cartier divisor
  on scheme over QQ covered with 3 patches
with coefficients in integer ring
defined by the formal sum of
  1 * effective cartier divisor on scheme over QQ covered with 3 patches
```
"""
cartier_divisor(E::EffectiveCartierDivisor) = CartierDivisor(E)

function CartierDivisor(C::EffectiveCartierDivisor)
  return CartierDivisor(ambient_scheme(C), ZZ, IdDict([C => one(ZZ)]))
end

function CartierDivisor(X::AbsCoveredScheme, kk::Ring)
  return CartierDivisor(X, kk, IdDict{EffectiveCartierDivisor, elem_type(kk)}())
end

function *(a::RingElem, C::EffectiveCartierDivisor)
  return CartierDivisor(ambient_scheme(C), parent(a), IdDict{EffectiveCartierDivisor, typeof(a)}([C => a]))
end
function *(a::Integer, C::EffectiveCartierDivisor)
  return CartierDivisor(ambient_scheme(C), ZZ, IdDict{EffectiveCartierDivisor, elem_type(ZZ)}([C => ZZ(a)]))
end

function ==(C::CartierDivisor, D::CartierDivisor)
  return iszero(C-D)
end

@doc raw"""
    effective_cartier_divisor(I::IdealSheaf;
                              trivializing_covering::Covering = default_covering(scheme(I)))
                                                                -> EffectiveCartierDivisor

Return the effective Cartier divisor defined by the ideal sheaf `I`, given
that `I` is principal in the given covering of the scheme on which it is defined.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P);

julia> II = IdealSheaf(Y, I);

julia> effective_cartier_divisor(II)
Effective cartier divisor
  on scheme over QQ covered with 3 patches
    1: [(y//x), (z//x)]   affine 2-space
    2: [(x//y), (z//y)]   affine 2-space
    3: [(x//z), (y//z)]   affine 2-space
defined by
  sheaf of ideals with restrictions
    1: Ideal (-(y//x)^2*(z//x) + 1)
    2: Ideal ((x//y)^3 - (z//y))
    3: Ideal ((x//z)^3 - (y//z)^2)
```
"""
effective_cartier_divisor(I::AbsIdealSheaf; trivializing_covering::Covering = default_covering(scheme(I)), check::Bool = true) = EffectiveCartierDivisor(I, trivializing_covering=trivializing_covering, check=check)

@doc raw"""
    effective_cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
    
Return the effective Cartier divisor on the projective scheme ``X`` defined by the homogeneous 
polynomial ``f``. 
"""
function effective_cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
  parent(f) === homogeneous_coordinate_ring(IP) || error("element does not belong to the correct ring")
  d = degree(f)
  X = covered_scheme(IP)
  triv_dict = IdDict{AbsAffineScheme, RingElem}()
  for U in affine_charts(X)
    triv_dict[U] = dehomogenization_map(IP, U)(f)
  end
  C = EffectiveCartierDivisor(X, triv_dict, trivializing_covering=default_covering(X))
  return C
end

@doc raw"""
    cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
    
Return the (effective) Cartier divisor on the projective scheme ``X`` defined by the homogeneous 
polynomial ``f``. 
"""
function cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
  return one(ZZ)*effective_cartier_divisor(IP, f)
end

### Decomposition of an effective Cartier Divisor into irreducible components
### (specialized variant of associated_points, using pure codimension 1
###  and taking multiplicities into account)
@doc raw"""
    irreducible_decomposition(C::EffectiveCartierDivisor; check::Bool=true) -> AbsWeilDivisor

Return $C$ as a linear combination of prime Weil divisors.

Assumes that the ambient scheme is integral and locally noetherian.
"""
function irreducible_decomposition(C::EffectiveCartierDivisor; check::Bool=true)
  X = ambient_scheme(C)
  @check is_integral(X) "ambient scheme must be integral"
  cov = default_covering(X)
  OOX = OO(X)

  charts_todo = copy(patches(cov))
  I = ideal_sheaf(C)
  associated_primes_temp = Vector{Tuple{typeof(I), Int}}()  ## already identified components

  # run through all charts and collect further irreducible components
  while length(charts_todo) > 0
    U = pop!(charts_todo)
    !is_one(I(U)) || continue                                ## supp(C) might not meet all charts
    I_temp=I(U)

    for (J,_) in associated_primes_temp
      !is_one(J(U)) || continue
      I_temp=saturation(I_temp,J(U))                         ## kick out known components
      !is_one(I_temp) || break                               ## break if nothing left
    end

    !is_one(I_temp) || break                                 ## break if nothing left
    components_here = minimal_primes(I_temp)
    for comp in components_here
      I_temp, saturation_index = saturation_with_index(I_temp, comp)
      I_sheaf_temp = PrimeIdealSheafFromChart(X, U, comp, check=false)
      push!(associated_primes_temp, (I_sheaf_temp, saturation_index))
    end
  end
  D = WeilDivisor(X, ZZ, IdDict(associated_primes_temp); check=false)
  return D
end

### Conversion into WeilDivisors
function weil_divisor(C::EffectiveCartierDivisor;
    is_prime::Bool=false # Indicate whether this divisor is already prime
  )
  return WeilDivisor(ideal_sheaf(C), ZZ, check=is_prime)
end

function weil_divisor(C::CartierDivisor)
  X = ambient_scheme(C)
  kk = coefficient_ring(C)
  result = WeilDivisor(X, kk)
  for c in components(C)
    result = result + C[c]*weil_divisor(c)
  end
  return result
end

@doc raw"""
    intersect(W::AbsWeilDivisor, C::EffectiveCartierDivisor; check::Bool=true)
    
Computes the intersection of ``W`` and ``C`` as in [Ful98](@cite) and 
returns an `AbsAlgebraicCycle` of codimension ``2``.
"""
function intersect(W::AbsWeilDivisor, C::EffectiveCartierDivisor; check::Bool=true)
  X = ambient_scheme(W)
  result = zero(W)
  for I in components(irreducible_decomposition(W))
    inc_Y = CoveredClosedEmbedding(X, I, check=false)
    Y = domain(inc_Y)
    pbC = pullback(inc_Y)(C) # Will complain if the defining equation of C is vanishing identically on Y
    W_sub = weil_divisor(pbC)
    result = result + W[I] * pushforward(inc_Y)(W_sub)
  end
  return result
end

@doc raw"""
    intersect(W::AbsWeilDivisor, C::CartierDivisor; check::Bool=true)

Computes the intersection of ``W`` and ``C`` as in [Ful98](@cite) and 
returns an `AbsAlgebraicCycle` of codimension ``2``.
"""
function intersect(W::AbsWeilDivisor, C::CartierDivisor; check::Bool=true)
  result = zero(W)
  for c in components(C)
    result = result + C[c] * intersect(W, c, check=check)
  end
  return result
end

function intersect(D::EffectiveCartierDivisor, C::EffectiveCartierDivisor)
  return intersect(irreducible_decomposition(weil_divisor(D)), C)
end

function intersect(D::EffectiveCartierDivisor, C::CartierDivisor)
  return intersect(irreducible_decomposition(weil_divisor(D)), C)
end

function intersect(D::CartierDivisor, C::EffectiveCartierDivisor)
  return intersect(irreducible_decomposition(weil_divisor(D)), C)
end

function intersect(D::CartierDivisor, C::CartierDivisor)
  return intersect(irreducible_decomposition(weil_divisor(D)), C)
end

dim(C::EffectiveCartierDivisor) = dim(ambient_scheme(C))-1
dim(C::CartierDivisor) = dim(ambient_scheme(C))-1

###########################################################################
## show functions for Cartier divisors
########################################################################### 
function Base.show(io::IO, C::EffectiveCartierDivisor)
  io = pretty(io)
  if get(io, :show_semi_compact, false)
    cov = Oscar._covering_for_printing(io, ambient_scheme(C))
    n = get(io, :label, "")
    _show_semi_compact(io, C, cov, n)
  elseif is_terse(io)
    print(io, "Effective cartier divisor")
  elseif has_attribute(C, :name)
    print(io, get_attribute(C, :name))
  else
    print(io, "Effective cartier divisor on ", Lowercase())
    show(io, ambient_scheme(C))
  end
end

# We keep track of the covering, so that we have more flexibility and
# consistency
function Base.show(io::IO, ::MIME"text/plain", C::EffectiveCartierDivisor)
  io = pretty(io)
  I = ideal_sheaf(C)
  X = ambient_scheme(C)
  cov = Oscar._covering_for_printing(io, X)

  print(io, "Effective cartier divisor")
  if has_attribute(C, :name)
    print(io, " ", get_attribute(C, :name))
  end
  println(io)
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => cov), X)
  println(io, Dedent())
  println(io, "defined by", Lowercase())
  print(io, Indent())
  show(IOContext(io, :show_semi_compact => true), I)
  print(io, Dedent())
end

# Use for nested printings: we omit the ambient variety, but we keep track of
# covering used in the nested printing, and we use `cov`
#
# For nested printings in morphisms, we need to distinguish labels from charts
# of the domain and of the codomain, to pass it to the description of ideal
# sheaves - this is done via the string `n`
#
# We usually use "a" for the domain and "b" for the codomain
function _show_semi_compact(io::IO, C::EffectiveCartierDivisor, cov::Covering, n::String)
  io = pretty(io)
  X = ambient_scheme(C)
  print(io, "Effective cartier divisor")
  if has_attribute(C, :name)
    print(io, " ", get_attribute(C, :name))
  end
  println(io, " defined by")
  print(io, Indent(), Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => cov, :label => n), ideal_sheaf(C))
  print(io, Dedent())
end

function Base.show(io::IO, C::CartierDivisor)
  io = pretty(io)
  if get(io, :show_semi_compact, false)
    cov = Oscar._covering_for_printing(io, ambient_scheme(C))
    n = get(io, :label, "")
    _show_semi_compact(io, C, cov, n)
  elseif is_terse(io)
    print(io, "Cartier divisor")
  elseif has_attribute(C, :name)
    print(io, get_attribute(C, :name))
  else
    print(io, "Cartier divisor on ", Lowercase())
    show(io, ambient_scheme(C))
  end
end

# We keep track of the covering, so that we have more flexibility and
# consistency
#
# We have to take care of some offsets to have our coefficients aligned on the
# right.
function Base.show(io::IO, ::MIME"text/plain", C::CartierDivisor)
  io = pretty(io)
  X = ambient_scheme(C)
  cov = Oscar._covering_for_printing(io, X)
  cc = components(C)
  if length(cc) == 0
    print(io, "Zero cartier divisor ")
    print(io, Indent(), "on ", Lowercase())
    show(IOContext(io, :covering => cov), ambient_scheme(C))
    print(io, Dedent())
  else
    print(io, "Cartier divisor")
    if has_attribute(C, :name)
      print(io, " ", get_attribute(C, :name))
    end
    println(io)
    print(io, Indent(), "on ", Lowercase())
    show(IOContext(io, :covering => cov), ambient_scheme(C))
    println(io)
    println(io, Dedent(), "with coefficients in ", Lowercase(), coefficient_ring(C))
    print(io, "defined by the formal sum of")
    println(io, Indent())
    co_str = String["$(C[I])" for I in cc]
    k = max(length.(co_str)...)
    for i in 1:length(components(C))
      I = cc[i]
      kI = length(co_str[i])
      print(io, " "^(k-kI)*"$(C[I]) * ")
      print(io, Lowercase())
      show(IOContext(io, :show_scheme => false), I)
      #show(IOContext(io, :show_scheme => false), ideal_sheaf(I))
      print(io, "\n")
    end
    print(io, Dedent())
  end
end

# Use for nested printings: we omit the ambient variety, but we keep track of
# covering used in the nested printing, and we use `cov`
#
# For nested printings in morphisms, we need to distinguish labels from charts
# of the domain and of the codomain, to pass it to the description of ideal
# sheaves - this is done via the string `n`
#
# We usually use "a" for the domain and "b" for the codomain
function _show_semi_compact(io::IO, C::CartierDivisor, cov::Covering, n::String)
  io = pretty(io)
  X = ambient_scheme(C)
  cc = components(C)
  if length(cc) == 0
    print(io, "Zero cartier divisor")
  else
    print(io, "Cartier divisor")
    if has_attribute(C, :name)
      print(io, " ", get_attribute(C, :name))
    end
    println(io, " defined by the formal sum of")
    print(io, Indent())
    co_str = String["$(C[I])" for I in cc]
    k = max(length.(co_str)...)
    for i in 1:length(components(C))
      I = cc[i]
      kI = length(co_str[i])
      print(io, " "^(k-kI)*"$(C[I]) * ")
      print(io, Lowercase())
      show(IOContext(io, :show_semi_compact => true, :covering => cov, :label => n), ideal_sheaf(I))
      if i != length(components(C))
        println(io, "--------------------------------------------------------------------------------")
      end
    end
    print(io, Dedent())
  end
end

function Base.hash(X::CartierDivisor, u::UInt)
  return u
end

function Base.hash(X::EffectiveCartierDivisor, u::UInt)
  return u
end

