export CartierDivisor
export cartier_divisor
export EffectiveCartierDivisor
export effective_cartier_divisor
export trivializing_covering

@attributes mutable struct EffectiveCartierDivisor{
                                                   CoveredSchemeType<:AbsCoveredScheme
                                                  }
  X::CoveredSchemeType
  I::IdealSheaf
  C::Covering

  function EffectiveCartierDivisor(
      X::AbsCoveredScheme, 
      D::IdDict{<:AbsAffineScheme, <:RingElem};
      trivializing_covering::Covering=begin
        C = Covering(collect(keys(D)), IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}())
        inherit_gluings!(C, default_covering(X))
        C
      end,
      check::Bool=true
    )
    for U in patches(trivializing_covering)
      U in keys(D) || error("the divisor must be prescribed on all patches of its trivializing covering")
    end
    ID = IdDict{AbsAffineScheme, Ideal}()
    for U in keys(D)
      ID[U] = ideal(OO(U), D[U])
    end
    I = IdealSheaf(X, ID, check=check)
    @check begin
      for U in keys(D)
        is_zero_divisor(OO(U)(D[U])) && error("local elements must not be zero divisors")
      end
      # TODO: 
      # - Check that every affine chart is covered
    end
    return new{typeof(X)}(X, I, trivializing_covering)
  end
end

function (C::EffectiveCartierDivisor)(U::AbsAffineScheme)
  return gens(C.I(U))
end

ideal_sheaf(C::EffectiveCartierDivisor) = C.I

scheme(C::EffectiveCartierDivisor) = C.X
trivializing_covering(C::EffectiveCartierDivisor) = C.C

function EffectiveCartierDivisor(I::IdealSheaf; 
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


@attributes mutable struct CartierDivisor{
                                          CoveredSchemeType<:AbsCoveredScheme,
                                          CoeffType<:RingElem
                                         }
  X::CoveredSchemeType
  R::Ring
  coeff_dict::IdDict{EffectiveCartierDivisor, CoeffType}

  function CartierDivisor(X::AbsCoveredScheme, R::Ring, coeff_dict::IdDict{<:EffectiveCartierDivisor, CoeffType}) where {CoeffType<:RingElem}
    all(x->(scheme(x)===X), keys(coeff_dict)) || error("all effective divisors must be defined over the same scheme")
    all(x->(parent(x) === R), values(coeff_dict)) || error("all coefficients must belong to the same parent")
    return new{typeof(X), CoeffType}(X, R, coeff_dict)
  end
end

scheme(C::CartierDivisor) = C.X
coefficient_ring(C::CartierDivisor) = C.R
coefficient_dict(C::CartierDivisor) = C.coeff_dict
getindex(C::CartierDivisor, k::EffectiveCartierDivisor) = coefficient_dict(C)[k]
components(C::CartierDivisor) = collect(keys(coefficient_dict(C)))

function +(C::CartierDivisor, D::CartierDivisor) 
  scheme(C) === scheme(D) || error("divisors must be defined over the same scheme")
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
  return CartierDivisor(scheme(C), coefficient_ring(C), coeff_dict)
end

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
  return CartierDivisor(scheme(C), coefficient_ring(C), coeff_dict)
end

function *(a::Integer, C::CartierDivisor)
  return coefficient_ring(C)(a)*C
end

function -(C::CartierDivisor, D::CartierDivisor) 
  return C + (-one(coefficient_ring(D)))*D
end

function iszero(C::CartierDivisor)
  return iszero(length(keys(coefficient_dict(C)))) || all(k->iszero(C[k]), components(C))
end

@doc raw"""
    cartier_divisor(E::EffectiveCartierDivisor) -> CartierDivisor

Given an effective cartier divisor `E`, return the cartier divisor
$1*E$.

Mathematically both objects are the same, this function is a coercion method
to see effective cartier divisors as irreducible cartier divisor with coefficient
1.

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
  1 * sheaf of ideals
```
"""
cartier_divisor(E::EffectiveCartierDivisor) = CartierDivisor(E)

function CartierDivisor(C::EffectiveCartierDivisor)
  return CartierDivisor(scheme(C), ZZ, IdDict([C => one(ZZ)]))
end

function CartierDivisor(X::AbsCoveredScheme, kk::Ring)
  return CartierDivisor(X, kk, IdDict{EffectiveCartierDivisor, elem_type(kk)}())
end

function *(a::RingElem, C::EffectiveCartierDivisor)
  return CartierDivisor(scheme(C), parent(a), IdDict{EffectiveCartierDivisor, typeof(a)}([C => a]))
end
function *(a::Integer, C::EffectiveCartierDivisor)
  return CartierDivisor(scheme(C), ZZ, IdDict{EffectiveCartierDivisor, elem_type(ZZ)}([C => ZZ(a)]))
end

function ==(C::CartierDivisor, D::CartierDivisor)
  C === D && return true
  for k in components(C)
    iszero(C[k]) || (haskey(coefficient_dict(D), k) && D[k] == C[k]) || error("equality check not implemented in this complicated case")
  end
  for k in components(D) 
    iszero(D[k]) || (haskey(coefficient_dict(C), k) && D[k] == C[k]) || error("equality check not implemented in this complicated case")
  end
  return true
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
effective_cartier_divisor(I::IdealSheaf; trivializing_covering::Covering = default_covering(scheme(I)), check::Bool = true) = EffectiveCartierDivisor(I, trivializing_covering=trivializing_covering, check=check)

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

function cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
  return one(ZZ)*effective_cartier_divisor(IP, f)
end

### Decomposition of an effective Cartier Divisor into irreducible components
### (specialized variant of associated_points, using pure codimension 1
###  and taking multiplicities into account)
@doc raw"""
    irreducible_decomposition(C::EffectiveCartierDivisor)

Return a `Vector` of pairs ``(I,k)`` corresponding to the irreducible components of ``C``. More precisely,  each ``I`` is a prime  `IdealSheaf` corresponding to an irreducible component of ``C`` and ``k``is the multiplicity of this component in ``C``.
"""
function irreducible_decomposition(C::EffectiveCartierDivisor)
  X = scheme(C)
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
      temp_dict=IdDict{AbsAffineScheme,Ideal}()
      temp_dict[U] = comp
      I_sheaf_temp = IdealSheaf(X, extend!(cov, temp_dict), check=false)
      push!(associated_primes_temp, (I_sheaf_temp, saturation_index))
    end
  end
  return(associated_primes_temp)
end

### Conversion into WeilDivisors
function weil_divisor(C::EffectiveCartierDivisor;
    is_prime::Bool=false # Indicate whether this divisor is already prime
  )
  return WeilDivisor(ideal_sheaf(C), ZZ, check=is_prime)

  # TODO: See what we can recycle from the code below.
  X = scheme(C)
  OOX = OO(X)

  decomp = Vector{Tuple{typeof(ideal_sheaf(C)), Int}}()
  if is_prime
    push!(decomp, (ideal_sheaf(C), 1))
  else
    decomp = irreducible_decomposition(C)
  end
  result = WeilDivisor(X, ZZ)

  for (I,k) in decomp
    result = result + k*WeilDivisor(I,ZZ, check=false)
  end

  return result
end

function weil_divisor(C::CartierDivisor)
  X = scheme(C)
  kk = coefficient_ring(C)
  result = WeilDivisor(X, kk)
  for c in components(C)
    result = result + C[c]*weil_divisor(c)
  end
  return result
end

function intersect(W::WeilDivisor, C::EffectiveCartierDivisor; check::Bool=true)
  X = scheme(W)
  result = zero(W)
  for I in components(W)
    @check is_prime(I) "all components of the first argument must be sheaves of prime ideals"
    inc_Y = CoveredClosedEmbedding(X, I, check=false)
    #inc_Y = CoveredClosedEmbedding(X, I, covering=trivializing_covering(C), check=false)
    Y = domain(inc_Y)
    pbC = pullback(inc_Y)(C) # Will complain if the defining equation of C is vanishing identically on Y
    W_sub = weil_divisor(pbC)
    result = result + W[I] * pushforward(inc_Y)(W_sub)
  end
  return result
end

@doc raw"""
    intersect(W::WeilDivisor, C::CartierDivisor; check::Bool=true)

Computes the intersection of ``W`` and ``C`` as in [Ful98](@cite) and 
returns an `AbsAlgebraicCycle` of codimension ``2``.

!!! note
  The `components` of ``W`` must be sheaves of prime ideals; use `irreducible_decomposition(W)` to achieve this. The check for primality can be switched off using `check=false`. 
"""
function intersect(W::WeilDivisor, C::CartierDivisor; check::Bool=true)
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


function pushforward(inc::CoveredClosedEmbedding, W::WeilDivisor)
  X = domain(inc)
  Y = codomain(inc)
  X === scheme(W) || error("divisor not defined on the domain")
  kk = coefficient_ring(W)
  ideal_dict = IdDict{IdealSheaf, elem_type(kk)}()
  for I in components(W)
    pfI = pushforward(inc)(I)
    ideal_dict[pfI] = W[I]
  end
  return WeilDivisor(Y, kk, ideal_dict, check=false)
end

dim(C::EffectiveCartierDivisor) = dim(scheme(C))-1
dim(C::CartierDivisor) = dim(scheme(C))-1

###########################################################################
## show functions for Cartier divisors
########################################################################### 
function Base.show(io::IO, C::EffectiveCartierDivisor)
  io = pretty(io)
  if get(io, :show_semi_compact, false)
    cov = Oscar._covering_for_printing(io, scheme(C))
    n = get(io, :label, "")
    _show_semi_compact(io, C, cov, n)
  elseif get(io, :supercompact, false)
    print(io, "Effective cartier divisor")
  elseif has_attribute(C, :name)
    print(io, get_attribute(C, :name))
  else
    print(io, "Effective cartier divisor on ", Lowercase())
    show(io, scheme(C))
  end
end

# We keep track of the covering, so that we have more flexibility and
# consistency
function Base.show(io::IO, ::MIME"text/plain", C::EffectiveCartierDivisor)
  io = pretty(io)
  I = ideal_sheaf(C)
  X = scheme(C)
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
  X = scheme(C)
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
    cov = Oscar._covering_for_printing(io, scheme(C))
    n = get(io, :label, "")
    _show_semi_compact(io, C, cov, n)
  elseif get(io, :supercompact, false)
    print(io, "Cartier divisor")
  elseif has_attribute(C, :name)
    print(io, get_attribute(C, :name))
  else
    print(io, "Cartier divisor on ", Lowercase())
    show(io, scheme(C))
  end
end

# We keep track of the covering, so that we have more flexibility and
# consistency
#
# We have to take care of some offsets to have our coefficients aligned on the
# right.
function Base.show(io::IO, ::MIME"text/plain", C::CartierDivisor)
  io = pretty(io)
  X = scheme(C)
  cov = Oscar._covering_for_printing(io, X)
  cc = components(C)
  if length(cc) == 0
    print(io, "Zero cartier divisor ")
    print(io, Indent(), "on ", Lowercase())
    show(IOContext(io, :covering => cov), scheme(C))
    print(io, Dedent())
  else
    print(io, "Cartier divisor")
    if has_attribute(C, :name)
      print(io, " ", get_attribute(C, :name))
    end
    println(io)
    print(io, Indent(), "on ", Lowercase())
    show(IOContext(io, :covering => cov), scheme(C))
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
      show(IOContext(io, :show_scheme => false), ideal_sheaf(I))
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
  X = scheme(C)
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
