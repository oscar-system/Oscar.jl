export CartierDivisor
export trivializing_covering

@attributes mutable struct EffectiveCartierDivisor{
                                                   CoveredSchemeType<:AbsCoveredScheme
                                                  }
  X::CoveredSchemeType
  I::IdealSheaf
  C::Covering

  function EffectiveCartierDivisor(
      X::AbsCoveredScheme, 
      D::IdDict{<:AbsSpec, <:RingElem};
      trivializing_covering::Covering=begin
        C = Covering(collect(keys(D)), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}())
        inherit_glueings!(C, default_covering(X))
        C
      end,
      check::Bool=true
    )
    for U in patches(trivializing_covering)
      U in keys(D) || error("the divisor must be prescribed on all patches of its trivializing covering")
    end
    ID = IdDict{AbsSpec, Ideal}()
    for U in keys(D)
      ID[U] = ideal(OO(U), D[U])
    end
    I = IdealSheaf(X, ID, check=check)
    if check
      for U in keys(D)
        is_zero_divisor(OO(U)(D[U])) && error("local elements must not be zero divisors")
      end
      # TODO: 
      # - Check that every affine chart is covered
    end
    return new{typeof(X)}(X, I, trivializing_covering)
  end
end

function (C::EffectiveCartierDivisor)(U::AbsSpec)
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
  eq_dict = IdDict{AbsSpec, RingElem}()
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

  function CartierDivisor(X::AbsCoveredScheme, R::Ring, coeff_dict::IdDict{EffectiveCartierDivisor, CoeffType}) where {CoeffType<:RingElem}
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


function effective_cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
  parent(f) === graded_coordinate_ring(IP) || error("element does not belong to the correct ring")
  d = degree(f)
  X = covered_scheme(IP)
  triv_dict = IdDict{AbsSpec, RingElem}()
  for U in affine_charts(X)
    triv_dict[U] = dehomogenize(IP, U)(f)
  end
  C = EffectiveCartierDivisor(X, triv_dict, trivializing_covering=default_covering(X))
  return C
end

function cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
  return one(ZZ)*effective_cartier_divisor(IP, f)
end

### Conversion into WeilDivisors
function weil_divisor(C::EffectiveCartierDivisor)
  X = scheme(C)
  OOX = OO(X)

  decomp = primary_decomposition(ideal_sheaf(C))

  n = length(decomp)
  primary_components = [a for (a, _) in decomp]
  prime_components = [b for (_, b) in decomp]
  
  result = WeilDivisor(X, ZZ)
  for i in 1:n
    P = prime_components[i]
    Q = primary_components[i]
    k = ZZ(0)
    for U in affine_charts(X)
      isone(P(U)) && continue
      R = base_ring(P(U))
      L, phi = localization(R, complement_of_prime_ideal(P(U)))
      F = FreeMod(L, 1)
      QF, _ = phi(Q(U))*F
      M, _ = quo(F, QF)
      k = length(M)
      break
    end
    result = result + k*WeilDivisor(P, ZZ)
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

function intersect(W::WeilDivisor, C::EffectiveCartierDivisor)
  X = scheme(W)
  result = zero(W)
  for I in components(W)
    inc_Y = CoveredClosedEmbedding(X, I, check=false)
    #inc_Y = CoveredClosedEmbedding(X, I, covering=trivializing_covering(C), check=false)
    Y = domain(inc_Y)
    pbC = pullback(inc_Y)(C) # Will complain if the defining equation of C is vanishing identically on Y
    W_sub = weil_divisor(pbC)
    result = result + W[I] * pushforward(inc_Y)(W_sub)
  end
  return result
end

function intersect(W::WeilDivisor, C::CartierDivisor)
  result = zero(W)
  for c in components(C)
    result = result + C[c] * intersect(W, c)
  end
  return result
end

function intersect(D::EffectiveCartierDivisor, C::EffectiveCartierDivisor)
  return intersect(weil_divisor(D), C)
end

function intersect(D::EffectiveCartierDivisor, C::CartierDivisor)
  return intersect(weil_divisor(D), C)
end

function intersect(D::CartierDivisor, C::EffectiveCartierDivisor)
  return intersect(weil_divisor(D), C)
end

function intersect(D::CartierDivisor, C::CartierDivisor)
  return intersect(weil_divisor(D), C)
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
  I = ideal_sheaf(C)
  X = C.X
  covering = C.C
  n = npatches(covering)

  println(io,"Effective Cartier Divisor on Covered Scheme with ",n," Charts")
end

function show_details(C::EffectiveCartierDivisor)
   show_details(stdout,C)
end

function show_details(io::IO,C::EffectiveCartierDivisor)
  I = ideal_sheaf(C)
  X = C.X

  covering = C.C
  n = npatches(covering)

  println(io,"Effective Cartier Divisor on Covered Scheme with ",n," Charts:\n")

  for (i,U) in enumerate(patches(covering))
    println(io,"Chart $i:")
    println(io,"   $(I(U))")
    println(io," ")
  end
end
