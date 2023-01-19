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

variety(C::EffectiveCartierDivisor) = C.X
scheme(C::EffectiveCartierDivisor) = C.X
trivializing_covering(C::EffectiveCartierDivisor) = C.C

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
gens(C::CartierDivisor) = collect(keys(coefficient_dict(C)))

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
  return iszero(length(keys(coefficient_dict(C)))) || all(k->iszero(C[k]), gens(C))
end

function CartierDivisor(C::EffectiveCartierDivisor)
  return CartierDivisor(scheme(C), ZZ, IdDict([C => one(ZZ)]))
end

function *(a::RingElem, C::EffectiveCartierDivisor)
  return CartierDivisor(scheme(C), parent(a), IdDict{EffectiveCartierDivisor, typeof(a)}([C => a]))
end
function *(a::Integer, C::EffectiveCartierDivisor)
  return CartierDivisor(scheme(C), ZZ(a), IdDict{EffectiveCartierDivisor, elem_type(ZZ)}([C => a]))
end

function ==(C::CartierDivisor, D::CartierDivisor)
  for k in gens(C)
    iszero(C[k]) || (haskey(coefficient_dict(D), k) && D[k] == C[k]) || return false
  end
  for k in gens(D) 
    haskey(coefficient_dict(C), k) && continue
    iszero(D[k]) || return false
  end
  return true
end


function cartier_divisor(IP::AbsProjectiveScheme, f::MPolyElem_dec)
  parent(f) === ambient_coordinate_ring(IP) || error("element does not belong to the correct ring")
  d = total_degree(f)
  X = covered_scheme(IP)
  triv_dict = IdDict{AbsSpec, RingElem}()
  for U in affine_charts(X)
    triv_dict[U] = dehomogenize(IP, U)(f)
  end
  C = EffectiveCartierDivisor(X, triv_dict, trivializing_covering=default_covering(X))
  return one(ZZ)*C
end
