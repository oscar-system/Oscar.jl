########################################################################
# Constructors from ideals                                             #
########################################################################

# TODO: Write one dummy constructor for the documentation with an ideal.
@Markdown.doc """
    function SpecOpen(X::AbsSpec, I::MPolyLocalizedIdeal)

Return the complement of the zero locus of ``I`` in ``X``.
"""
function SpecOpen(X::AbsSpec, I::MPolyLocalizedIdeal; check::Bool=true)
  base_ring(I) === OO(X) || error("Ideal does not belong to the correct ring")
  g = [numerator(a) for a in gens(I) if !iszero(numerator(a))]
  return SpecOpen(X, g, check=check)
end

function SpecOpen(X::AbsSpec, I::MPolyQuoLocalizedIdeal; check::Bool=true)
  base_ring(I) === OO(X) || error("Ideal does not belong to the correct ring")
  g = [lifted_numerator(a) for a in gens(I) if !iszero(numerator(a))]
  return SpecOpen(X, g, check=check)
end

function SpecOpen(X::AbsSpec, I::MPolyIdeal; check::Bool=true)
  return SpecOpen(X, [g for g in gens(I) if !iszero(OO(X)(g))], check=check)
end

########################################################################
# Constructors from closed subvarieties                                #
########################################################################
@doc Markdown.doc"""
  complement(X::Scheme, Y::Scheme) -> Scheme

Return the complement ``X \ Y`` of ``Y`` in ``X``.

Since we want the complement `U = X \ Y` to have a well defined scheme structure, 
we require that `Y` is closed in `X`.
"""
complement(X::Scheme,Y::Scheme)

function complement(X::AbsSpec, Z::AbsSpec{<:Ring, <:MPolyRing})
  ambient_ring(X) == ambient_ring(Z) || error("X and Z do not compare")
  return EmptyScheme(base_ring(X))
end

function complement(X::AbsSpec, Z::AbsSpec{<:Ring, <:MPolyQuo})
  ambient_ring(X) == ambient_ring(Z) || error("X and Z do not compare")
  return SpecOpen(X, modulus(OO(Z)))
end

function complement(X::AbsSpec, 
    Z::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing};
    check::Bool=true
  )
  check && (is_closed_embedding(Z, X) || error("not a closed embedding"))
  return SpecOpen(Y, modulus(quotient_ring(OO(Z))))
end

########################################################################
# Conversion from AbsSpec                                              #
########################################################################
SpecOpen(X::AbsSpec) = SpecOpen(X, [one(ambient_ring(X))], check=false)


########################################################################
# Additional constructors                                              #
########################################################################
function product(U::SpecOpen, Y::AbsSpec)
  X = ambient(U)
  P, pX, pY = product(X, Y)
  V = SpecOpen(P, lifted_numerator.(pullback(pX).(gens(U))))
  res_pX = restrict(pX, V, U, check=false)
  res_pY = restrict(pY, V, SpecOpen(Y), check=false)
  return V, res_pX, res_pY
end
  
function subscheme(U::SpecOpen, I::Ideal)
  Z = subscheme(ambient(U), I) #Takes care of coercion and complains if necessary
  return SpecOpen(Z, [g for g in gens(U) if !iszero(OO(Z)(g))])
end

