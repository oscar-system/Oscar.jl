
########################################################################
# Attributes of PrincipalOpenSubsets                                   #
########################################################################

underlying_scheme(U::PrincipalOpenSubset) = U.U

@doc raw"""
    ambient_scheme(U::PrincipalOpenSubset) -> AffineScheme

For ``U = D(f) ⊆ X`` a principal open subset of ``X`` this returns ``X``.

This is not to be confused with the ambient affine space ``X ⊆ 𝔸 ⁿ``.
"""
ambient_scheme(U::PrincipalOpenSubset) = U.X

function complement_equation(U::PrincipalOpenSubset) 
  if !isdefined(U, :h)
    U.h = prod(OO(U).(U.f); init=one(OO(U)))
  end
  return U.h  # ::elem_type(OO(ambient_scheme(U)))
end

# An internal getter to get factors of polynomial representatives of the complement equation's numerator.
function poly_complement_equations(U::PrincipalOpenSubset)
  if !isdefined(U, :f)
    U.f = lifted_numerator.(complement_equation(U))
  end
  return U.f::Vector{elem_type(ambient_coordinate_ring(U))}
end

### assure compatibility with AffineSchemeOpenSubscheme
complement_equations(U::PrincipalOpenSubset) = [lifted_numerator(complement_equation(U))]
number_of_complement_equations(U::PrincipalOpenSubset) = 1
number_of_generators(U::PrincipalOpenSubset) = 1
gens(U::PrincipalOpenSubset) = [U]
gen(U::PrincipalOpenSubset, i::Int) = (i == 1 ? U : error("index out of range"))
getindex(U::PrincipalOpenSubset, i::Int) = (i == 1 ? U : error("index out of range"))

function inclusion_morphism(U::PrincipalOpenSubset; check::Bool=false)
  if !isdefined(U, :inc)
    X = ambient_scheme(U)
    inc = morphism(U, X, hom(OO(X), OO(U), gens(OO(U)), check=check), check=check)
    U.inc = PrincipalOpenEmbedding(inc, [complement_equation(U)], check=check)
  end
  return U.inc::PrincipalOpenEmbedding
end
