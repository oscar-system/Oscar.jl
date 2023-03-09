
########################################################################
# Attributes of PrincipalOpenSubsets                                   #
########################################################################

underlying_scheme(U::PrincipalOpenSubset) = U.U

@doc Markdown.doc"""
    ambient_scheme(U::PrincipalOpenSubset) -> Spec

For ``U = D(f) ⊆ X`` a principal open subset of ``X`` this returns ``X``. 

This is not to be confused with the ambient affine space ``X ⊆ 𝔸 ⁿ``.
"""
ambient_scheme(U::PrincipalOpenSubset) = U.X
complement_equation(U::PrincipalOpenSubset) = U.f::elem_type(OO(ambient_scheme(U)))

### assure compatibility with SpecOpen 
gens(U::PrincipalOpenSubset) = [lifted_numerator(complement_equation(U))]
ngens(U::PrincipalOpenSubset) = 1
getindex(U::PrincipalOpenSubset, i::Int) = (i == 1 ? U : error("index out of range"))

function inclusion_morphism(U::PrincipalOpenSubset; check::Bool=false) 
  if !isdefined(U, :inc)
    X = ambient_scheme(U)
    inc = SpecMor(U, X, hom(OO(X), OO(U), gens(OO(U)), check=check))
    U.inc = OpenInclusion(inc, ideal(OO(X), complement_equation(U)), check=check)
  end
  return U.inc
end
