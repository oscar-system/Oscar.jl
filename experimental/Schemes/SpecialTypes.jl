export PrincipalOpenSubset
export ambient_scheme, complement_equation, inclusion_morphism

export underlying_morphism, complement_ideal, complement_scheme

export image_ideal

export ideal_type

########################################################################
# Methods for PrincipalOpenSubset                                      #
########################################################################
underlying_scheme(U::PrincipalOpenSubset) = U.U
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

PrincipalOpenSubset(X::AbsSpec) = PrincipalOpenSubset(X, one(OO(X)))

function preimage(f::AbsSpecMor, U::PrincipalOpenSubset; check::Bool=true) 
  if ambient_scheme(U) != codomain(f) 
    Z = preimage(f, ambient_scheme(U), check=check)
    return PrincipalOpenSubset(Z, OO(Z)(pullback(f)(complement_equation(U)), check=false))
  end
  return PrincipalOpenSubset(domain(f), pullback(f)(complement_equation(U)))
end

@Markdown.doc """
    generic_fraction(a::MPolyLocalizedRingElem, U::PrincipalOpenSubset)

Given a regular function ``a ∈ 𝒪(U)`` on a principal open 
subset ``U ⊂ X`` of an affine scheme ``X``, return a 
fraction ``p/q`` in `Quot(P)` (where ``P`` is the `ambient_ring` of 
the `ambient` scheme ``X`` of ``U``) which represents ``a``
in the sense that the maximal extension of its restriction 
to ``U`` returns ``a``.
"""
function generic_fraction(a::MPolyLocalizedRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

function generic_fraction(a::MPolyQuoLocalizedRingElem, U::PrincipalOpenSubset)
  X = ambient_scheme(U)
  parent(a) == OO(U) || error("domains are not compatible")
  return lifted_numerator(a)//lifted_denominator(a)
end

@attr function is_dense(U::PrincipalOpenSubset)
  return !is_zero_divisor(complement_equation(U))
end

function is_constant(a::MPolyLocalizedRingElem) 
  reduce_fraction(a)
  return is_constant(numerator(a)) && is_constant(denominator(a))
end

### Already implemented in AA -- but probably buggy?
#function is_zero_divisor(f::MPolyElem)
#  iszero(f) && return true
#  if is_constant(f)
#    c = coefficients(f)[1]
#    return is_zero_divisor(c)
#  end
#  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
#end

function is_zero_divisor(f::MPolyLocalizedRingElem)
  iszero(f) && return true
  if is_constant(f)
    c = first(coefficients(numerator(f)))
    return is_zero_divisor(c)
  end
  return is_zero_divisor(numerator(f))
end

### The following method is only implemented when the coefficient ring is a field.
# The code should be valid generically, but the Singular backend needed for the 
# ideal quotient is probably buggy for non-fields.
function is_zero_divisor(f::MPolyQuoElem{<:MPolyElem{<:FieldElem}})
  iszero(f) && return true
  b = simplify(f)
  # The next block is basically useless when the coefficient ring is 
  # a field, because it is merely another `is_zero`-check. However, 
  # once more functionality is working, it will actually do stuff and 
  # the above signature can be widened.
  if is_constant(lift(b))
    c = first(coefficients(lift(b)))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end

function is_zero_divisor(f::MPolyQuoLocalizedRingElem{<:Field})
  iszero(f) && return true
  # The next block is basically useless when the coefficient ring is 
  # a field, because it is merely another `is_zero`-check. However, 
  # once more functionality is working, it will actually do stuff and 
  # the above signature can be widened.
  if is_constant(lifted_numerator(f)) && is_constant(lifted_denominator(f))
    c = first(coefficients(lift(numerator(f))))
    return is_zero_divisor(c)
  end
  return !is_zero(quotient(ideal(parent(f), zero(f)), ideal(parent(f), f)))
end

########################################################################
# Methods for OpenInclusion                                            #
########################################################################
underlying_morphism(f::OpenInclusion) = f.inc
complement_ideal(f::OpenInclusion) = f.I
complement_scheme(f::OpenInclusion) = f.Z


########################################################################
# Methods for ClosedEmbedding                                          #
########################################################################
underlying_morphism(f::ClosedEmbedding) = f.inc

@Markdown.doc """
    image_ideal(f::ClosedEmbedding)

For a closed embedding ``f : X → Y`` of affine schemes ``X = Spec(S)`` 
into ``Y = Spec(R)`` such that ``S ≅ R/I`` via ``f`` for some ideal 
``I ⊂ R`` this returns ``I``.
"""
image_ideal(f::ClosedEmbedding) = f.I::ideal_type(OO(codomain(f)))

function complement(f::ClosedEmbedding)
  if !isdefined(f, :U)
    U = SpecOpen(codomain(f), image_ideal(f))
    f.U = U
  end
  return f.U
end

ideal_type(::Type{RT}) where {RT<:MPolyRing} = MPolyIdeal{elem_type(RT)}
ideal_type(::Type{RT}) where {PolyType, RT<:MPolyQuo{PolyType}} = MPolyQuoIdeal{PolyType}
ideal_type(R::Ring) = ideal_type(typeof(R))


########################################################################
# Methods for SimpleGlueing                                            #
########################################################################
patches(G::SimpleGlueing) = (G.X, G.Y)
glueing_morphisms(G::SimpleGlueing) = (G.f, G.g)
glueing_domains(G::SimpleGlueing) = (G.U, G.V)

@attr SimpleGlueing function inverse(G::SimpleGlueing)
  Ginv = SimpleGlueing(G.Y, G.X, G.g, G.f, check=false)
  set_attribute!(Ginv, :inverse, G)
  return Ginv
end

function compose(G::GT, H::GT) where {GT<:SimpleGlueing}
  if patches(G)[2] === patches(H)[2]
    return compose(G, inverse(H))
  elseif patches(G)[1] === patches(H)[1]
    return compose(inverse(G), H)
  elseif patches(G)[1] === patches(H)[2]
    return compose(inverse(G), inverse(H))
  end
  X, Y = patches(G)
  Y === patches(H)[1] || error("Glueings not compatible")
  Z = patches(H)[2]
  f, f_inv = glueing_morphisms(G)
  g, g_inv = glueing_morphisms(H)
  U_new = PrincipalOpenSubset(ambient_scheme(domain(f)), 
                              [complement_equation(domain(f)), 
                               lifted_numerator(pullback(f)(complement_equation(domain(g))))
                              ])
  W_new = PrincipalOpenSubset(ambient_scheme(domain(g_inv)), 
                              [complement_equation(domain(g_inv)), 
                               lifted_numerator(pullback(g_inv)(complement_equation(domain(f_inv))))
                              ])
  V_new = PrincipalOpenSubset(ambient_scheme(domain(g)), 
                              [complement_equation(domain(g)), complement_equation(domain(f_inv))]
                             )
  h = compose(restrict(f, U_new, V_new, check=false), 
              restrict(g, V_new, W_new, check=false))
  h_inv = compose(restrict(g_inv, W_new, V_new, check=false),
                  restrict(f_inv, V_new, U_new, check=false))
  set_attribute!(h, :inverse, h_inv)
  set_attribute!(h_inv, :inverse, h)
  return SimpleGlueing(X, Z, h, h_inv)
end

function restrict(G::SimpleGlueing, X::AbsSpec, Y::AbsSpec; check::Bool=true)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  UX = PrincipalOpenSubset(X, OO(X)(lifted_numerator(complement_equation(U))))
  VY = PrincipalOpenSubset(Y, OO(Y)(lifted_numerator(complement_equation(V))))
  f_res = restrict(f, UX, VY, check=check)
  g_res = restrict(g, VY, UX, check=check)
  return SimpleGlueing(X, Y, f_res, g_res)
end

function Glueing(
    X::AbsSpec, Y::AbsSpec, 
    f::AbsSpecMor{<:PrincipalOpenSubset}, 
    g::AbsSpecMor{<:PrincipalOpenSubset};
    check::Bool=true)
  return SimpleGlueing(X, Y, f, g, check=check)
end

### Conversion of a SimpleGlueing to a sophisticated one
function Glueing(G::SimpleGlueing)
  X, Y = patches(G)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  incY = inclusion_morphism(V, Y)
  incX = inclusion_morphism(U, X)
  Uo = SpecOpen(U)
  Vo = SpecOpen(V)
  fo = SpecOpenMor(Uo, Vo, [compose(f, incY)], check=false)
  go = SpecOpenMor(Vo, Uo, [compose(g, incX)], check=false)
  return Glueing(X, Y, fo, go, check=false)
end
