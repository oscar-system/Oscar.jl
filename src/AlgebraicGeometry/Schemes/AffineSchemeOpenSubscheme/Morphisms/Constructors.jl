########################################################################
# Special constructors for morphisms of AffineSchemeOpenSubschemes                      #
########################################################################
### Construct a morphisms by the restriction of coordinate functions
# to every affine patch
function AffineSchemeOpenSubschemeMor(U::AffineSchemeOpenSubscheme, V::AffineSchemeOpenSubscheme, f::Vector{<:RingElem}; check::Bool=true)
  Y = ambient_scheme(V)
  maps = [morphism(W, Y, f, check=check) for W in affine_patches(U)]
  return AffineSchemeOpenSubschemeMor(U, V, [morphism(W, Y, f, check=check) for W in affine_patches(U)], check=check)
end

function AffineSchemeOpenSubschemeMor(X::AffineSchemeType, d::RET, Y::AffineSchemeType, e::RET, f::Vector{RET}; check::Bool=true) where {AffineSchemeType<:AffineScheme, RET<:RingElem}
  U = AffineSchemeOpenSubscheme(X, [d], check=check)
  V = AffineSchemeOpenSubscheme(Y, [e], check=check)
  return AffineSchemeOpenSubschemeMor(U, V, [morphism(U[1], Y, OO(U[1]).(f), check=check)], check=check)
end

########################################################################
# Constructors for canonical maps                                      #
########################################################################
function inclusion_morphism(U::AffineSchemeOpenSubscheme, V::AffineSchemeOpenSubscheme; check::Bool=true)
  X = ambient_scheme(U)
  @check is_subscheme(U, V) "method not implemented"
  return AffineSchemeOpenSubschemeMor(U, V, gens(ambient_coordinate_ring(X)), check=false)
end

inclusion_morphism(X::AffineSchemeOpenSubscheme, Y::AbsAffineScheme; check::Bool=true) = inclusion_morphism(X, AffineSchemeOpenSubscheme(Y), check=check)

function id_hom(U::AffineSchemeOpenSubscheme)
  phi = AffineSchemeOpenSubschemeMor(U, U,
                    [morphism(V, ambient_scheme(U), gens(OO(V)), check=false) for V in affine_patches(U)],
                    check=false
                   )
  return phi
end

########################################################################
# Restrictions of morphisms to AffineSchemeOpenSubschemes                               #
########################################################################
@doc raw"""
    restrict(f::SchemeMor, U::Scheme, V::Scheme; check::Bool=true)

Return the restriction ``g: U → V`` of ``f`` to ``U`` and ``V``.
"""
function restrict(f::SchemeMor, U::Scheme, V::Scheme; check::Bool)
  error("method not implemented")
end

function restrict(f::AffineSchemeMor, U::AffineSchemeOpenSubscheme, V::AffineSchemeOpenSubscheme; check::Bool=true)
  @check begin
    is_subscheme(U, domain(f)) || error("$U is not contained in the domain of $f")
    is_subscheme(V, codomain(f)) || error("$V is not contained in the codomain of $f")
    all(x->is_subscheme(preimage(f, x), U), affine_patches(V)) || error("preimage of $V is not contained in $U")
  end
  return AffineSchemeOpenSubschemeMor(U, V, [restrict(f, W, ambient_scheme(V), check=check) for W in affine_patches(U)])
end

function restrict(
    f::AffineSchemeOpenSubschemeMor,
    U::AffineSchemeOpenSubscheme,
    V::AffineSchemeOpenSubscheme;
    check::Bool=true
  )
  @check begin
    is_subscheme(U, domain(f)) || error("the given open is not an open subset of the domain of the map")
    is_subscheme(V, codomain(f)) || error("the given open is not an open subset of the codomain of the map")
    is_subscheme(preimage(f,V), U) || error("f(U) is not contained in V")
  end
  inc = inclusion_morphism(U, domain(f), check=check)
  help_map = compose(inc, f, check=check)
  Y = ambient_scheme(V)
  h = [restrict(g, domain(g), Y, check=check) for g in maps_on_patches(help_map)]
  return AffineSchemeOpenSubschemeMor(U, V, h, check=check)
end

function restrict(f::AffineSchemeOpenSubschemeMor, W::AbsAffineScheme, Y::AbsAffineScheme; check::Bool=true)
  @check begin
    is_subscheme(W, domain(f)) || error("$U is not contained in the domain of $f")
    is_subscheme(W, preimage(f, Y)) || error("image of $W is not contained in $Y")
  end
  phi = restriction_map(domain(f), W)
  fy = [phi(pullback(f)(y)) for y in OO(codomain(f)).(gens(ambient_coordinate_ring(Y)))]
  return morphism(W, Y, fy, check=false)
end


### the restriction of a morphism to closed subsets in domain and codomain
#function restrict(
#    f::AffineSchemeOpenSubschemeMor,
#    X::AbsAffineScheme,
#    Y::AbsAffineScheme;
#    check::Bool=true
#  )
#  U = intersect(X, domain(f), check=check)
#  V = intersect(Y, codomain(f), check=check)
#
#  new_maps_on_patches = [restrict(f[i], U[i], Y, check=check) for i in 1:n_patches(U)]
#
#  return AffineSchemeOpenSubschemeMor(U, V, new_maps_on_patches, check=check)
#end

########################################################################
# Maximal extensions of rational maps given by rational coordinate     #
# functions on the affine patches.                                     #
########################################################################
@doc raw"""
    maximal_extension(X::AbsAffineScheme, Y::AbsAffineScheme, f::AbstractAlgebra.Generic.FracFieldElem)

Given a rational map ``ϕ : X ---> Y ⊂ Spec 𝕜[y₁,…,yₙ]`` of affine schemes
determined by ``ϕ*(yⱼ) = fⱼ = aⱼ/bⱼ``, find the maximal open subset ``U⊂ X``
to which ``ϕ`` can be extended to a regular map ``g : U → Y`` and return ``g``.
"""
function maximal_extension(
    X::AbsAffineScheme,
    Y::AbsAffineScheme,
    f::Vector{AbstractAlgebra.Generic.FracFieldElem{RET}}
  ) where {RET<:RingElem}
  U, g = maximal_extension(X, f)
  n = length(affine_patches(U))
  maps = Vector{AbsAffineSchemeMor}()
  for i in 1:n
    push!(maps, morphism(affine_patches(U)[i], Y, [restrictions(a)[i] for a in g]))
  end
  h = AffineSchemeOpenSubschemeMor(U, AffineSchemeOpenSubscheme(Y), maps)
  return h
end

function maximal_extension(
    X::AbsAffineScheme,
    Y::AbsAffineScheme,
    f::Vector{<:RingElem}
  )
  h = maximal_extension(X, Y, fraction_field(ambient_coordinate_ring(X)).(f))
  return h
end

