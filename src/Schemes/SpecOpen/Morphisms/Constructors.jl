export inclusion_morphism, identity_map, restrict, maximal_extension
########################################################################
# Special constructors for morphisms of SpecOpens                      #
########################################################################
### Construct a morphisms by the restriction of coordinate functions 
# to every affine patch
function SpecOpenMor(U::SpecOpen, V::SpecOpen, f::Vector{<:RingElem}; check::Bool=true)
  Y = ambient(V)
  maps = [SpecMor(W, Y, f) for W in affine_patches(U)]
  return SpecOpenMor(U, V, [SpecMor(W, Y, f) for W in affine_patches(U)], check=check) 
end

function SpecOpenMor(X::SpecType, d::RET, Y::SpecType, e::RET, f::Vector{RET}; check::Bool=true) where {SpecType<:Spec, RET<:RingElem}
  U = SpecOpen(X, [d], check=check)
  V = SpecOpen(Y, [e], check=check)
  return SpecOpenMor(U, V, [SpecMor(U[1], Y, OO(U[1]).(f), check=check)], check=check)
end

########################################################################
# Constructors for canonical maps                                      #
########################################################################
function inclusion_morphism(U::SpecOpen, V::SpecOpen; check::Bool=true)
  X = ambient(U)
  if check 
    issubset(U, V) || error("method not implemented")
  end
  return SpecOpenMor(U, V, gens(ambient_ring(X)), check=false)
end

inclusion_morphism(X::SpecOpen, Y::AbsSpec; check::Bool=true) = inclusion_morphism(X, SpecOpen(Y), check=check)

function identity_map(U::SpecOpen) 
  phi = SpecOpenMor(U, U, 
                    [SpecMor(V, ambient(U), gens(OO(V)), check=false) for V in affine_patches(U)], 
                    check=false
                   )
  return phi
end

########################################################################
# Restrictions of morphisms to SpecOpens                               #
########################################################################
@doc Markdown.doc"""
    restrict(f::SchemeMor, U::Scheme, V::Scheme; check::Bool=true)

Return the restriction ``g: U → V`` of ``f`` to ``U`` and ``V``.
"""
restrict(f::SchemeMor, U::Scheme, V::Scheme; check::Bool)

function restrict(f::SpecMor, U::SpecOpen, V::SpecOpen; check::Bool=true)
  if check
    issubset(U, domain(f)) || error("$U is not contained in the domain of $f")
    all(x->issubset(preimage(f, x), U), affine_patches(V)) || error("preimage of $V is not contained in $U")
  end
  return SpecOpenMor(U, V, [restrict(f, W, ambient(V), check=check) for W in affine_patches(U)])
end

function restrict(
    f::SpecOpenMor,
    U::SpecOpen,
    V::SpecOpen;
    check::Bool=true
  )
  if check
    issubset(U, domain(f)) || error("the given open is not an open subset of the domain of the map")
    issubset(V, codomain(f)) || error("the given open is not an open subset of the codomain of the map")
    issubset(preimage(f,V), U) || error("f(U) is not contained in V")
  end
  inc = inclusion_morphism(U, domain(f), check=check)
  help_map = compose(inc, f, check=check)
  Y = ambient(V)
  h = [restrict(g, domain(g), Y, check=check) for g in maps_on_patches(help_map)]
  return SpecOpenMor(U, V, h, check=check)
end

function restrict(f::SpecOpenMor, W::AbsSpec, Y::AbsSpec; check::Bool=true)
  if check
    issubset(W, domain(f)) || error("$U is not contained in the domain of $f")
    issubset(W, preimage(f, Y)) || error("image of $W is not contained in $Y")
  end
  phi = restriction_map(domain(f), W)
  fy = [phi(pullback(f)(y)) for y in OO(codomain(f)).(gens(ambient_ring(Y)))]
  return SpecMor(W, Y, fy, check=false)
end


### the restriction of a morphism to closed subsets in domain and codomain
#function restrict(
#    f::SpecOpenMor,
#    X::AbsSpec,
#    Y::AbsSpec;
#    check::Bool=true
#  )
#  U = intersect(X, domain(f), check=check)
#  V = intersect(Y, codomain(f), check=check)
#
#  new_maps_on_patches = [restrict(f[i], U[i], Y, check=check) for i in 1:npatches(U)]
#
#  return SpecOpenMor(U, V, new_maps_on_patches, check=check)
#end

########################################################################
# Maximal extensions of rational maps given by rational coordinate     #
# functions on the affine patches.                                     #
########################################################################
@Markdown.doc """
    maximal_extension(X::AbsSpec, Y::AbsSpec, f::AbstractAlgebra.Generic.Frac)

Given a rational map ``ϕ : X ---> Y ⊂ Spec 𝕜[y₁,…,yₙ]`` of affine schemes 
determined by ``ϕ*(yⱼ) = fⱼ = aⱼ/bⱼ``, find the maximal open subset ``U⊂ X`` 
to which ``ϕ`` can be extended to a regular map ``g : U → Y`` and return ``g``.
"""
function maximal_extension(
    X::AbsSpec,
    Y::AbsSpec,
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {RET<:RingElem}
  U, g = maximal_extension(X, f)
  n = length(affine_patches(U))
  maps = Vector{AbsSpecMor}()
  for i in 1:n
    push!(maps, SpecMor(affine_patches(U)[i], Y, [restrictions(a)[i] for a in g]))
  end
  h = SpecOpenMor(U, SpecOpen(Y), maps)
  return h
end

function maximal_extension(
    X::AbsSpec, 
    Y::AbsSpec, 
    f::Vector{<:RingElem}
  )
  h = maximal_extension(X, Y, FractionField(ambient_ring(X)).(f))
  return h
end

