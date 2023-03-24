export ==
export ProjectiveScheme
export ProjectiveSchemeMor
export affine_cone
export affine_patch_type
export affine_patch_type
export ambient_coordinate_ring
export base_ring
export base_ring_type
export base_scheme
export base_scheme_type
export codomain
export covered_projection_to_base
export covered_scheme
export covered_scheme_morphism
export dehomogenize
export domain
export relative_ambient_dimension
export fiber_product
export frac_to_homog_pair
export gens
export getindex
export homog_to_frac
export homogeneous_coordinates
export identity_map
export images_of_variables
export inclusion_morphism
export is_well_defined
export map_on_affine_cones
export morphism_type
export poly_to_homog
export projection_to_base
export projective_scheme_type
export projective_space
export set_base_scheme!
export subscheme

@doc Markdown.doc"""
    homog_to_frac(X::ProjectiveScheme) 

Returns a map that converts a polynomial in the 
`ambient_coordinate_ring` of `X` into a function on the
`affine_cone` of `X`.
"""
function homog_to_frac(P::AbsProjectiveScheme)
  return homog_to_frac(underlying_scheme(P))
end

@doc Markdown.doc"""
    poly_to_homog(X::ProjectiveScheme)

Return a map that converts an element of the `ambient_coordinate_ring` 
of the `affine_cone` `CX` of `X` to an element of the 
`ambient_coordinate_ring` of `X`.
"""
function poly_to_homog(P::AbsProjectiveScheme)
  return poly_to_homog(underlying_scheme(P))
end

@doc Markdown.doc"""
    frac_to_homog_pair(X::ProjectiveScheme)

Return a map that converts an element ``f = p/q`` of the ring of 
functions `OO` of the `affine_cone` of `X` into a pair 
``(a, b)`` of elements of the `ambient_coordinate_ring` of `X`
corresponding to ``p`` and ``q``, respectively.
"""
function frac_to_homog_pair(P::AbsProjectiveScheme)
  return frac_to_homog_pair(underlying_scheme(P))
end

########################################################################
# Methods for ProjectiveScheme                                         #
########################################################################


function homog_to_frac(X::ProjectiveScheme) 
  if !has_attribute(X, :homog_to_frac)
    affine_cone(X)
  end
  return get_attribute(X, :homog_to_frac)
end

function poly_to_homog(X::ProjectiveScheme)
  if !has_attribute(X, :poly_to_homog)
    affine_cone(X)
  end
  return get_attribute(X, :poly_to_homog)
end

function frac_to_homog_pair(X::ProjectiveScheme)
  if !has_attribute(X, :frac_to_homog_pair)
    affine_cone(X)
  end
  return get_attribute(X, :frac_to_homog_pair)
end

    
### This is a temporary fix that needs to be addressed in AbstractAlgebra, issue #1105
Generic.ordering(S::MPolyDecRing) = :degrevlex

function (f::MPolyAnyMap{<:MPolyRing, <:AbstractAlgebra.NCRing})(I::MPolyIdeal)
  return ideal(codomain(f), [f(g) for g in gens(I)])
end

# assure compatibility with generic code for MPolyQuos:
lift(f::MPolyRingElem) = f

########################################################################
# Methods for ProjectiveSchemeMor                                      #
########################################################################

### getters 
domain(phi::ProjectiveSchemeMor) = phi.domain
codomain(phi::ProjectiveSchemeMor) = phi.codomain
pullback(phi::ProjectiveSchemeMor) = phi.pullback
function base_ring_morphism(phi::ProjectiveSchemeMor) 
  if isdefined(phi, :base_ring_morphism)
    return phi.base_ring_morphism
  end
  return identity_map(base_ring(domain(phi)))
end

### additional constructors
function ProjectiveSchemeMor(
    X::AbsProjectiveScheme, 
    Y::AbsProjectiveScheme, 
    a::Vector{<:RingElem}
  )
  base_ring(X) === base_ring(Y) || error("projective schemes must be defined over the same base ring")
  Q = graded_coordinate_ring(X)
  all(x->parent(x)===Q, a) || return ProjectiveSchemeMor(X, Y, Q.(a))
  P = graded_coordinate_ring(Y)
  return ProjectiveSchemeMor(X, Y, hom(P, Q, a))
end

# in case we have honest base schemes, also make the map of schemes available
function base_map(phi::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:MPolyQuoLocRing}})
  if !isdefined(phi, :map_on_base_schemes)
    phi.map_on_base_schemes = SpecMor(base_scheme(domain(phi)), base_scheme(codomain(phi)), coefficient_map(pullback(phi)))
  end
  return phi.map_on_base_schemes::SchemeMor
end

function map_on_affine_cones(
    phi::ProjectiveSchemeMor{
                             <:AbsProjectiveScheme{<:Union{MPolyRing, MPolyQuoRing, 
                                                           MPolyQuoLocRing, MPolyLocRing
                                                          }},
                             <:AbsProjectiveScheme{<:Union{MPolyRing, MPolyQuoRing, 
                                                           MPolyQuoLocRing, MPolyLocRing
                                                          }}
                            };
    check::Bool=true
  )
  @show "trigger 2"
  @show typeof(phi)<:ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Field}, <:AbsProjectiveScheme{<:Field}}
  @show typeof(phi)
  if !isdefined(phi, :map_on_affine_cones)
    pb_phi = pullback(phi)
    C_dom, flat_dom = affine_cone(domain(phi))
    C_cod, flat_cod = affine_cone(codomain(phi))
    @show gens(OO(C_cod))
    @show inverse(flat_cod).(gens(OO(C_cod)))
    @show pb_phi
    @show typeof(pb_phi)
    @show pb_phi.(inverse(flat_cod).(gens(OO(C_cod))))
    @show flat_dom.(pb_phi.(inverse(flat_cod).(gens(OO(C_cod)))))
    pb_res = hom(OO(C_cod), OO(C_dom), flat_dom.(pb_phi.(inverse(flat_cod).(gens(OO(C_cod))))), check=check) # TODO: Set check=false
    phi.map_on_affine_cones = SpecMor(C_dom, C_cod, pb_res)
  end
  return phi.map_on_affine_cones

  if !isdefined(phi, :map_on_affine_cones)
    # The diagram is 
    #   P  →  Q
    #   ↓     ↓
    #   X  →  Y
    # corresponding to a map of rings 
    #
    #   T  ←  S
    #   ↑     ↑
    #   B  ←  A
    #
    #
    X = base_scheme(domain(phi))
    B = OO(X)
    Y = base_scheme(codomain(phi))
    A = OO(Y)
    S = ambient_coordinate_ring(codomain(phi))
    T = ambient_coordinate_ring(domain(phi))
    P = domain(phi)
    Q = codomain(phi)
    pb_P = pullback(projection_to_base(P))
    pb_Q = pullback(projection_to_base(Q))
    imgs_base = pb_P.(base_ring_morphism(phi).(gens(A)))
    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), vcat(imgs_fiber, imgs_base), check=check)
  end
  return phi.map_on_affine_cones::AbsSpecMor
end

function map_on_affine_cones(
    phi::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Field}, <:AbsProjectiveScheme{<:Field}};
    check::Bool=true
  )
  @show "trigger"
  if !isdefined(phi, :map_on_affine_cones)
    pb_phi = pullback(phi)
    C_dom, flat_dom = affine_cone(domain(phi))
    C_cod, flat_cod = affine_cone(codomain(phi))
    pb_res = hom(OO(C_cod), OO(C_dom), flat_dom.(pb_phi.(gens(graded_coordinate_ring(codomain(phi))))), check=false)
    phi.map_on_affine_cones = SpecMor(C_dom, C_cod, pb_res)
  end
  return phi.map_on_affine_cones::AbsSpecMor
end

### This method is for the case of a morphism of 
# projective schemes over a common and unchanged base ring/field.
# So it will be the most usual case. 
#function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:AbstractAlgebra.Ring}})
#  if !isdefined(phi, :map_on_affine_cones)
#    S = ambient_coordinate_ring(codomain(phi))
#    T = ambient_coordinate_ring(domain(phi))
#    P = domain(phi)
#    Q = codomain(phi)
#    imgs_fiber = [homog_to_frac(P)(g) for g in pullback(phi).(gens(S))]
#    phi.map_on_affine_cones = SpecMor(affine_cone(P), affine_cone(Q), imgs_fiber)
#  end
#  return phi.map_on_affine_cones::AbsSpecMor
#end
#
function map_on_affine_cones(phi::ProjectiveSchemeMor{<:ProjectiveScheme{<:SpecOpenRing}, <:ProjectiveScheme{<:SpecOpenRing}})
  if !isdefined(phi, :map_on_affine_cones)
    X = domain(phi)
    CX, map_X = affine_cone(X)
    P = ambient_scheme(CX)
    BX = base_scheme(X)
    BP = ambient_scheme(BX)
    Y = codomain(phi)
    CY, map_Y = affine_cone(Y)
    Q = ambient_scheme(CY)
    BY = base_scheme(Y)
    BQ = ambient_scheme(BY)
    @show codomain(map_X) === OO(CX)
    @show parent(map_X(pullback(phi)(first(gens(ambient_coordinate_ring(Y)))))) === OO(CX)
    @show parent(map_X(pullback(phi)(first(gens(ambient_coordinate_ring(Y))))))
    @show OO(CX)
    fiber_coord_imgs = map_X.(pullback(phi).(gens(ambient_coordinate_ring(Y)))) # elements in OO(CX)
    #@show pullback(phi).(pullback(projection_to_base(Y)).(gens(OO(BY))))
    base_coord_imgs = map_X.(pullback(phi).(ambient_coordinate_ring(Y).(gens(OO(BY)))))
    coord_imgs = vcat(base_coord_imgs, fiber_coord_imgs)

    list = [restriction_map(CX, CX[i]).(coord_imgs) for i in 1:ngens(CX)]
    list = [SpecMor(CX[i], Q, restriction_map(CX, CX[i]).(coord_imgs)) for i in 1:ngens(CX)]
    phi.map_on_affine_cones = SpecOpenMor(CX, CY, list, check=false)
  end
  return phi.map_on_affine_cones::SpecOpenMor
end

function is_well_defined(phi::ProjectiveSchemeMor) 
  CP, _ = affine_cone(domain(phi))
  CQ, _ = affine_cone(codomain(phi))
  return issubset(CP, preimage(map_on_affine_cones(phi), CQ))
end

function compose(f::ProjectiveSchemeMor, g::ProjectiveSchemeMor)
  return ProjectiveSchemeMor(domain(f), codomain(g), compose(pullback(g), pullback(f)))
end

function ==(f::ProjectiveSchemeMor, g::ProjectiveSchemeMor) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for s in gens(ambient_coordinate_ring(codomain(f)))
    pullback(f)(s) - pullback(g)(s) in defining_ideal(domain(f)) || return false
  end
  return true
end

function ==(f::ProjectiveSchemeMor{<:ProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}}, 
            g::ProjectiveSchemeMor{<:ProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}}) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return map_on_affine_cones(f) == map_on_affine_cones(g)
end

### additional constructors

function fiber_product(f::Hecke.Map{DomType, CodType}, P::ProjectiveScheme{DomType}) where {DomType<:Ring, CodType<:Ring}
  R = base_ring(P) 
  R == domain(f) || error("rings not compatible")
  Rnew = codomain(f)
  S = ambient_coordinate_ring(P)
  Qambient = projective_space(Rnew, symbols(S))
  Snew = ambient_coordinate_ring(Qambient)
  phi = hom(S, Snew, f, gens(Snew))
  Q = subscheme(Qambient, phi(defining_ideal(P)))
  return Q, ProjectiveSchemeMor(Q, P, hom(S, ambient_coordinate_ring(Q), f, gens(ambient_coordinate_ring(Q))))
end

function fiber_product(
    f::AbsSpecMor, 
    P::ProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
  )
  codomain(f) == base_scheme(P) || error("codomain and base_scheme are incompatible")
  X = domain(f)
  Y = codomain(f)
  Q_ambient = projective_space(X, symbols(ambient_coordinate_ring(P)))
  help_map = hom(ambient_coordinate_ring(P),
                 ambient_coordinate_ring(Q_ambient),
                 pullback(f),
                 gens(ambient_coordinate_ring(Q_ambient))
                )
  I = help_map(defining_ideal(P))
  Q = subscheme(Q_ambient, I)
  return Q, ProjectiveSchemeMor(Q, P, 
                                hom(ambient_coordinate_ring(P),
                                    ambient_coordinate_ring(Q),
                                    pullback(f),
                                    gens(ambient_coordinate_ring(Q))
                                   )
                               )
end

fiber_product(X::AbsSpec, 
              P::ProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
             ) = fiber_product(inclusion_morphism(X, base_scheme(P)), P)

### canonical map constructors

@doc Markdown.doc"""
    inclusion_morphism(P::T, Q::T)

Assuming that ``P ⊂ Q`` is a subscheme, both proper over an inclusion of 
their base schemes, this returns the associated `ProjectiveSchemeMor`.
"""
function inclusion_morphism(
    P::AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}},
    Q::AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}
  )
  X = base_scheme(P)
  Y = base_scheme(Q)
  f = inclusion_morphism(X, Y) # will throw if X and Y are not compatible
  return ProjectiveSchemeMor(P, Q, 
                             hom(ambient_coordinate_ring(Q),
                                 ambient_coordinate_ring(P),
                                 pullback(f), 
                                 gens(ambient_coordinate_ring(P))
                                )
                            )
end

function inclusion_morphism(P::T, Q::T) where {T<:AbsProjectiveScheme{<:AbstractAlgebra.Ring}}
  A = base_ring(Q)
  B = base_ring(P)
  A === B || error("can not compare schemes for non-equal base rings") # TODO: Extend by check for canonical maps, once they are available
  return ProjectiveSchemeMor(P, Q, 
                             hom(ambient_coordinate_ring(Q),
                                 ambient_coordinate_ring(P),
                                 gens(ambient_coordinate_ring(P))
                                )
                            )
end

identity_map(P::ProjectiveScheme) = ProjectiveSchemeMor(P, P, 
                                                        hom(ambient_coordinate_ring(P),
                                                            ambient_coordinate_ring(P),
                                                            gens(ambient_coordinate_ring(P))
                                                           )
                                                       )

@doc Markdown.doc"""
    covered_scheme(P::ProjectiveScheme)
    
Return a `CoveredScheme` ``X`` isomorphic to `P` with standard affine charts given by dehomogenization. 

Use `dehomogenize(P, U)` with `U` one of the `affine_charts` of ``X`` to 
obtain the dehomogenization map from the `ambient_coordinate_ring` of `P` 
to the `coordinate_ring` of `U`.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> Pcov = covered_scheme(P)
covered scheme with 3 affine patches in its default covering
```
"""
@attr AbsCoveredScheme function covered_scheme(P::ProjectiveScheme)
    C = standard_covering(P) 
    X = CoveredScheme(C)
    return X
end

@attr function covered_projection_to_base(X::ProjectiveScheme{<:Union{<:MPolyQuoLocRing, <:MPolyLocRing, <:MPolyQuoRing, <:MPolyRing}})
  if !has_attribute(X, :covering_projection_to_base) 
    C = standard_covering(X)
  end
  covering_projection = get_attribute(X, :covering_projection_to_base)::CoveringMorphism
  projection = CoveredSchemeMorphism(covered_scheme(X), CoveredScheme(codomain(covering_projection)), covering_projection)
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    U::AbsSpec
  ) where {
    CRT<:Ring
  }
  charts = affine_charts(covered_scheme(X))
  any(x->(x===U), charts) || error("second argument is not an affine chart of the first")
  i = findfirst(k->(charts[k] === U), 1:relative_ambient_dimension(X)+1) - 1
  S = ambient_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  return hom(S, OO(U), s)
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    i::Int
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = ambient_coordinate_ring(X)
  C = standard_covering(X)
  U = C[i+1]
  p = covered_projection_to_base(X)
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  return hom(S, OO(U), pullback(p[U]), s)
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    U::AbsSpec
  ) where {
    CRT<:Union{MPolyQuoLocRing, MPolyLocRing, MPolyRing, MPolyQuoRing}
  }
  return dehomogenize(X, X[U][2]-1)
end

@doc Markdown.doc"""
    homogenize(P::AbsProjectiveScheme, U::AbsSpec)

Given an affine chart ``U ⊂ P`` of an `AbsProjectiveScheme` 
``P``, return a method ``h`` for the homogenization of elements 
``a ∈ 𝒪(U)``. 

That means ``h(a)`` returns a pair ``(p, q)`` representing a fraction 
``p/q ∈ S`` of the `ambient_coordinate_ring` of ``P`` such 
that ``a`` is the dehomogenization of ``p/q``.

**Note:** For the time being, this only works for affine 
charts which are of the standard form ``sᵢ ≠ 0`` for ``sᵢ∈ S``
one of the homogeneous coordinates of ``P``.

**Note:** Since this map returns representatives only, it 
is not a mathematical morphism and, hence, in particular 
not an instance of `Hecke.Map`.

**Note:** Since `fraction_field` relies on some implementation 
of division for the elements, we can not return the fraction 
directly. 
"""
function homogenize(P::AbsProjectiveScheme, U::AbsSpec)
  # TODO: Ideally, one needs to provide this function 
  # only once for every pair (P, U), so we should think of 
  # some internal way for caching. The @attr macro is not 
  # suitable for this, because it is not sensitive for 
  # which U is put in. 

  # Find the chart where a belongs to
  X = covered_scheme(P)
  i = findfirst(V->(U===V), affine_charts(X))
  i === nothing && error("the given affine scheme is not one of the standard affine charts")
  
  # Determine those variables which come from the homogeneous 
  # coordinates
  S = ambient_coordinate_ring(P)
  n = ngens(S)
  R = ambient_coordinate_ring(U)
  x = gens(R)
  s = x[1:n-1]
  x = x[n:end]
  B = base_ring(P)
  y = gens(B)
  t = gens(S)

  w = vcat([1 for j in 1:n-1], [0 for j in n:ngens(R)])
  v = gens(S)
  # prepare a vector of elements on which to evaluate the lifts
  popat!(v, i)
  v = vcat(v, S.(gens(B)))
  function my_dehom(a::RingElem)
    parent(a) === OO(U) || error("element does not belong to the correct ring")
    p = lifted_numerator(a)
    q = lifted_denominator(a)
    deg_p = total_degree(p, w)
    deg_q = total_degree(q, w)
    deg_a = deg_p - deg_q
    ss = S[i] # the homogenization variable
    
    # preliminary lifts, not yet homogenized!
    pp = evaluate(p, v) 
    qq = evaluate(q, v)

    # homogenize numerator and denominator
    pp = sum([c*m*ss^(deg_p - total_degree(m)) for (c, m) in zip(coefficients(pp), monomials(pp))])
    qq = sum([c*m*ss^(deg_q - total_degree(m)) for (c, m) in zip(coefficients(qq), monomials(qq))])

    if deg_a > 0
      return (pp, qq*ss^deg_a)
    elseif deg_a <0
      return (pp * ss^(-deg_a), qq)
    end
    return (pp, qq)
  end
  return my_dehom
end

@doc Markdown.doc"""
    total_degree(f::MPolyRingElem, w::Vector{Int})

Given a multivariate polynomial `f` and a weight vector `w` 
return the total degree of `f` with respect to the weights `w`.
"""
function total_degree(f::MPolyRingElem, w::Vector{Int})
  x = gens(parent(f))
  n = length(x)
  n == length(w) || error("weight vector does not have the correct length")
  vals = [sum([degree(m, j)*w[j] for j in 1:n]) for m in monomials(f)]
  return maximum(vals)
end

function getindex(X::ProjectiveScheme, U::AbsSpec)
  Xcov = covered_scheme(X)
  for C in coverings(Xcov)
    for j in 1:npatches(C)
      if U === C[j] 
        return C, j
      end
    end
  end
  return nothing, 0
end

function dehomogenize(
    X::ProjectiveScheme{CRT}, 
    U::AbsSpec
  ) where {
    CRT<:MPolyQuoLocRing
  }
  # look up U in the coverings of X
  cover_of_U, index_of_U = X[U]
  Xcov = covered_scheme(X)
  S = ambient_coordinate_ring(X)

  s = Vector{elem_type(OO(U))}()
  if cover_of_U === standard_covering(X)
    S = ambient_coordinate_ring(X)
    C = standard_covering(X)
    p = covered_projection_to_base(X)
    s = vcat(gens(OO(U))[1:index_of_U-1], [one(OO(U))], gens(OO(U))[index_of_U:relative_ambient_dimension(X)])
    return hom(S, OO(U), pullback(p[U]), s)
  else
    ref = Xcov[cover_of_U, standard_covering(X)]
    V = codomain(ref[U])
    index_of_V = standard_covering(X)[V]
    t = vcat(gens(OO(V))[1:index_of_V-1], [one(OO(V))], gens(OO(V))[index_of_V:relative_ambient_dimension(X)])
    s = pullback(ref[U]).(t)
    pb = compose(ref[U], covered_projection_to_base(X)[V])
    return hom(S, OO(U), pullback(pb), s)
  end
end

function dehomogenize(
    X::ProjectiveScheme{CRT},
    i::Int
  ) where {
    CRT<:AbstractAlgebra.Ring
  }
  i in 0:relative_ambient_dimension(X) || error("the given integer is not in the admissible range")
  S = ambient_coordinate_ring(X)
  C = default_covering(covered_scheme(X))
  U = C[i+1]
  s = vcat(gens(OO(U))[1:i], [one(OO(U))], gens(OO(U))[i+1:relative_ambient_dimension(X)])
  return hom(S, OO(U), s)
end

### Hack for a detour to speed up mapping of elements 
# This is terribly slow in all kinds of quotient rings 
# because of massive checks for `iszero` due to memory 
# management.
function (f::Oscar.MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocRing, <:Nothing})(a::MPolyRingElem)
  if !has_attribute(f, :lifted_map)
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    g = hom(S, L, lift.(f.img_gens))
    set_attribute!(f, :lifted_map, g)
  end
  g = get_attribute(f, :lifted_map)
  return codomain(f)(g(a), check=false)
end

function (f::Oscar.MPolyAnyMap{<:MPolyRing, <:MPolyQuoLocRing, <:MPolyQuoLocalizedRingHom})(a::MPolyRingElem)
  if !has_attribute(f, :lifted_map)
    S = domain(f)
    W = codomain(f)
    L = localized_ring(W)
    g = hom(S, L, x -> lift(f.coeff_map(x)), lift.(f.img_gens))
    set_attribute!(f, :lifted_map, g)
  end
  g = get_attribute(f, :lifted_map)
  return codomain(f)(g(a), check=false)
end

gens(A::MPolyDecRing, i::Int) = A[i]

@doc Markdown.doc"""
    covered_scheme_morphism(f::ProjectiveSchemeMor)

Given a morphism of `ProjectiveScheme`s ``f : X → Y``, construct and 
return the same morphism as a `CoveredSchemeMorphism` of the `covered_scheme`s 
of ``X`` and ``Y``, respectively.
"""
@attr function covered_scheme_morphism(f::ProjectiveSchemeMor)
  PX = domain(f)
  PY = codomain(f)
  SX = ambient_coordinate_ring(PX)
  SY = ambient_coordinate_ring(PY)
  pbf = pullback(f) # The pullback on the free polynomial rings, not the quotients

  X = covered_scheme(PX)
  Y = covered_scheme(PY)

  mor_dict = IdDict{AbsSpec, AbsSpecMor}()
  U = affine_charts(X)
  for i in 1:ngens(SX)
    U_i = U[i]
    dehom = dehomogenize(PX, U_i) # the dehomogenization map SX → 𝒪(Uᵢ)
    for j in 1:ngens(SY)
      y = gens(SY, j)
      denom = dehom(pbf(y))
      V_j = affine_charts(Y)[j]
      U_ij = PrincipalOpenSubset(U_i, denom)
      u = inv(OO(U_ij)(denom))
      mor_dict[U_ij] = SpecMor(U_ij, V_j, 
                               hom(OO(V_j), OO(U_ij), 
                                   [OO(U_ij)(dehom(pbf(gens(SY, k))))*u for k in 1:ngens(SY) if k != j]
                                  )
                              )
    end
  end
  # We skip the glueings for the time being.
  # Eventually, they should be made lazy.
  CC = Covering(collect(keys(mor_dict)), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}())
  inherit_glueings!(CC, default_covering(X))
  phi = CoveringMorphism(CC, default_covering(Y), mor_dict)
  push!(coverings(X), CC)

  ff = CoveredSchemeMorphism(X, Y, phi)
  return ff
end

