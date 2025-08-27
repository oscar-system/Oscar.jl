### generic getters
domain(f::AbsProjectiveSchemeMorphism) = domain(underlying_morphism(f))
codomain(f::AbsProjectiveSchemeMorphism) = codomain(underlying_morphism(f))
pullback(f::AbsProjectiveSchemeMorphism) = pullback(underlying_morphism(f))
base_ring_morphism(f::AbsProjectiveSchemeMorphism) = base_ring_morphism(underlying_morphism(f))
base_map(f::AbsProjectiveSchemeMorphism) = base_map(underlying_morphism(f))
map_on_affine_cones(f::AbsProjectiveSchemeMorphism) = map_on_affine_cones(underlying_morphism(f))

underlying_morphism(f::AbsProjectiveSchemeMorphism) = error("no underlying morphism for morphisms of type $(typeof(f))")

### getters for the minimal concrete type
domain(phi::ProjectiveSchemeMor) = phi.domain
codomain(phi::ProjectiveSchemeMor) = phi.codomain
@doc raw"""
    pullback(phi::ProjectiveSchemeMor)

For a morphism `phi` of projective schemes, this returns the associated 
morphism of graded affine algebras.
"""
pullback(phi::ProjectiveSchemeMor) = phi.pullback

@doc raw"""
    base_ring_morphism(phi::ProjectiveSchemeMor) 

For a morphism `phi : P → Q` of relative projective spaces 
over `psi : Spec(A) → Spec(B)` this returns the associated 
map `B → A`.
"""
function base_ring_morphism(phi::ProjectiveSchemeMor) 
  if isdefined(phi, :base_ring_morphism)
    return phi.base_ring_morphism
  end
  return id_hom(base_ring(domain(phi)))
end

# in case we have honest base schemes, also make the map of schemes available
@doc raw"""
    base_map(phi::ProjectiveSchemeMor)

For a morphism `phi : P → Q` of relative projective spaces 
over `psi : Spec(A) → Spec(B)` this returns `psi`.
"""
function base_map(phi::ProjectiveSchemeMor)
  if !isdefined(phi, :map_on_base_schemes)
    phi.map_on_base_schemes = morphism(base_scheme(domain(phi)), base_scheme(codomain(phi)), coefficient_map(pullback(phi)))
  end
  return phi.map_on_base_schemes::SchemeMor
end

function base_map(phi::ProjectiveSchemeMor{<:Any, <:Any, <:Any, Nothing})
  error("there exists no associated base map")
end

# Map on affine cones over the same base scheme.
@doc raw"""
    map_on_affine_cones(phi::ProjectiveSchemeMor)

For a morphism `phi : X → Y` this returns the associated morphism 
of the `affine_cone`s ``C(X) → C(Y)``.
"""
function map_on_affine_cones(phi::ProjectiveSchemeMor)
  error("method not implemented for this type of input")
end

function map_on_affine_cones(
    phi::ProjectiveSchemeMor{
                             <:AbsProjectiveScheme{<:Union{MPolyRing, MPolyQuoRing, 
                                                           MPolyQuoLocRing, MPolyLocRing
                                                          }},
                             <:AbsProjectiveScheme{<:Union{MPolyRing, MPolyQuoRing, 
                                                           MPolyQuoLocRing, MPolyLocRing
                                                          }},
                             <:Map,
                             Nothing # This indicates the same base scheme for domain and codomain.
                            };
    check::Bool=true
  )
  if !isdefined(phi, :map_on_affine_cones)
    pb_phi = pullback(phi)
    C_dom, flat_dom = affine_cone(domain(phi))
    C_cod, flat_cod = affine_cone(codomain(phi))
    v = inverse(flat_cod).(gens(OO(C_cod)))
    pb_res = hom(OO(C_cod), OO(C_dom), flat_dom.(pb_phi.(inverse(flat_cod).(gens(OO(C_cod))))), check=check) # TODO: Set check=false
    phi.map_on_affine_cones = morphism(C_dom, C_cod, pb_res)
  end
  return phi.map_on_affine_cones
end

# Map on affine cones with a non-trivial map on base schemes.
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
    CP, mP = affine_cone(P)
    CQ, mQ = affine_cone(Q)
    imgs_fiber = [mP(g) for g in pullback(phi).(gens(S))]
    phi.map_on_affine_cones = morphism(CP, CQ, vcat(imgs_fiber, imgs_base), check=check)
  end
  return phi.map_on_affine_cones::AbsAffineSchemeMor
end

function map_on_affine_cones(
    phi::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Field}, <:AbsProjectiveScheme{<:Field},
                            <:Map, Nothing};
    check::Bool=true
  )
  if !isdefined(phi, :map_on_affine_cones)
    pb_phi = pullback(phi)
    C_dom, flat_dom = affine_cone(domain(phi))
    C_cod, flat_cod = affine_cone(codomain(phi))
    pb_res = hom(OO(C_cod), OO(C_dom), flat_dom.(pb_phi.(gens(homogeneous_coordinate_ring(codomain(phi))))), check=false)
    phi.map_on_affine_cones = morphism(C_dom, C_cod, pb_res)
  end
  return phi.map_on_affine_cones::AbsAffineSchemeMor
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:AffineSchemeOpenSubschemeRing}, <:AbsProjectiveScheme{<:AffineSchemeOpenSubschemeRing}})
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
    fiber_coord_imgs = map_X.(pullback(phi).(gens(ambient_coordinate_ring(Y)))) # elements in OO(CX)
    base_coord_imgs = map_X.(pullback(phi).(ambient_coordinate_ring(Y).(gens(OO(BY)))))
    coord_imgs = vcat(base_coord_imgs, fiber_coord_imgs)

    list = [restriction_map(CX, CX[i]).(coord_imgs) for i in 1:ngens(CX)]
    list = [morphism(CX[i], Q, restriction_map(CX, CX[i]).(coord_imgs)) for i in 1:ngens(CX)]
    phi.map_on_affine_cones = AffineSchemeOpenSubschemeMor(CX, CY, list, check=false)
  end
  return phi.map_on_affine_cones::AffineSchemeOpenSubschemeMor
end


########################################################################
# Special implementations for closed embeddings
########################################################################

underlying_morphism(f::ProjectiveClosedEmbedding) = f.underlying_morphism
image_ideal(f::ProjectiveClosedEmbedding) = f.ideal_of_image 

########################################################################
# Rational maps
########################################################################

### generic getters
domain(f::AbsRationalMap) = domain(underlying_rational_map(f))
codomain(f::AbsRationalMap) = codomain(underlying_rational_map(f))
pullback(f::AbsRationalMap) = pullback(underlying_rational_map(f))
graph_ring(f::AbsRationalMap) = graph_ring(underlying_rational_map(f))

underlying_rational_map(f::AbsRationalMap) = error("no underlying rational map for rational map of type $(typeof(f))")

### getters for the minimal concrete type
domain(phi::RationalMap) = phi.domain
codomain(phi::RationalMap) = phi.codomain
@doc raw"""
    pullback(phi::RationalMap)

For a rational map `phi` of projective varieties, this returns the associated 
morphism of graded affine algebras.
"""
pullback(phi::RationalMap) = phi.pullback

function graph_ring(f::RationalMap)
  if !isdefined(f, :graph_ring)
    pbf = pullback(f)
    SY = domain(pbf)
    SX = codomain(pbf)
    S, inc_SX, inc_SY = tensor_product(SX, SY)
    IG = ideal(S, elem_type(S)[inc_SY(y) - inc_SX(pbf(y)) for y in gens(SY)])
    SG, pr = quo(S, IG)
    f.graph_ring = SG, hom(SX, SG, pr.(inc_SX.(gens(SX))); check=false), hom(SY, SG, pr.(inc_SY.(gens(SY))); check=false)
  end
  return f.graph_ring
end


