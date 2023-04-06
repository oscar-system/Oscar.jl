
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

# in case we have honest base schemes, also make the map of schemes available
function base_map(phi::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:MPolyQuoLocRing}})
  if !isdefined(phi, :map_on_base_schemes)
    phi.map_on_base_schemes = SpecMor(base_scheme(domain(phi)), base_scheme(codomain(phi)), coefficient_map(pullback(phi)))
  end
  return phi.map_on_base_schemes::SchemeMor
end

# Map on affine cones over the same base scheme.
function map_on_affine_cones(
    phi::ProjectiveSchemeMor{
                             <:AbsProjectiveScheme{<:Union{MPolyRing, MPolyQuoRing, 
                                                           MPolyQuoLocRing, MPolyLocRing
                                                          }},
                             <:AbsProjectiveScheme{<:Union{MPolyRing, MPolyQuoRing, 
                                                           MPolyQuoLocRing, MPolyLocRing
                                                          }},
                             <:Hecke.Map,
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
    phi.map_on_affine_cones = SpecMor(C_dom, C_cod, pb_res)
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
    phi.map_on_affine_cones = SpecMor(CP, CQ, vcat(imgs_fiber, imgs_base), check=check)
  end
  return phi.map_on_affine_cones::AbsSpecMor
end

function map_on_affine_cones(
    phi::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Field}, <:AbsProjectiveScheme{<:Field},
                            <:Hecke.Map, Nothing};
    check::Bool=true
  )
  if !isdefined(phi, :map_on_affine_cones)
    pb_phi = pullback(phi)
    C_dom, flat_dom = affine_cone(domain(phi))
    C_cod, flat_cod = affine_cone(codomain(phi))
    pb_res = hom(OO(C_cod), OO(C_dom), flat_dom.(pb_phi.(gens(homogeneous_coordinate_ring(codomain(phi))))), check=false)
    phi.map_on_affine_cones = SpecMor(C_dom, C_cod, pb_res)
  end
  return phi.map_on_affine_cones::AbsSpecMor
end

function map_on_affine_cones(phi::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:SpecOpenRing}, <:AbsProjectiveScheme{<:SpecOpenRing}})
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
    list = [SpecMor(CX[i], Q, restriction_map(CX, CX[i]).(coord_imgs)) for i in 1:ngens(CX)]
    phi.map_on_affine_cones = SpecOpenMor(CX, CY, list, check=false)
  end
  return phi.map_on_affine_cones::SpecOpenMor
end

