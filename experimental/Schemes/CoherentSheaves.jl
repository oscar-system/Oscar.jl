export CoherentSheaf
export tautological_bundle
export twisting_sheaf
export cotangent_sheaf
export tangent_sheaf
export sheaf_of_rings
export dual
export LineBundle

abstract type AbsCoherentSheaf{
                               SpaceType, OpenType,
                               OutputType, RestrictionType
                              } <: AbsPreSheaf{
                                               SpaceType, OpenType,
                                               OutputType, RestrictionType
                                              } 
end

### Interface for coherent sheaves

@Markdown.doc """
    scheme(F::AbsCoherentSheaf)

Return the scheme on which this sheaf is defined.
"""
scheme(F::AbsCoherentSheaf) = space(underlying_presheaf(F))

@Markdown.doc """
    sheaf_of_rings(F::AbsCoherentSheaf) 

Return the sheaf of rings over which ``ℱ`` is defined.
"""
function sheaf_of_rings(F::AbsCoherentSheaf) 
  error("method not implemented for coherent sheaves of type $(typeof(F))")
end

### The following provides a function for the internal checks that 
# a given set U is open in and admissible for sheaves of modules on X.
#
# We allow the following cases:
#
#  * U::PrincipalOpenSubset in W===ambient_scheme(U) in the basic charts of X
#  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
#  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
#    and U and V contained in the glueing domains of their ambient schemes
#  * U::AbsSpec ⊂ U::AbsSpec in the basic charts of X
#  * U::AbsSpec ⊂ X for U in the basic charts
#  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
function _is_open_for_modules(X::AbsCoveredScheme)
  function is_open_func(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
    C = default_covering(X)
    A = ambient_scheme(U)
    A in C || return false
    B = ambient_scheme(V)
    B in C || return false
    if A === B
      is_subset(U, V) || return false
    else
      G = C[A, B] # Get the glueing
      f, g = glueing_morphisms(G)
      is_subset(U, domain(f)) || return false
      gU = preimage(g, U)
      is_subset(gU, V) || return false
    end
    return true
  end
  function is_open_func(U::PrincipalOpenSubset, Y::AbsCoveredScheme)
    return Y === X && ambient_scheme(U) in default_covering(X)
  end
  function is_open_func(U::AbsSpec, Y::AbsCoveredScheme)
    return Y === X && U in default_covering(X)
  end
  function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
    return X === Y === Z
  end
  function is_open_func(U::AbsSpec, V::AbsSpec)
    U in default_covering(X) || return false
    V in default_covering(X) || return false
    G = default_covering(X)[U, V]
    return issubset(U, glueing_domains(G)[1])
  end
  function is_open_func(U::PrincipalOpenSubset, V::AbsSpec)
    V in default_covering(X) || return false
    ambient_scheme(U) === V && return true
    W = ambient_scheme(U)
    W in default_covering(X) || return false
    G = default_covering(X)[W, V]
    return is_subset(U, glueing_domains(G)[1])
  end
  return is_open_func
end


########################################################################
# Coherent sheaves of modules on covered schemes                       #
########################################################################
@Markdown.doc """
    SheafOfModules <: AbsPreSheaf

A sheaf of modules ``ℳ`` on an `AbsCoveredScheme` ``X``.

Note that due to technical reasons, the admissible open subsets are restricted
to the following:
 * `U::AbsSpec` among the `basic_patches` of the `default_covering` of `X`;
 * `U::PrincipalOpenSubset` with `ambient_scheme(U)` in the `basic_patches` of the `default_covering` of `X`.

One can call the restriction maps of ``ℳ`` across charts implicitly using the
identifications given by the glueings in the `default_covering`.
"""
@attributes mutable struct SheafOfModules{SpaceType, OpenType, 
                                          OutputType,
                                          RestrictionType
                                         } <: AbsCoherentSheaf{
                                                               SpaceType, OpenType,
                                                               OutputType, RestrictionType
                                                              }
  ID::IdDict{AbsSpec, ModuleFP} # the modules on the basic patches of the default covering
  OOX::StructureSheafOfRings # the structure sheaf on X
  M::PreSheafOnScheme # the underlying presheaf of modules for caching

  ### Sheaves of modules on covered schemes
  function SheafOfModules(X::AbsCoveredScheme, 
      MD::IdDict{AbsSpec, ModuleFP}, # A dictionary of modules on the `affine_charts` of `X`
      MG::IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}; # A dictionary for pairs `(U, V)` of 
                                                       # `affine_charts` of `X` such that 
                                                       # A = MG[(U, V)] has entries aᵢⱼ with 
                                                       # gᵢ = ∑ⱼ aᵢⱼ ⋅ fⱼ on U ∩ V with gᵢ the 
                                                       # restrictions of the generators of M[U]
                                                       # and fⱼ the restrictions of the generators 
                                                       # of MD[V]. The aᵢⱼ are elements of 𝒪(U ∩ V)
                                                       # represented as a subset of V.
      check::Bool=true
    )
    OOX = OO(X)
    #OOX = StructureSheafOfRings(X)

    ### Production of the modules on open sets; to be cached
    function production_func(U::AbsSpec, object_cache::IdDict, restriction_cache::IdDict)
      haskey(MD, U) && return MD[U]
      error("module on $U was not found")
    end
    function production_func(U::PrincipalOpenSubset, object_cache::IdDict, restriction_cache::IdDict)
      V = ambient_scheme(U)
      MV = production_func(V, object_cache, restriction_cache)
      rho = OOX(V, U)
      MU, phi = change_base_ring(rho, MV)
      restriction_cache[(V, U)] = phi # Also cache the restriction map
      return MU
    end

    ### Production of the restriction maps; to be cached
#   function restriction_func(V::AbsSpec, MV::ModuleFP, U::AbsSpec, MU::ModuleFP)
#     # This will only be called in case V and U are `affine_charts` of `X`.
#     # In particular, we may assume that U is contained in the glueing 
#     # overlap U ∩ V of U and V (by virtue of the is_open_func check).
#     rho = OOX(V, U) 
#     A = MG[(V, U)] # The transition matrix
#     return hom(MV, MU, [sum([A[i, j] * MU[j] for j in 1:ngens(MU)]) for j in 1:ngens(MV)], rho)
#   end
#   function restriction_func(V::AbsSpec, MV::ModuleFP, U::PrincipalOpenSubset, MU::ModuleFP)
#     # There are two cases: Either U is a PrincipalOpenSubset of V, or U 
#     # is a PrincipalOpenSubset of another affine_chart W. In both cases, 
#     # we need to find W (which equals V in the first case) and use the transition 
#     # matrix for the changes between these two sets. 
#     W = U
#     while W isa PrincipalOpenSubset
#       W = ambient(W)
#     end
#     rho = OOX(V, U) 
#     A = MG[(V, W)] # The transition matrix
#     return hom(MV, MU, [sum([A[i, j] * MU[j] for j in 1:ngens(MU)]) for j in 1:ngens(MV)], rho)
#   end
    function restriction_func(V::AbsSpec, U::AbsSpec, 
        object_cache::IdDict, restriction_cache::IdDict
      )
      MV = haskey(object_cache, V) ? object_cache[V] : production_func(V, object_cache, restriction_cache)
      MU = haskey(object_cache, U) ? object_cache[U] : production_func(U, object_cache, restriction_cache)
      # There are two cases: Either U is a PrincipalOpenSubset of V, or U 
      # is a PrincipalOpenSubset of another affine_chart W. In both cases, 
      # we need to find W (which equals V in the first case) and use the transition 
      # matrix for the changes between these two sets. 
      VV = V
      while VV isa PrincipalOpenSubset
          VV = ambient_scheme(VV)
      end
      UU = U
      while UU isa PrincipalOpenSubset
          UU = ambient_scheme(UU)
      end
      rho = OOX(V, U) 
      A = MG[(VV, UU)] # The transition matrix
      return hom(MV, MU, [sum([A[i, j] * MU[j] for j in 1:ngens(MU)]) for i in 1:ngens(MV)], rho)
    end

    Mpre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=ModuleFP,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_for_modules(X)
                     )
    M = new{typeof(X), AbsSpec, ModuleFP, Hecke.Map}(MD, OOX, Mpre)
    if false
      # Check that all sheaves of modules are compatible on the overlaps.
      # TODO: eventually replace by a check that on every basic
      # affine patch, the ideal sheaf can be inferred from what is
      # given on one dense open subset.
      C = default_covering(X)
      for U in basic_patches(default_covering(X))
        for V in basic_patches(default_covering(X))
          G = C[U, V]
          A, B = glueing_domains(G)
          for i in 1:ngens(A)
            I(A[i]) == ideal(OOX(A[i]), I(V, A[i]).(gens(I(V)))) || error("ideals do not coincide on overlap")
          end
          for i in 1:ngens(B)
            I(B[i]) == ideal(OOX(B[i]), I(U, B[i]).(gens(I(U)))) || error("ideals do not coincide on overlap")
          end
        end
      end
    end
    return M
  end
end

underlying_presheaf(M::SheafOfModules) = M.M
sheaf_of_rings(M::SheafOfModules) = M.OOX

function twisting_sheaf(IP::ProjectiveScheme{<:Field}, d::Int)
    X = covered_scheme(IP)
    MD = IdDict{AbsSpec, ModuleFP}()
    S = ambient_coordinate_ring(IP)
    n = ngens(S)-1
    for i in 1:n+1
        U = affine_charts(X)[i]
        MD[U] = FreeMod(OO(U), ["$(symbols(S)[i])^$d"])
    end

    MG = IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}()
    C = default_covering(X)
    for G in values(glueings(C))
        (U, V) = patches(G)
        (UU, VV) = glueing_domains(G)
        h_U = complement_equation(UU)
        h_V = complement_equation(VV)
        MG[(U, V)] = (d>= 0 ? OO(VV)(h_V^d) : (inv(OO(VV)(h_V))^(-d)))*one(MatrixSpace(OO(VV), 1, 1))
        MG[(V, U)] = (d>= 0 ? OO(UU)(h_U^d) : (inv(OO(UU)(h_U))^(-d)))*one(MatrixSpace(OO(UU), 1, 1))
    end

    M = SheafOfModules(X, MD, MG)
    return M
end

function tautological_bundle(IP::ProjectiveScheme{<:Field})
    return twisting_sheaf(IP, -1)
end

function cotangent_sheaf(X::AbsCoveredScheme)
    MD = IdDict{AbsSpec, ModuleFP}()
    for U in affine_charts(X)
        MD[U] = cotangent_module(U)
    end
    MG = IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}()
    C = default_covering(X)
    for G in values(glueings(C))
        (U, V) = patches(G)
        (UU, VV) = glueing_domains(G)
        (f, g) = glueing_morphisms(G)
        MG[(U, V)] = transpose(jacobi_matrix(pullback(g).(gens(OO(UU)))))
        MG[(V, U)] = transpose(jacobi_matrix(pullback(f).(gens(OO(VV)))))
    end

    
    M = SheafOfModules(X, MD, MG)
    return M
end

function cotangent_module(X::AbsSpec{<:Field, <:MPolyRing})
    R = OO(X)
    F = FreeMod(R, ["d$(x)" for x in symbols(R)])
    return F
end

@attr SheafOfModules function cotangent_module(X::AbsSpec{<:Field, <:MPolyLocalizedRing})
    R = OO(X)
    P = base_ring(R)
    F = FreeMod(R, ["d$(x)" for x in symbols(P)])
    return F
end

@attr ModuleFP function cotangent_module(X::AbsSpec{<:Field, <:MPolyQuo})
    R = OO(X)
    P = base_ring(R)
    F = FreeMod(R, ["d$(x)" for x in symbols(P)])
    rels, _ = sub(F, transpose(change_base_ring(R, jacobi_matrix(gens(modulus(R))))))
    M, _ = quo(F, rels)
    return M
end

function cotangent_module(X::AbsSpec{<:Field, <:MPolyQuoLocalizedRing})
    R = OO(X)
    P = base_ring(R)
    F = FreeMod(R, ["d$(x)" for x in symbols(P)])
    rels, _ = sub(F, transpose(change_base_ring(R, jacobi_matrix(gens(modulus(R))))))
    M, _ = quo(F, rels)
    return M
end

########################################################################
# Hom-Sheaves                                                          #
########################################################################
#= 
  Hom sheaves ℋom(ℱ, 𝒢) are special. 

  First of all, they can be made completely lazy, as their modules 
  on U ⊂ X can be created from ℱ(U) and 𝒢(U) on the fly and the same
  holds for their transition- and restriction functions.

  Second, Hom sheaves come with a domain, a codomain, and an 
  interpretation mapping of their sections as homomorphisms.

  We realize hom sheaves in this way, taking only ℱ and 𝒢 as input 
  in the constructor.
=#

@attributes mutable struct HomSheaf{SpaceType, OpenType, OutputType,
                                    RestrictionType
                                   } <: AbsCoherentSheaf{
                                                         SpaceType, OpenType,
                                                         OutputType, RestrictionType
                                                        }
  domain::AbsCoherentSheaf{SpaceType, OpenType, OutputType, RestrictionType}
  codomain::AbsCoherentSheaf{SpaceType, OpenType, OutputType, RestrictionType}
  OOX::StructureSheafOfRings
  M::PreSheafOnScheme

  function HomSheaf(F::AbsCoherentSheaf, G::AbsCoherentSheaf)
    X = scheme(F)
    X === scheme(G) || error("sheaves must be defined over the same scheme")
    OOX = sheaf_of_rings(F)
    OOX === sheaf_of_rings(G) || error("sheaves must be defined over the same sheaves of rings")

    ### Production of the modules on open sets; to be cached
    function production_func(U::AbsSpec, object_cache::IdDict, restriction_cache::IdDict)
      return hom(F(U), G(U))[1]
    end

    function restriction_func(V::AbsSpec, U::AbsSpec,
        object_cache::IdDict, restriction_cache::IdDict
      )
      MV = haskey(object_cache, V) ? object_cache[V] : production_func(V, object_cache, restriction_cache)
      MU = haskey(object_cache, U) ? object_cache[U] : production_func(U, object_cache, restriction_cache)
      dom_res = F(V, U)
      cod_res = G(V, U)
      f = gens(F(V))
      rf = dom_res.(f)
      # The following two lines will work, because a set of generators for ℱ(V) 
      # always restricts to a set of generators for ℱ(U). Due to changes of 
      # charts, this might be a non-trivial change of bases, however.
      dom_sub, inc = sub(F(U), rf)
      B = [coordinates(e, dom_sub) for e in ambient_representatives_generators(F(U))]
      images = elem_type(MU)[]
      for phi in gens(MV)
        phi_map = element_to_homomorphism(phi)
        images_f = [sum([B[i][j]*cod_res(phi_map(f[j])) for j in 1:length(f)]) for i in 1:length(B)]
        psi = hom(F(U), G(U), images_f)
        push!(images, homomorphism_to_element(MU, psi))
      end

      return hom(MV, MU, images, OOX(V, U)) # TODO: Set check=false?
    end
      
    Mpre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=ModuleFP,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_for_modules(X)
                     )
    M = new{typeof(X), AbsSpec, ModuleFP, Hecke.Map}(F, G, OOX, Mpre)

    return M
  end
end

underlying_presheaf(M::HomSheaf) = M.M
sheaf_of_rings(M::HomSheaf) = M.OOX
domain(M::HomSheaf) = M.domain
codomain(M::HomSheaf) = M.codomain

function free_module(R::StructureSheafOfRings, n::Int)
  return free_module(R, ["e_$i" for i in 1:n])
end

function free_module(R::StructureSheafOfRings, gen_names::Vector{String})
  return free_module(R, Symbol.(gen_names))
end

function free_module(R::StructureSheafOfRings, gen_names::Vector{Symbol})
  X = space(R)
  n = length(gen_names)
  MD = IdDict{AbsSpec, ModuleFP}()
  for U in affine_charts(X)
    MD[U] = FreeMod(OO(U), gen_names)
  end

  MG = IdDict{Tuple{AbsSpec, AbsSpec}, MatrixElem}()
  C = default_covering(X)
  for G in values(glueings(C))
    (U, V) = patches(G)
    (UU, VV) = glueing_domains(G)
    MG[(U, V)] = one(MatrixSpace(OO(VV), n, n))
    MG[(V, U)] = one(MatrixSpace(OO(UU), n, n))
  end

  M = SheafOfModules(X, MD, MG)
  return M
end

function dual(M::SheafOfModules)
  OOX = sheaf_of_rings(M)
  F = free_module(OOX, ["1"])
  return HomSheaf(M, F)
end

@attr HomSheaf function tangent_sheaf(X::AbsCoveredScheme)
  return dual(cotangent_sheaf(X))
end

@attributes mutable struct LineBundle{SpaceType, OpenType, OutputType,
                                      RestrictionType
                                     } <: AbsCoherentSheaf{
                                                           SpaceType, OpenType,
                                                           OutputType, RestrictionType
                                                          }
  OOX::StructureSheafOfRings # the structure sheaf on X
  M::PreSheafOnScheme # the underlying presheaf of modules for caching
  
  function LineBundle(C::CartierDivisor;
      check::Bool=true
    )
    X = variety(C) 
    OOX = OO(X)

    U = default_covering(X)[1]
    ### Production of the modules on open sets; to be cached
    function production_func(U::AbsSpec, object_cache::IdDict, restriction_cache::IdDict)
        return FreeMod(OO(U), 1)
    end
    
    function restriction_func(V::AbsSpec, U::AbsSpec, 
        object_cache::IdDict, restriction_cache::IdDict
      )
      MV = haskey(object_cache, V) ? object_cache[V] : production_func(V, object_cache, restriction_cache)
      MU = haskey(object_cache, U) ? object_cache[U] : production_func(U, object_cache, restriction_cache)
      # There are two cases: Either U is a PrincipalOpenSubset of V, or U 
      # is a PrincipalOpenSubset of another affine_chart W. In both cases, 
      # we need to find W (which equals V in the first case) and use the transition 
      # matrix for the changes between these two sets. 
      VV = V
      while VV isa PrincipalOpenSubset
          VV = ambient_scheme(VV)
      end
      UU = U
      while UU isa PrincipalOpenSubset
          UU = ambient_scheme(UU)
      end
      rho = OOX(V, U) 
      VVtoV = OOX(VV, V)
      UUtoU = OOX(UU, U)
      cond, c = divides(rho(VVtoV(C(VV))), UUtoU(C(UU)))
      cond || error("invalid transition")
      return hom(MV, MU, [c*MU[1]], rho)
    end

    Mpre = PreSheafOnScheme(X, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=FreeMod,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_for_modules(X)
                     )
    M = new{typeof(X), AbsSpec, ModuleFP, Hecke.Map}(OOX, Mpre)
    return M
  end
end

underlying_presheaf(L::LineBundle) = L.M
sheaf_of_rings(L::LineBundle) = L.OOX

########################################################################
# Pushforwards of sheaves for closed embeddings                        #
########################################################################

#=
# Let f : X ↪ Y be a closed embedding with ideal sheaf ℐ on Y and ℳ 
# a sheaf of modules on X. For an open set U ⊂ Y we have 
# f_* ℳ (U) to be the 𝒪_X(f⁻¹(U))-module ℳ (f⁻¹(U)), but seen as 
# an 𝒪_Y(U)-module via the natural restriction of functions. 
#
# Mathematically, this is almost an implicit operation. But since we 
# do not have natural bi-module structures, we need to set up a new 
# sheaf of modules on Y, together with canonical identifications 
# with the modules on X. 
#
# It is clear that this can and should be made lazy.
#                                                                     =#

@attributes mutable struct PushforwardSheaf{SpaceType, OpenType, OutputType,
                                            RestrictionType
                                           } <: AbsCoherentSheaf{
                                                                 SpaceType, OpenType,
                                                                 OutputType, RestrictionType
                                                                }
  inc::CoveredClosedEmbedding
  OOX::StructureSheafOfRings
  OOY::StructureSheafOfRings
  M::AbsCoherentSheaf
  ident::IdDict{AbsSpec, Union{Hecke.Map, Nothing}} # a dictionary caching the identifications
  F::PreSheafOnScheme

  function PushforwardSheaf(inc::CoveredClosedEmbedding, M::AbsCoherentSheaf)
    X = domain(inc)
    X === scheme(M) || error("sheaf must be defined over the domain of the embedding")
    OOX = sheaf_of_rings(M)
    Y = codomain(inc)
    OOY = OO(Y)

    ### Production of the modules on open sets; to be cached
    function production_func(U::AbsSpec, object_cache::IdDict, restriction_cache::IdDict)
      # In case X was empty, return the zero module and store nothing in the identifications.
      if isempty(X) 
        ident[U] = nothing
        return FreeMod(OOY(U), 0)
      end

      # Check whether U ⊂ Y has a nontrivial preimage in X
      f = maps_with_given_codomain(inc, U) # there should be at most one!
      if iszero(length(f))
        ident[U] = nothing
        return FreeMod(OOY(U), 0)
      end
      ff = first(f)
      UX = domain(ff)
      MU, ident_map = _pushforward(pullback(ff), image_ideal(ff), M(UX))
      ident[U] = ident_map
      return MU
    end

    function restriction_func(V::AbsSpec, U::AbsSpec, 
        object_cache::IdDict, restriction_cache::IdDict
      )
      MYV = haskey(object_cache, V) ? object_cache[V] : production_func(V, object_cache, restriction_cache)
      MYU = haskey(object_cache, U) ? object_cache[U] : production_func(U, object_cache, restriction_cache)
      incV_list = maps_with_given_codomain(inc, V)
      incU_list = maps_with_given_codomain(inc, U)
      # return the zero homomorphism in case one of the two sets has 
      # empty preimage.
      if iszero(length(incV_list)) || iszero(length(incU_list)) 
        return hom(MYV, MYU, elem_type(MYU)[zero(MYU) for i in 1:ngens(MYV)], OOY(V, U))
      end
      incV = first(incV_list)
      incU = first(incU_list)
      res_orig = M(domain(incV), domain(incU))
      img_gens = res_orig.(gens(M(domain(incV))))
      return hom(MYV, MYU, (x->preimage(ident[U], x)).(img_gens), OOY(V, U))
    end
    
    ident = IdDict{AbsSpec, Union{Hecke.Map, Nothing}}()

    Blubber = PreSheafOnScheme(Y, production_func, restriction_func,
                      OpenType=AbsSpec, OutputType=ModuleFP,
                      RestrictionType=Hecke.Map,
                      is_open_func=_is_open_for_modules(Y)
                     )
    MY = new{typeof(Y), AbsSpec, ModuleFP, Hecke.Map}(inc, OOX, OOY, M, ident, Blubber)
    return MY
  end
end

underlying_presheaf(M::PushforwardSheaf) = M.F
sheaf_of_rings(M::PushforwardSheaf) = M.OOY

########################################################################
# pushforward of modules                                               #
########################################################################
#
# It is assumed that f : R → S is a map of rings such that S ≅ R/I and 
# M is an S-module. We transform M into an R-module by adding the 
# necessary relations. The return value is that new module M', together 
# with its identification map M' → M. Note that we can not give the 
# inverse of this map, since there is no well-defined underlying ring 
# homomorphism.
function _pushforward(f::Hecke.Map{<:Ring, <:Ring}, I::Ideal, M::FreeMod)
  R = domain(f)
  S = codomain(f)
  base_ring(I) === R || error("ideal is not defined over the correct ring")
  base_ring(M) === S || error("module is not defined over the correct ring")
  FR = FreeMod(R, M.S) # M.S are the symbols of M
  #FRtoM = hom(FR, M, gens(M), f)
  MR, res = quo(FR, (I*FR)[1])
  ident = hom(MR, M, gens(M), f)
  return MR, ident
end

function _pushforward(f::Hecke.Map{<:Ring, <:Ring}, I::Ideal, M::SubQuo)
  R = domain(f)
  S = codomain(f)
  base_ring(I) === R || error("ideal is not defined over the correct ring")
  base_ring(M) === S || error("module is not defined over the correct ring")
  FS = ambient_free_module(M)
  Mgens = ambient_representatives_generators(M)
  Mrels = relations(M)
  FR, identF = _pushforward(f, I, FS)
  MRgens = [preimage(identF, v) for v in ambient_representatives_generators(M)]
  MRrels = elem_type(FR)[preimage(identF, w) for w in relations(M)]
  MRrels_ext = vcat(MRrels, elem_type(FR)[g*e for g in gens(I) for e in gens(FR)])
  MR = quo(sub(FR, MRgens)[1], sub(FR, MRrels_ext)[1])[1]
  ident = hom(MR, M, gens(M), f)
  return MR, ident
end


      


