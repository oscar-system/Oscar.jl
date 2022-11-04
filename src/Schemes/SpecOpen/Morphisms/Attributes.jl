export morphism_type, maps_on_patches, pullback
########################################################################
# Type getters                                                         #
########################################################################
morphism_type(U::S, V::T) where {S<:SpecOpen, T<:SpecOpen} = SpecOpenMor{S, T}
morphism_type(::Type{DomainType}, ::Type{CodomainType}) where {DomainType<:SpecOpen, CodomainType<:SpecOpen} = SpecOpenMor{DomainType, CodomainType}

########################################################################
# Basic getters                                                        #
########################################################################
domain(f::SpecOpenMor) = f.domain
codomain(f::SpecOpenMor) = f.codomain
maps_on_patches(f::SpecOpenMor) = f.maps_on_patches
getindex(f::SpecOpenMor, i::Int) = maps_on_patches(f)[i]

function getindex(f::SpecOpenMor, D::AbsSpec)
  U = affine_patches(domain(f))
  D in U || error("affine patch not found in the domain")
  for i = 1:length(U)
    D === U[i] && return f[i]
  end
end

function pullback(f::SpecOpenMor)
  if !isdefined(f, :pullback)
    U = codomain(f)
    V = domain(f)
    pbs_from_ambient = [pullback(g) for g in maps_on_patches(f)]
    W = [SpecOpen(V[i], ideal(OO(V[i]), pullback(f[i]).(gens(U))), check=false) for i in 1:ngens(V)]
    pb_res = [[pullback(restrict(f[i], W[i][j], U[j], check=false)) for j in 1:ngens(U)] for i in 1:ngens(V)]
    lift_maps = [restriction_map(W[i], V[i], one(ambient_ring(V[i])), check=false) for i in 1:ngens(V)]
    function mymap(a::SpecOpenRingElem)
      b = [lift_maps[i](
              SpecOpenRingElem(
                  OO(W[i]), 
                  [pb_res[i][j](a[j]) for j in 1:ngens(U)],
                  check=false)
             ) for i in 1:ngens(V)
          ]
      return SpecOpenRingElem(OO(V), b, check=false)
    end
    f.pullback = Hecke.MapFromFunc(mymap, OO(U), OO(V))
  end
  return f.pullback::Hecke.Map{typeof(OO(codomain(f))), typeof(OO(domain(f)))}
end

