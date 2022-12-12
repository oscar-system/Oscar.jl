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
    # First check whether we are in an easy case where we only have a restriction
    if ambient_coordinate_ring(domain(f)) === ambient_coordinate_ring(codomain(f)) && gens(V) == gens(U)
      function my_restr(a::SpecOpenRingElem)
        return SpecOpenRingElem(OO(V), [OO(V[i])(a[i]) for i in 1:ngens(V)])
      end
      f.pullback = Hecke.MapFromFunc(my_restr, OO(U), OO(V))
      return f.pullback::Hecke.Map{typeof(OO(codomain(f))), typeof(OO(domain(f)))}
    end

    pbs_from_ambient = [pullback(g) for g in maps_on_patches(f)]
    d = [pullback(f[i]).(gens(U)) for i in 1:ngens(V)]
    W = [SpecOpen(V[i], ideal(OO(V[i]), d[i]), check=false) for i in 1:ngens(V)]
    # Not every element of d needs to be non-zero. But zero elements will be 
    # discarded by the constructor for the SpecOpen! So we need to manually 
    # relate them for the time being. This should only be a temporary fix, however.
    a = Vector{Vector{Int}}()
    for i in 1:length(d)
      b = Vector{Int}()
      j = 1
      for g in d[i]
        if iszero(g) 
          j = j + 1
        else
          push!(b, j)
          j = j + 1
        end
      end
      push!(a, b)
    end
    @show a
    @show ngens(V)
    @show ngens(U)
    pb_res = [[pullback(restrict(f[i], W[i][j], U[a[i][j]], check=false)) for j in 1:ngens(W[i])] for i in 1:ngens(V)]
    lift_maps = [restriction_map(W[i], V[i], one(ambient_coordinate_ring(V[i])), check=false) for i in 1:ngens(V)]
    function mymap(a::SpecOpenRingElem)
      b = [lift_maps[i](
              SpecOpenRingElem(
                  OO(W[i]), 
                  [pb_res[i][j](a[j]) for j in 1:ngens(U)],
                  check=false)
             ) for i in 1:ngens(V)
          ]
      return SpecOpenRingElem(OO(V), b, check=false)

      # shortcut for regular elements on the ambient variety
      # Does not work because of parent check done by the Hecke.Map
      # function mymap(a::RingElem)
      #   parent(a) === OO(ambient(U)) || return mymap(OO(ambient(U))(a))
      #   return SpecOpenRingElem(OO(V), 
      #                           [f[i](a) for i in 1:ngens(V)], 
      #                           check=false)
      # end
    end
    f.pullback = Hecke.MapFromFunc(mymap, OO(U), OO(V))
  end
  return f.pullback::Hecke.Map{typeof(OO(codomain(f))), typeof(OO(domain(f)))}
end

