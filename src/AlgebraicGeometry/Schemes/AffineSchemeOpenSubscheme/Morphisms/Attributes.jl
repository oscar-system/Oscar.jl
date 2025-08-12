########################################################################
# Type getters                                                         #
########################################################################
morphism_type(U::S, V::T) where {S<:AffineSchemeOpenSubscheme, T<:AffineSchemeOpenSubscheme} = AffineSchemeOpenSubschemeMor{S, T}
morphism_type(::Type{DomainType}, ::Type{CodomainType}) where {DomainType<:AffineSchemeOpenSubscheme, CodomainType<:AffineSchemeOpenSubscheme} = AffineSchemeOpenSubschemeMor{DomainType, CodomainType}

########################################################################
# Basic getters                                                        #
########################################################################
domain(f::AffineSchemeOpenSubschemeMor) = f.domain
codomain(f::AffineSchemeOpenSubschemeMor) = f.codomain
maps_on_patches(f::AffineSchemeOpenSubschemeMor) = f.maps_on_patches
getindex(f::AffineSchemeOpenSubschemeMor, i::Int) = maps_on_patches(f)[i]

function getindex(f::AffineSchemeOpenSubschemeMor, D::AbsAffineScheme)
  U = affine_patches(domain(f))
  D in U || error("affine patch not found in the domain")
  for i = 1:length(U)
    D === U[i] && return f[i]
  end
end

function pullback(f::AffineSchemeOpenSubschemeMor)
  if !isdefined(f, :pullback)
    U = codomain(f)
    V = domain(f)
    # First check whether we are in an easy case where we only have a restriction
    if ambient_coordinate_ring(domain(f)) === ambient_coordinate_ring(codomain(f)) && complement_equations(V) == complement_equations(U)
      function my_restr(a::AffineSchemeOpenSubschemeRingElem)
        return AffineSchemeOpenSubschemeRingElem(OO(V), [OO(V[i])(a[i]) for i in 1:ngens(V)])
      end
      f.pullback = MapFromFunc(OO(U), OO(V), my_restr)
      return f.pullback::Map{typeof(OO(codomain(f))), typeof(OO(domain(f)))}
    end

    pbs_from_ambient = [pullback(g) for g in maps_on_patches(f)]
    d = [pullback(f[i]).(complement_equations(U)) for i in 1:ngens(V)]
    W = [AffineSchemeOpenSubscheme(V[i], ideal(OO(V[i]), d[i]), check=false) for i in 1:ngens(V)]
    # Not every element of d needs to be non-zero. But zero elements will be 
    # discarded by the constructor for the AffineSchemeOpenSubscheme! So we need to manually 
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
    pb_res = [[pullback(restrict(f[i], W[i][j], U[a[i][j]], check=false)) for j in 1:ngens(W[i])] for i in 1:ngens(V)]
    lift_maps = [restriction_map(W[i], V[i], one(ambient_coordinate_ring(V[i])), check=false) for i in 1:ngens(V)]
    function mymap(a::AffineSchemeOpenSubschemeRingElem)
      b = [lift_maps[i](
              AffineSchemeOpenSubschemeRingElem(
                  OO(W[i]), 
                  [pb_res[i][j](a[j]) for j in 1:ngens(U)],
                  check=false)
             ) for i in 1:ngens(V)
          ]
      return AffineSchemeOpenSubschemeRingElem(OO(V), b, check=false)

      # shortcut for regular elements on the ambient variety
      # Does not work because of parent check done by the Map
      # function mymap(a::RingElem)
      #   parent(a) === OO(ambient(U)) || return mymap(OO(ambient(U))(a))
      #   return AffineSchemeOpenSubschemeRingElem(OO(V), 
      #                           [f[i](a) for i in 1:ngens(V)], 
      #                           check=false)
      # end
    end
    f.pullback = MapFromFunc(OO(U), OO(V), mymap)
  end
  return f.pullback::Map{typeof(OO(codomain(f))), typeof(OO(domain(f)))}
end

