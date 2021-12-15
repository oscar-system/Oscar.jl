mutable struct ProjectiveSpace{RingType, RingElemType}
  A::RingType	# the base ring
  n::Int	# the dimension
  S::MPolyRing{RingElemType}

  function ProjectiveSpace(R::RingType, n::Int; var_name::String="s") where {RingType<:Ring}
    return new{RingType, elem_type{R}}(R, n, PolynomialRing(R, [ var_name*"$i" for i in 0:n ]))
  end
end

base_ring(P::ProjectiveSpace) = P.A
fiber_dimension(P::ProjectiveSpace) = P.n
homogeneous_coordinate_ring(P::ProjectiveSpace) = P.S

function change_of_rings(P::ProjectiveSpace{RingType}, phi) where {RingType<:Ring} 
  error("not implemented")
end

mutable struct MorphismOfProjectiveSpaces{
    DomainRingType, DomainRingElemType, 
    CodomainRingType, CodomainRingElemType
  }
  domain::ProjectiveSpace
  codomain::ProjectiveSpace
  phi # the map between the base rings
  var_images # the images of the variables
end
