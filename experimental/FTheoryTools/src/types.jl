################################################
# 1: The julia types for the requires geometries
################################################

@attributes mutable struct FamilyOfSpaces
  coordinate_ring::MPolyRing
  grading::Matrix{Int64}
  dim::Int
  FamilyOfSpaces(coordinate_ring::MPolyRing, grading::Matrix{Int64}, dim::Int) = new(coordinate_ring, grading, dim)
end
const FTheorySpace = Union{AbsCoveredScheme, FamilyOfSpaces}


################################################
# 2: The julia types for the F-Theory models
################################################

abstract type AbstractFTheoryModel end

@attributes mutable struct ClosedSubschemeModel <: AbstractFTheoryModel
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  # total_space
  # fiber_ambient_space
  # singular_loci
  # chow_ring
  # is_cy
  # generic_fiber
  # modell_weil_group
  # weil_chatelet_group
  # gauge_group
  # fiber_diagram (cod. 1,2,3)
  # is_flat
  # Q_factorial_terminal
  # toric_Q_factorial_terminal
  # jacbian_fibration
  # mirror_dual
  ClosedSubschemeModel(base_space::FTheorySpace, ambient_space::FTheorySpace) = new(base_space, ambient_space)
end


@attributes mutable struct CompleteIntersectionModel <: AbstractFTheoryModel
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  defining_ideal::MPolyIdeal
  # zero section
  CompleteIntersectionModel(base_space::FTheorySpace, ambient_space::FTheorySpace, defining_ideal::MPolyIdeal) = new(base_space, ambient_space, defining_ideal)
end


@attributes mutable struct HypersurfaceModel <: AbstractFTheoryModel
  explicit_model_sections::Dict{String, <: MPolyRingElem}
  hypersurface_equation_parametrization::MPolyRingElem
  hypersurface_equation::MPolyRingElem
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  fiber_ambient_space::AbsCoveredScheme
  function HypersurfaceModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                             hypersurface_equation_parametrization::MPolyRingElem,
                             hypersurface_equation::MPolyRingElem,
                             base_space::FTheorySpace,
                             ambient_space::FTheorySpace,
                             fiber_ambient_space::AbsCoveredScheme)
    return new(explicit_model_sections, hypersurface_equation_parametrization, hypersurface_equation, base_space, ambient_space, fiber_ambient_space)
  end
end


@attributes mutable struct WeierstrassModel <: AbstractFTheoryModel
  explicit_model_sections::Dict{String, <: MPolyRingElem}
  defining_section_parametrization::Dict{String, <: MPolyRingElem}
  weierstrass_polynomial::MPolyRingElem
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  fiber_ambient_space::AbsCoveredScheme
  function WeierstrassModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                            defining_section_parametrization::Dict{String, <: MPolyRingElem},
                            weierstrass_polynomial::MPolyRingElem,
                            base_space::FTheorySpace,
                            ambient_space::FTheorySpace)
    fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    return new(explicit_model_sections, defining_section_parametrization, weierstrass_polynomial, base_space, ambient_space, fiber_ambient_space)
  end
end


@attributes mutable struct GlobalTateModel <: AbstractFTheoryModel
  explicit_model_sections::Dict{String, <: MPolyRingElem}
  defining_section_parametrization::Dict{String, <: MPolyRingElem}
  tate_polynomial::MPolyRingElem
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  fiber_ambient_space::AbsCoveredScheme
  function GlobalTateModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                          defining_section_parametrization::Dict{String, <: MPolyRingElem},
                          tate_polynomial::MPolyRingElem,
                          base_space::FTheorySpace,
                          ambient_space::FTheorySpace)
    fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    return new(explicit_model_sections, defining_section_parametrization, tate_polynomial, base_space, ambient_space, fiber_ambient_space)
  end
end
