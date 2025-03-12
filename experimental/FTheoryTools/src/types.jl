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
  # mordell_weil_group
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
  defining_classes::Dict{String, ToricDivisorClass}
  model_section_parametrization::Dict{String, <: MPolyRingElem}
  function HypersurfaceModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                             hypersurface_equation_parametrization::MPolyRingElem,
                             hypersurface_equation::MPolyRingElem,
                             base_space::FTheorySpace,
                             ambient_space::FTheorySpace,
                             fiber_ambient_space::AbsCoveredScheme)
    result = new(explicit_model_sections, hypersurface_equation_parametrization, hypersurface_equation, base_space, ambient_space, fiber_ambient_space)
    result.defining_classes = Dict{String, ToricDivisorClass}()
    result.model_section_parametrization = Dict{String, MPolyRingElem}()
    return result
  end
end


@attributes mutable struct WeierstrassModel <: AbstractFTheoryModel
  explicit_model_sections::Dict{String, <: MPolyRingElem}
  model_section_parametrization::Dict{String, <: MPolyRingElem}
  weierstrass_ideal_sheaf::AbsIdealSheaf
  weierstrass_polynomial::MPolyRingElem
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  fiber_ambient_space::AbsCoveredScheme
  defining_classes::Dict{String, ToricDivisorClass}

  function WeierstrassModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                            model_section_parametrization::Dict{String, <: MPolyRingElem},
                            weierstrass_ideal_sheaf::AbsIdealSheaf,
                            base_space::FTheorySpace,
                            ambient_space::FTheorySpace)
    result = new(explicit_model_sections, model_section_parametrization)
    result.weierstrass_ideal_sheaf = weierstrass_ideal_sheaf
    result.base_space = base_space
    result.ambient_space = ambient_space
    result.fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end

  function WeierstrassModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                          model_section_parametrization::Dict{String, <: MPolyRingElem},
                          weierstrass_polynomial::MPolyRingElem,
                          base_space::FTheorySpace,
                          ambient_space::FTheorySpace)
    result = new(explicit_model_sections, model_section_parametrization)
    result.weierstrass_polynomial = weierstrass_polynomial
    result.base_space = base_space
    result.ambient_space = ambient_space
    result.fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end

end


@attributes mutable struct GlobalTateModel <: AbstractFTheoryModel
  explicit_model_sections::Dict{String, <: MPolyRingElem}
  model_section_parametrization::Dict{String, <: MPolyRingElem}
  tate_ideal_sheaf::AbsIdealSheaf
  tate_polynomial::MPolyRingElem
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  fiber_ambient_space::AbsCoveredScheme
  defining_classes::Dict{String, ToricDivisorClass}

  function GlobalTateModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                          model_section_parametrization::Dict{String, <: MPolyRingElem},
                          tate_ideal_sheaf::AbsIdealSheaf,
                          base_space::FTheorySpace,
                          ambient_space::FTheorySpace)
    result = new(explicit_model_sections, model_section_parametrization)
    result.tate_ideal_sheaf = tate_ideal_sheaf
    result.base_space = base_space
    result.ambient_space = ambient_space
    result.fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end

  function GlobalTateModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                          model_section_parametrization::Dict{String, <: MPolyRingElem},
                          tate_polynomial::MPolyRingElem,
                          base_space::FTheorySpace,
                          ambient_space::FTheorySpace)
    result = new(explicit_model_sections, model_section_parametrization)
    result.tate_polynomial = tate_polynomial
    result.base_space = base_space
    result.ambient_space = ambient_space
    result.fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end

end



################################################
# 3: The julia type for G4-fluxes
################################################

@attributes mutable struct G4Flux
  model::AbstractFTheoryModel
  class::CohomologyClass
  G4Flux(model::AbstractFTheoryModel, class::CohomologyClass) = new(model, class)
end



################################################
# 4 The julia type for a family of G4-fluxes
################################################

@attributes mutable struct FamilyOfG4Fluxes
  model::AbstractFTheoryModel
  mat_int::QQMatrix
  mat_rat::QQMatrix
  FamilyOfG4Fluxes(model::AbstractFTheoryModel, mat_int::QQMatrix, mat_rat::QQMatrix) = new(model, mat_int, mat_rat)
end
