################################################
# 1: The Julia types for FTheoryTools
################################################

abstract type AbstractFTheoryModel end

ToricOrNonToricCoveredScheme = Union{ToricCoveredScheme{BRT}, CoveredScheme{BRT}} where {BRT<:Ring}

@attributes mutable struct ClosedSubschemeModel <: AbstractFTheoryModel
  base_space::ToricOrNonToricCoveredScheme
  ambient_space::ToricOrNonToricCoveredScheme
  #@req typeof(base_space) === typeof(ambient_space) "Base and ambient space must be of the same type"
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
  ClosedSubschemeModel(base_space::ToricOrNonToricCoveredScheme, ambient_space::ToricOrNonToricCoveredScheme) = new(base_space, ambient_space)
end


@attributes mutable struct CompleteIntersectionModel <: AbstractFTheoryModel
  base_space::ToricOrNonToricCoveredScheme
  ambient_space::ToricOrNonToricCoveredScheme
  defining_ideal::MPolyIdeal
  # zero section
  CompleteIntersectionModel(base_space::ToricOrNonToricCoveredScheme, ambient_space::ToricOrNonToricCoveredScheme, defining_ideal::MPolyIdeal) = new(base_space, ambient_space, defining_ideal)
end


@attributes mutable struct HypersurfaceModel <: AbstractFTheoryModel
  base_space::ToricOrNonToricCoveredScheme
  ambient_space::ToricOrNonToricCoveredScheme
  hypersurface_equation::MPolyRingElem{QQFieldElem}
  HypersurfaceModel(base_space::ToricOrNonToricCoveredScheme, ambient_space::ToricOrNonToricCoveredScheme, hypersurface_equation::MPolyRingElem{QQFieldElem}) = new(base_space, ambient_space, hypersurface_equation)
end


@attributes mutable struct GlobalWeierstrassModel
  weierstrass_f::MPolyRingElem{QQFieldElem}
  weierstrass_g::MPolyRingElem{QQFieldElem}
  weierstrass_polynomial::MPolyRingElem{QQFieldElem}
  base_space::ToricOrNonToricCoveredScheme
  ambient_space::ToricOrNonToricCoveredScheme
  fiber_ambient_space::ToricOrNonToricCoveredScheme
  function GlobalWeierstrassModel(weierstrass_f::MPolyRingElem{QQFieldElem},
                            weierstrass_g::MPolyRingElem{QQFieldElem},
                            weierstrass_polynomial::MPolyRingElem{QQFieldElem},
                            base_space::ToricOrNonToricCoveredScheme,
                            ambient_space::ToricOrNonToricCoveredScheme)
    return new(weierstrass_f, weierstrass_g, weierstrass_polynomial, base_space, ambient_space, weighted_projective_space(ToricCoveredScheme, [2,3,1]))
  end
end


@attributes mutable struct GlobalTateModel
  tate_a1::MPolyRingElem{QQFieldElem}
  tate_a2::MPolyRingElem{QQFieldElem}
  tate_a3::MPolyRingElem{QQFieldElem}
  tate_a4::MPolyRingElem{QQFieldElem}
  tate_a6::MPolyRingElem{QQFieldElem}
  tate_polynomial::MPolyRingElem{QQFieldElem}
  base_space::ToricOrNonToricCoveredScheme
  ambient_space::ToricOrNonToricCoveredScheme
  fiber_ambient_space::ToricOrNonToricCoveredScheme
  function GlobalTateModel(tate_a1::MPolyRingElem{QQFieldElem},
                          tate_a2::MPolyRingElem{QQFieldElem},
                          tate_a3::MPolyRingElem{QQFieldElem},
                          tate_a4::MPolyRingElem{QQFieldElem},
                          tate_a6::MPolyRingElem{QQFieldElem},
                          tate_polynomial::MPolyRingElem{QQFieldElem},
                          base_space::ToricOrNonToricCoveredScheme,
                          ambient_space::ToricOrNonToricCoveredScheme)
    return new(tate_a1, tate_a2, tate_a3, tate_a4, tate_a6, tate_polynomial, base_space, ambient_space, weighted_projective_space(ToricCoveredScheme, [2,3,1]))
  end
end


################################################
# 2: conceptual parents
################################################

@attr ClosedSubschemeModel function conceptual_parent(cim::CompleteIntersectionModel)
  # do something
end

@attr CompleteIntersectionModel function conceptual_parent(hm::HypersurfaceModel)
  # do something
end

@attr HypersurfaceModel function conceptual_parent(gtm::GlobalTateModel)
  # do something
end

@attr HypersurfaceModel function conceptual_parent(gwm::GlobalWeierstrassModel)
  # do something
end
