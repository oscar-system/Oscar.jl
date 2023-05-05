################################################
# 1: The Julia types for FTheoryTools
################################################

abstract type AbstractFTheoryModel end

space_type = Union{ToricCoveredScheme{BRT}, CoveredScheme{BRT}} where {BRT<:Ring}

@attributes mutable struct ClosedSubschemeModel <: AbstractFTheoryModel
  base_space::space_type
  ambient_space::space_type
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
  ClosedSubschemeModel(base_space::space_type, ambient_space::space_type) = new(base_space, ambient_space)
end


@attributes mutable struct CompleteIntersectionModel <: AbstractFTheoryModel
  base_space::space_type
  ambient_space::space_type
  defining_ideal::MPolyIdeal
  # zero section
  CompleteIntersectionModel(base_space::space_type, ambient_space::space_type, defining_ideal::MPolyIdeal) = new(base_space, ambient_space, defining_ideal)
end


@attributes mutable struct HypersurfaceModel <: AbstractFTheoryModel
  base_space::space_type
  ambient_space::space_type
  hypersurface_equation::MPolyRingElem{QQFieldElem}
  HypersurfaceModel(base_space::space_type, ambient_space::space_type, hypersurface_equation::MPolyRingElem{QQFieldElem}) = new(base_space, ambient_space, hypersurface_equation)
end


@attributes mutable struct GlobalWeierstrassModel
  poly_f::MPolyRingElem{QQFieldElem}
  poly_g::MPolyRingElem{QQFieldElem}
  pw::MPolyRingElem{QQFieldElem}
  toric_base_space::AbstractNormalToricVariety
  toric_ambient_space::AbstractNormalToricVariety
  calabi_yau_hypersurface::ClosedSubvarietyOfToricVariety
  function GlobalWeierstrassModel(poly_f::MPolyRingElem{QQFieldElem},
                            poly_g::MPolyRingElem{QQFieldElem},
                            pw::MPolyRingElem{QQFieldElem},
                            toric_base_space::AbstractNormalToricVariety,
                            toric_ambient_space::AbstractNormalToricVariety,
                            calabi_yau_hypersurface::ClosedSubvarietyOfToricVariety)
    return new(poly_f, poly_g, pw, toric_base_space, toric_ambient_space, calabi_yau_hypersurface)
  end
end


@attributes mutable struct GlobalTateModel
  a1::MPolyRingElem{QQFieldElem}
  a2::MPolyRingElem{QQFieldElem}
  a3::MPolyRingElem{QQFieldElem}
  a4::MPolyRingElem{QQFieldElem}
  a6::MPolyRingElem{QQFieldElem}
  pt::MPolyRingElem{QQFieldElem}
  toric_base_space::AbstractNormalToricVariety
  toric_ambient_space::AbstractNormalToricVariety
  calabi_yau_hypersurface::ClosedSubvarietyOfToricVariety
  function GlobalTateModel(a1::MPolyRingElem{QQFieldElem},
                          a2::MPolyRingElem{QQFieldElem},
                          a3::MPolyRingElem{QQFieldElem},
                          a4::MPolyRingElem{QQFieldElem},
                          a6::MPolyRingElem{QQFieldElem},
                          pt::MPolyRingElem{QQFieldElem},
                          toric_base_space::AbstractNormalToricVariety,
                          toric_ambient_space::AbstractNormalToricVariety,
                          calabi_yau_hypersurface::ClosedSubvarietyOfToricVariety)
    return new(a1, a2, a3, a4, a6, pt, toric_base_space, toric_ambient_space, calabi_yau_hypersurface)
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


################################################
# 3: Constructors
################################################

#=
function global_tate_model(base::AbstractNormalToricVariety)
  toric_ambient_space = _ambient_space_from_base(base)
  (a1, a2, a3, a4, a6) = _tate_sections(base)
  pt = _tate_polynomial([a1, a2, a3, a4, a6], cox_ring(toric_ambient_space))
  calabi_yau_hypersurface = closed_subvariety_of_toric_variety(toric_ambient_space, [pt])
  model = GlobalTateModel(a1, a2, a3, a4, a6, pt, base, toric_ambient_space, calabi_yau_hypersurface)
  set_attribute!(model, :base_fully_specified, true)
  return model
end


################################################
# 3: Display methods
################################################

function Base.show(io::IO, t::GlobalTateModel)
  if base_fully_specified(t)
    print(io, "Global Tate model over a concrete base")
  else
    print(io, "Global Tate model over a not fully specified base")
  end
end
=#
