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
  defining_classes::Dict{String, ToricDivisorClass}
  function HypersurfaceModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                             hypersurface_equation_parametrization::MPolyRingElem,
                             hypersurface_equation::MPolyRingElem,
                             base_space::FTheorySpace,
                             ambient_space::FTheorySpace,
                             fiber_ambient_space::AbsCoveredScheme)
    result = new(explicit_model_sections, hypersurface_equation_parametrization, hypersurface_equation, base_space, ambient_space, fiber_ambient_space)
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end
end


@attributes mutable struct WeierstrassModel <: AbstractFTheoryModel
  explicit_model_sections::Dict{String, <: MPolyRingElem}
  defining_section_parametrization::Dict{String, <: MPolyRingElem}
  weierstrass_ideal_sheaf::AbsIdealSheaf
  weierstrass_polynomial::MPolyRingElem
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  fiber_ambient_space::AbsCoveredScheme
  defining_classes::Dict{String, ToricDivisorClass}

  function WeierstrassModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                            defining_section_parametrization::Dict{String, <: MPolyRingElem},
                            weierstrass_ideal_sheaf::AbsIdealSheaf,
                            base_space::FTheorySpace,
                            ambient_space::FTheorySpace)
    result = new(explicit_model_sections, defining_section_parametrization)
    result.weierstrass_ideal_sheaf = weierstrass_ideal_sheaf
    result.base_space = base_space
    result.ambient_space = ambient_space
    result.fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end

  function WeierstrassModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                          defining_section_parametrization::Dict{String, <: MPolyRingElem},
                          weierstrass_polynomial::MPolyRingElem,
                          base_space::FTheorySpace,
                          ambient_space::FTheorySpace)
    result = new(explicit_model_sections, defining_section_parametrization)
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
  defining_section_parametrization::Dict{String, <: MPolyRingElem}
  tate_ideal_sheaf::AbsIdealSheaf
  tate_polynomial::MPolyRingElem
  base_space::FTheorySpace
  ambient_space::FTheorySpace
  fiber_ambient_space::AbsCoveredScheme
  defining_classes::Dict{String, ToricDivisorClass}

  function GlobalTateModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                          defining_section_parametrization::Dict{String, <: MPolyRingElem},
                          tate_ideal_sheaf::AbsIdealSheaf,
                          base_space::FTheorySpace,
                          ambient_space::FTheorySpace)
    result = new(explicit_model_sections, defining_section_parametrization)
    result.tate_ideal_sheaf = tate_ideal_sheaf
    result.base_space = base_space
    result.ambient_space = ambient_space
    result.fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end

  function GlobalTateModel(explicit_model_sections::Dict{String, <: MPolyRingElem},
                          defining_section_parametrization::Dict{String, <: MPolyRingElem},
                          tate_polynomial::MPolyRingElem,
                          base_space::FTheorySpace,
                          ambient_space::FTheorySpace)
    result = new(explicit_model_sections, defining_section_parametrization)
    result.tate_polynomial = tate_polynomial
    result.base_space = base_space
    result.ambient_space = ambient_space
    result.fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
    result.defining_classes = Dict{String, ToricDivisorClass}()
    return result
  end

end


##############################################################################
# 3 Struct for Quadrillion F-theory Standard Models when read in from database
##############################################################################

struct QSMModel
  # Information about the polytope underlying the F-theory QSM.
  vertices::Vector{Vector{QQFieldElem}}
  poly_index::Int

  # We build toric 3-fold from triangulating the lattice points in said polytope.
  # Oftentimes, there are a lot of such triangulations (up to 10^15 for the case at hand), which we cannot
  # hope to enumerate in a reasonable time in a computer. The following gives us metadata, to gauge how hard
  # this triangulation task is. First, the boolean triang_quick tells if we can hope to enumerate all
  # triangulations in a reasonable time. This in turn is linked to the question if we can find
  # fine regular triangulations of all facets, the difficulty of which scales primarily with the number of
  # lattice points. Hence, we also provide the maximal number of lattice points in a facet of the polytope in question
  # in the integer max_lattice_pts_in_facet. On top of this, an estimate for the total number of triangulations
  # is provided by the big integer estimated_number_oftriangulations. This estimate is exact if triang_quick = true.
  triang_quick::Bool
  max_lattice_pts_in_facet::Int
  estimated_number_of_triangulations::Int

  # We select one of the many triangulations, construct a 3d toric base B3 and thereby the hypersurface model in question,
  # that is then the key object of study of this F-theory construction.
  hs_model::HypersurfaceModel

  # As per usual, topological data of this geometry is important. Key is the triple intersection number of the
  # anticanonical divisor of the 3-dimensional toric base, as well as its Hodge numbers.
  Kbar3::Int
  h11::Int
  h12::Int
  h13::Int
  h22::Int

  # Recall that B3 is 3-dimensional toric variety. Let s in H^0(B3, Kbar_B3), then V(s) is a K3-surface.
  # Moreover, let xi the coordinates of the Cox ring of B3. Then V(xi) is a divisor in B3.
  # Furthermore, Ci = V(xi) cap V(s) is a divisor in the K3-surface V(s). We study these curves Ci in large detail.
  # Here is some information about these curves:
  genus_ci::Dict{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}, Int}
  degree_of_Kbar_of_tv_restricted_to_ci::Dict{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}, Int}
  intersection_number_among_ci_cj::Matrix{Int}
  index_facet_interior_divisors::Vector{Int}
  intersection_number_among_nontrivial_ci_cj::Matrix{Int}
  
  # The collection of the Ci form a nodal curve. To every nodal curve one can associate a (dual) graph. In this graph,
  # every irreducible component of the nodal curve becomes a node/vertex of the dual graph, and every
  # nodal singularity of the nodal curve turns into an edge of the dual graph. Here this is rather simple.
  # Every Ci above is an irreducible component of the nodal curve in question and the topological intersection numbers
  # among the Ci tell us how many nodal singularities link the Ci. Hence, we construct the dual graph as follows:
  # 1. View the Ci as nodes of an undirected graph G.
  # 2. If the top. intersection number of Ci and Cj is zero, there is no edge between the nodes of G corresponding to Ci and Cj.
  # 3. If the top. intersection number of Ci and Cj is n (> 0), then there are n edges between the nodes of G corresponding to Ci and Cj.
  # The following lists the information about this dual graph.
  # Currently, we cannot label the nodes/vertices of a OSCAR graph. However, it is important to remember what vertex/node in the
  # dual graph corresponds to which geometric locus V(xi, s). Therefore, we keep the labels that link the node of the OSCAR graph
  # to the geometric loci V(xi, s) in the vector components_of_dual_graph::Vector{String}. At least for now.
  # Should it ever be possible (favorable?) to directly attach these labels to the graph, one can remove components_of_dual_graph.
  dual_graph::Graph{Undirected}
  components_of_dual_graph::Vector{String}
  degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph::Dict{String, Int64}
  genus_of_components_of_dual_graph::Dict{String, Int64}

  # In our research, we conduct certain combinatoric computations based on this graph, the data in
  # degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph, genus_of_components_of_dual_graph and a bit more meta data
  # that is not currently included (yet). These computations are hard. It turns out, that one can replace the graph with another
  # graph, so that the computations are easier (a.k.a. the runtimes are a lot shorter). In a nutshell, this means to remove
  # a lot of nodes, and adjust the edges accordingly. Let me not go into more details here. A full description can e.g. be found in
  # https://arxiv.org/abs/2104.08297 and the follow-up papers thereof. Here we collect the information of said simplified graph,
  # by mirroring the strategy for the above dual graph.
  simplified_dual_graph::Graph{Undirected}
  components_of_simplified_dual_graph::Vector{String}
  degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph::Dict{String, Int64}
  genus_of_components_of_simplified_dual_graph::Dict{String, Int64}
  
end


################################################
# 4: The julia type for G4-fluxes
################################################

@attributes mutable struct G4Flux
  model::AbstractFTheoryModel
  class::CohomologyClass
  G4Flux(model::AbstractFTheoryModel, class::CohomologyClass) = new(model, class)
end
