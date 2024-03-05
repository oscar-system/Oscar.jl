# Deprecated after 0.14.*
Base.@deprecate_binding is_isomorphic_with_symmetric_group is_isomorphic_to_symmetric_group
Base.@deprecate_binding has_is_isomorphic_with_symmetric_group has_is_isomorphic_to_symmetric_group
Base.@deprecate_binding set_is_isomorphic_with_symmetric_group set_is_isomorphic_to_symmetric_group

Base.@deprecate_binding is_isomorphic_with_alternating_group is_isomorphic_to_alternating_group
Base.@deprecate_binding has_is_isomorphic_with_alternating_group has_is_isomorphic_to_alternating_group
Base.@deprecate_binding set_is_isomorphic_with_alternating_group set_is_isomorphic_to_alternating_group

Base.@deprecate_binding proj_space projective_space

Base.@deprecate_binding are_algebraically_independent is_algebraically_independent_with_relations

Base.@deprecate ambient_ring(U::AbsMultSet) ring(U)

# Deprecated after 0.15
Base.@deprecate_binding _compute_glueing_base_change _compute_gluing_base_change
Base.@deprecate_binding _compute_inherited_glueing _compute_inherited_gluing
Base.@deprecate_binding _compute_toric_glueing _compute_toric_gluing
Base.@deprecate_binding add_glueing! add_gluing!
Base.@deprecate_binding base_glueing base_gluing
Base.@deprecate_binding glueing_domains gluing_domains
Base.@deprecate_binding glueing_graph gluing_graph
Base.@deprecate_binding glueing_morphisms gluing_morphisms
Base.@deprecate_binding glueings gluings
Base.@deprecate_binding inherit_glueings! inherit_gluings!
Base.@deprecate_binding is_connected_glueing is_connected_gluing
Base.@deprecate_binding pruned_glueing_graph pruned_gluing_graph
Base.@deprecate_binding underlying_glueing underlying_gluing
Base.@deprecate_binding update_glueing_graph update_gluing_graph
Base.@deprecate_binding AbsGlueing AbsGluing
Base.@deprecate_binding AbsProjectiveGlueing AbsProjectiveGluing
Base.@deprecate_binding BaseChangeGlueingData BaseChangeGluingData
Base.@deprecate_binding CoveredProjectiveGlueingData CoveredProjectiveGluingData
Base.@deprecate_binding Glueing Gluing
Base.@deprecate_binding InheritGlueingData InheritGluingData
Base.@deprecate_binding LazyProjectiveGlueing LazyProjectiveGluing
Base.@deprecate_binding ProjectiveGlueing ProjectiveGluing
Base.@deprecate_binding ProjectiveGlueingData ProjectiveGluingData
Base.@deprecate_binding SimpleGlueing SimpleGluing
Base.@deprecate_binding ToricGlueingData ToricGluingData

Base.@deprecate_binding jacobi_matrix jacobian_matrix
Base.@deprecate_binding jacobi_ideal jacobian_ideal

@deprecate embedding(G::DirectProductGroup, j::Int) canonical_injection(G, j)
@deprecate projection(G::DirectProductGroup, j::Int) canonical_projection(G, j)
@deprecate embedding(G::SemidirectProductGroup, n::Int) canonical_injection(G, n)
@deprecate projection(G::SemidirectProductGroup) canonical_projection(G)
@deprecate embedding(W::WreathProductGroup, n::Int) canonical_injection(W, n)
@deprecate projection(W::WreathProductGroup) canonical_projection(W)

@deprecate num_partitions number_of_partitions
@deprecate num_positive_roots number_of_positive_roots
@deprecate num_roots number_of_roots
@deprecate num_simple_roots number_of_simple_roots
@deprecate num_standard_tableaux number_of_standard_tableaux
@deprecate number_atlas_groups number_of_atlas_groups
@deprecate number_conjugacy_classes number_of_conjugacy_classes
@deprecate has_number_conjugacy_classes has_number_of_conjugacy_classes
@deprecate set_number_conjugacy_classes set_number_of_conjugacy_classes
@deprecate number_moved_points number_of_moved_points
@deprecate has_number_moved_points has_number_of_moved_points
@deprecate set_number_moved_points set_number_of_moved_points
@deprecate number_perfect_groups number_of_perfect_groups
@deprecate has_number_perfect_groups has_number_of_perfect_groups
@deprecate number_primitive_groups number_of_primitive_groups
@deprecate has_number_primitive_groups has_number_of_primitive_groups
@deprecate number_small_groups number_of_small_groups
@deprecate has_number_small_groups has_number_of_small_groups
@deprecate number_transitive_groups number_of_transitive_groups
@deprecate has_number_transitive_groups has_number_of_transitive_groups

@deprecate factorisations factorizations
@deprecate centraliser centralizer

@deprecate hall_subgroup_reps(G::T, P::AbstractVector{<:IntegerUnion}) where T <: Union{GAPGroup, FinGenAbGroup} map(representative, hall_subgroup_classes(G, P))
@deprecate hall_subgroups_representatives(G::GAPGroup, P::AbstractVector{<:IntegerUnion}) map(representative, hall_subgroup_classes(G, P))

function hall_subgroup(G::T, P::AbstractVector{<:IntegerUnion}) where T <: Union{GAPGroup, FinGenAbGroup}
  Base.depwarn("The function hall_subgroup is deprecated. Please use hall_subgroup_classes.", :hall_subgroup)
  @req is_solvable(G) "The group is not solvable"
  return representative(hall_subgroup_classes(G, P)[1])
end

@deprecate low_index_subgroup_reps(G::T, n::Int) where T <: Union{GAPGroup, FinGenAbGroup} map(representative, low_index_subgroup_classes(G, n))

@deprecate complement_class_reps(G::T, N::T) where T <: GAPGroup map(representative, complement_classes(G, N))

@deprecate maximal_subgroup_reps(G::T) where T <: Union{GAPGroup, FinGenAbGroup} map(representative, maximal_subgroup_classes(G))
@deprecate subgroup_reps(G::T) where T <: Union{GAPGroup, FinGenAbGroup} map(representative, subgroup_classes(G))
@deprecate conjugacy_classes_maximal_subgroups(G::T) where T <: Union{GAPGroup, FinGenAbGroup} maximal_subgroup_classes(G)
@deprecate conjugacy_classes_subgroups(G::T) where T <: Union{GAPGroup, FinGenAbGroup} subgroup_classes(G)

@deprecate labelled_matrix_formatted labeled_matrix_formatted

@deprecate Spec AffineScheme
@deprecate proj(E::ToricLineBundle...) projectivization
@deprecate proj(E::ToricDivisor...) projectivization
@deprecate AbsSpec AbsAffineScheme
@deprecate SpecMor AffineSchemeMor
@deprecate AbsSpecMor AbsAffineSchemeMor
@deprecate SimplifiedSpec SimplifiedAffineScheme
@deprecate SpecOpen AffineSchemeOpenSubscheme
@deprecate SpecOpenMor AffineSchemeOpenSubschemeMor
@deprecate SpecOpenRing AffineSchemeOpenSubschemeRing
@deprecate SpecOpenRingElem AffineSchemeOpenSubschemeRingElem
@deprecate SpecSubset AffineSchemeSubset
@deprecate StdSpec StdAffineScheme
