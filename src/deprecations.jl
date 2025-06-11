#=
Short explanation of the different ways to deprecate things in julia:
1. Rename a type.
   You only need a single deprecation, but this needs to be a `@deprecate_binding` to be able to use OldType in signatures and type annotations.

Base.@deprecate_binding OldType NewType

2. Rename a function from `old_func` to `new_func` and deprecate `old_func` while the interface stays the same.
   You only need to add one single deprecation that takes care of deprecating and redirecting all methods.

@deprecate old_func new_func

3. Rename a function from `old_func` to `new_func` and deprecate `old_func` while the interface changes.
   For each removed method of `old_func` you need to add a deprecation. As an example, consider a change in the argument order.

@deprecate old_func(a::SomeType, b::SomeOtherType) new_func(b, a)

4. More complex function deprecations.
   Add a method for `old_func` with the old interface that first calls `Base.depwarn` and then calls `new_func` with the new interface.
=#

# Deprecated after 0.14.*
@deprecate is_isomorphic_with_symmetric_group is_isomorphic_to_symmetric_group
@deprecate has_is_isomorphic_with_symmetric_group has_is_isomorphic_to_symmetric_group
@deprecate set_is_isomorphic_with_symmetric_group set_is_isomorphic_to_symmetric_group

@deprecate is_isomorphic_with_alternating_group is_isomorphic_to_alternating_group
@deprecate has_is_isomorphic_with_alternating_group has_is_isomorphic_to_alternating_group
@deprecate set_is_isomorphic_with_alternating_group set_is_isomorphic_to_alternating_group

@deprecate proj_space projective_space

@deprecate are_algebraically_independent is_algebraically_independent_with_relations

@deprecate ambient_ring(U::AbsMultSet) ring(U)

# Deprecated for 1.0
@deprecate _compute_glueing_base_change _compute_gluing_base_change
@deprecate _compute_toric_glueing _compute_toric_gluing
@deprecate add_glueing! add_gluing!
@deprecate glueing_domains gluing_domains
@deprecate glueing_graph gluing_graph
@deprecate glueing_morphisms gluing_morphisms
@deprecate glueings gluings
@deprecate is_connected_glueing is_connected_gluing
@deprecate pruned_glueing_graph pruned_gluing_graph
@deprecate underlying_glueing underlying_gluing
@deprecate update_glueing_graph update_gluing_graph
Base.@deprecate_binding AbsGlueing AbsGluing
Base.@deprecate_binding BaseChangeGlueingData BaseChangeGluingData
Base.@deprecate_binding Glueing Gluing
Base.@deprecate_binding SimpleGlueing SimpleGluing
Base.@deprecate_binding ToricGlueingData ToricGluingData

@deprecate jacobi_matrix jacobian_matrix
@deprecate jacobi_ideal jacobian_ideal

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

Base.@deprecate_binding Spec AffineScheme
@deprecate proj(E::ToricLineBundle...) projectivization
@deprecate proj(E::ToricDivisor...) projectivization
Base.@deprecate_binding AbsSpec AbsAffineScheme
Base.@deprecate_binding SpecMor AffineSchemeMor
Base.@deprecate_binding AbsSpecMor AbsAffineSchemeMor
Base.@deprecate_binding SimplifiedSpec SimplifiedAffineScheme
Base.@deprecate_binding SpecOpen AffineSchemeOpenSubscheme
Base.@deprecate_binding SpecOpenMor AffineSchemeOpenSubschemeMor
Base.@deprecate_binding SpecOpenRing AffineSchemeOpenSubschemeRing
Base.@deprecate_binding SpecOpenRingElem AffineSchemeOpenSubschemeRingElem
Base.@deprecate_binding SpecSubset AffineSchemeSubset
Base.@deprecate_binding StdSpec StdAffineScheme

# deprecated for 1.1
@deprecate morphism_of_projective_schemes morphism

function Base.getindex(r::Hecke.SRow, R::AbstractAlgebra.Ring, u::AbstractUnitRange)
  Base.depwarn("`getindex(::SRow, ::Ring, ::AbstractUnitRange)` is deprecated, use `getindex(::SRow, ::AbstractUnitRange)` instead.", :getindex)
  @req base_ring(r) === R "Parent ring mismatch"
  return getindex(r, u)
end

@deprecate is_full_fp_group(G::FPGroup) GAPWrap.IsFpGroup(GapObj(G))

@deprecate minimal_generators minimal_generating_set
Base.@deprecate_binding MPolyRingElemLoc MPolyLocRingElem

# deprecated for 1.2
Base.@deprecate_binding QQAbElem QQAbFieldElem

Base.@deprecate_binding FreeAssAlgIdeal FreeAssociativeAlgebraIdeal

Base.@deprecate_binding in_linear_system is_in_linear_system
@deprecate scheme(W::AbsAlgebraicCycle) ambient_scheme(W)
@deprecate scheme(W::CartierDivisor) ambient_scheme(W)
@deprecate scheme(W::EffectiveCartierDivisor) ambient_scheme(W)

@deprecate mordell_weil_lattice(X::EllipticSurface) mordell_weil_sublattice(X) 
@deprecate minimal_generating_set(G::GAPGroup) minimal_size_generating_set(G)
@deprecate has_minimal_generating_set(G::GAPGroup) has_minimal_size_generating_set(G)
@deprecate set_minimal_generating_set(G::GAPGroup, v) set_minimal_size_generating_set(G, v)

# deprecated for 1.3
@deprecate acting_domain(C::GroupCoset) acting_group(C)
@deprecate acting_domain(Omega::GSet) acting_group(Omega)
@deprecate grid_morphism lattice_homomorphism

@deprecate gset(G::PermGroup) natural_gset(G)
@deprecate gset(G::MatrixGroup{T, MT}) where {MT, T <: FinFieldElem} natural_gset(G)
