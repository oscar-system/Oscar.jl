@alias has_isfinite has_is_finite
@alias SymmetricGroup symmetric_group

# make some Julia names compatible with our naming conventions
@alias is_subset issubset
@alias is_valid isvalid

# for backwards compatibility
@alias hall_subgroups_representatives hall_subgroup_reps
@alias hasrelshp has_relshp
@alias hastorusfactor has_torusfactor
@alias inner_automorphisms_group inner_automorphism_group
#@alias isabsolutely_irreducible is_absolutely_irreducible
@alias isaffine is_affine
@alias isalmostsimple is_almostsimple
@alias isalternating_form is_alternating_form
@alias isample is_ample
@alias isbicoset is_bicoset
@alias isbinomial is_binomial
@alias isbounded is_bounded
@alias iscartier is_cartier
@alias iscellular is_cellular
@alias iscomplete is_complete
@alias iscongruent is_congruent
@alias isconjugate_subgroup is_conjugate_subgroup
#@alias isdecomposable is_decomposable
@alias isdecorated is_decorated
@alias isdihedral_group is_dihedral_group
@alias isfano is_fano
@alias isfeasible is_feasible
@alias isfiltered is_filtered
@alias isfinite_order is_finiteorder
@alias isfinitelygenerated is_finitelygenerated
@alias isfull_direct_product is_full_direct_product
@alias isfull_semidirect_product is_full_semidirect_product
@alias isfull_wreath_product is_full_wreath_product
@alias isfulldimensional is_fulldimensional
@alias isgenerated_by_standard_unit_vectors is_generated_by_standard_unit_vectors
@alias isglobal is_global
@alias isgraded is_graded
@alias ishermitian_form is_hermitian_form
@alias ishermitian_matrix is_hermitian_matrix
@alias isinner_automorphism is_inner_automorphism
@alias isinvariant is_invariant
@alias isisomorphic_with_alternating_group is_isomorphic_with_alternating_group
@alias isisomorphic_with_symmetric_group is_isomorphic_with_symmetric_group
@alias isleft is_left
@alias islocal is_local
@alias ismixed is_mixed
@alias ismolien_series_implemented is_molien_series_implemented
@alias isnatural_alternating_group is_natural_alternating_group
@alias isnatural_symmetric_group is_natural_symmetric_group
@alias isnef is_nef
@alias isobviouslyabelian is_obviouslyabelian
@alias isorbifold is_orbifold
@alias isperfect is_perfect
@alias ispgroup is_pgroup
@alias ispointed is_pointed
@alias isprojective is_projective
@alias ispure is_pure
@alias isquadratic_form is_quadratic_form
@alias isquaternion_group is_quaternion_group
@alias isright is_right
@alias issemiregular is_semiregular
@alias issemisimple is_semisimple
@alias issimplicial is_simplicial
@alias issingular is_singular
@alias isskewsymmetric_matrix is_skewsymmetric_matrix
@alias issmooth_curve is_smooth_curve
@alias issolvable is_solvable
@alias issupersolvable is_supersolvable
@alias issymmetric_form is_symmetric_form
@alias istransitive is_transitive
@alias isunipotent is_unipotent
@alias isunital is_unital
@alias iswelldefined is_welldefined

# Allow backwards compatibility after removal of Oscar.Graphs module.
const Graphs = Oscar

# Compatibility with pre-0.12.x
@alias MPolyElem_dec MPolyDecRingElem
@alias MPolyRing_dec MPolyDecRing
@alias MPolyLocalizedRingElem MPolyLocRingElem
@alias MPolyLocalizedRing MPolyLocRing
@alias MPolyQuoElem MPolyQuoRingElem
@alias MPolyQuo MPolyQuoRing
@alias MPolyQuoLocalizedRingElem MPolyQuoLocRingElem
@alias MPolyQuoLocalizedRing MPolyQuoLocRing
@alias SubQuoElem SubquoModuleElem
@alias SubQuo SubquoModule
#@alias SubQuoElem_dec SubquoDecModuleElem
#@alias SubQuo_dec SubquoDecModule
@alias GradedPolynomialRing graded_polynomial_ring
