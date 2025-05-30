# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
tries to adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.0](https://github.com/oscar-system/Oscar.jl/releases/tag/v1.4.0) - 2025-05-30

The following gives an overview of the changes compared to the previous release. This list is not
complete, many more internal or minor changes were made, but we tried to only list those changes
which we think might affect some users directly.

### Algebraic Geometry

- [#4452](https://github.com/oscar-system/Oscar.jl/pull/4452) Add `graph_curve(::Graph)`

### Combinatorics

- [#4746](https://github.com/oscar-system/Oscar.jl/pull/4746) Import `Multipartition` functionality from JuLie
- [#4735](https://github.com/oscar-system/Oscar.jl/pull/4735) Add iterator for combinations
- [#4663](https://github.com/oscar-system/Oscar.jl/pull/4663) Add labelings for graphs
- [#3928](https://github.com/oscar-system/Oscar.jl/pull/3928) Add partially ordered sets
- [#3928](https://github.com/oscar-system/Oscar.jl/pull/3928) Add `indegree` and `outdegree` for graphs

### Commutative Algebra

- [#4890](https://github.com/oscar-system/Oscar.jl/pull/4890) Restrict `is_global`, `is_local` to total monomial orderings, add `is_(global|local)_block`
- [#4850](https://github.com/oscar-system/Oscar.jl/pull/4850) Add `quo` for `LaurentMPolyWrapRing`
- [#4810](https://github.com/oscar-system/Oscar.jl/pull/4810) Fix `tensor_product` for `SubquoModules`
- [#4765](https://github.com/oscar-system/Oscar.jl/pull/4765) New wrapper for Singular triangular decompositions
- [#4706](https://github.com/oscar-system/Oscar.jl/pull/4706) Improve comparison of morphisms of modules
- [#4683](https://github.com/oscar-system/Oscar.jl/pull/4683) Add hint for `eliminate` using a proper subring
- [#4626](https://github.com/oscar-system/Oscar.jl/pull/4626) Allow sorting w.r.t. monomial orderings / module orderings

### F-Theory Tools

- [#4876](https://github.com/oscar-system/Oscar.jl/pull/4876) Rename `global_gauge_quotients` to `global_gauge_group_quotient`
- [#4869](https://github.com/oscar-system/Oscar.jl/pull/4869) Update QSM Artifact
- [#4844](https://github.com/oscar-system/Oscar.jl/pull/4844) Improve printing of G4-flux and families thereof
- [#4844](https://github.com/oscar-system/Oscar.jl/pull/4844) Bug fix in computing smallest containing flux family for a given individual G4-flux
- [#4844](https://github.com/oscar-system/Oscar.jl/pull/4844) Update .mrdi-files (artifact) for model 1511.03209
- [#4779](https://github.com/oscar-system/Oscar.jl/pull/4779) Add exceptional classes and indices
- [#4728](https://github.com/oscar-system/Oscar.jl/pull/4728) Include 1/2 c2 in identification of flux families
- [#4711](https://github.com/oscar-system/Oscar.jl/pull/4711) Add convenience constructors for flux instances
- [#4694](https://github.com/oscar-system/Oscar.jl/pull/4694) Rename `is_vertical` to `passes_transversality_checks` and execute related checks

### Groups

- [#4921](https://github.com/oscar-system/Oscar.jl/pull/4921) Allow inducing G-sets along group homomorphisms
- [#4888](https://github.com/oscar-system/Oscar.jl/pull/4888) Document the relation between "abelian invariants" and "elementary divisors"
- [#4839](https://github.com/oscar-system/Oscar.jl/pull/4839) Add accessors `group` and `subgroup` for `SubgroupTransversal`
- [#4821](https://github.com/oscar-system/Oscar.jl/pull/4821) Export/document `character_table_complex_reflection_group`
- [#4771](https://github.com/oscar-system/Oscar.jl/pull/4771) Use `Set` for `GSet` block systems
- [#4698](https://github.com/oscar-system/Oscar.jl/pull/4698) Add `rank` and `torsion_free_rank` methods for permutation groups, pc groups, free groups
- [#4692](https://github.com/oscar-system/Oscar.jl/pull/4692) Add `GSet` methods for `blocks` and related functions
- [#4661](https://github.com/oscar-system/Oscar.jl/pull/4661) Add `dicyclic_group`, `is_dicyclic_group` and have `quaternion_group` and `is_quaternion_group` be aliases of those
- [#4659](https://github.com/oscar-system/Oscar.jl/pull/4659) Implement `transitivity` and `rank_action` for G-sets
- [#4628](https://github.com/oscar-system/Oscar.jl/pull/4628) Add local Schur indices for a character (Unger's algorithm)
- [#4609](https://github.com/oscar-system/Oscar.jl/pull/4609) Generalize G-sets to Weyl groups

### Lie Theory

- [#4878](https://github.com/oscar-system/Oscar.jl/pull/4878) Fix `symmetric_power` of a dim 0 module
- [#4807](https://github.com/oscar-system/Oscar.jl/pull/4807) Change default ordering in `universal_enveloping_algebra` to be admissible
- [#4789](https://github.com/oscar-system/Oscar.jl/pull/4789) Add `highest_root(R::RootSystem)` for convenience
- [#4729](https://github.com/oscar-system/Oscar.jl/pull/4729) Experimental: Add support for reducible types in `isomorphism(PermGroup, ::WeylGroup)`
- [#4687](https://github.com/oscar-system/Oscar.jl/pull/4687) Experimental: Add braid moves for words in Weyl groups
- [#4641](https://github.com/oscar-system/Oscar.jl/pull/4641) Experimental: Add `irreducible_factors` and `inner_direct_product` for Weyl groups

### Number Theory

- [#4905](https://github.com/oscar-system/Oscar.jl/pull/4905) Fix `pc_group_with_isomorphism(::FinGenAbGroup)`
- [#4903](https://github.com/oscar-system/Oscar.jl/pull/4903) Fix irreducibility test for `AbsSimpleNumFieldOrderElem`
- [#4879](https://github.com/oscar-system/Oscar.jl/pull/4879) Add `is_perfect` for algebraic closures of finite fields
- [#4837](https://github.com/oscar-system/Oscar.jl/pull/4837) Fix `maximal_order` for `NfNSGen`
- [#4740](https://github.com/oscar-system/Oscar.jl/pull/4740) Add `degree_of_character_field`

### Tropical Geometry

- [#4838](https://github.com/oscar-system/Oscar.jl/pull/4838) Fix `tropical_variety_zerodimensional`
- [#4781](https://github.com/oscar-system/Oscar.jl/pull/4781) Add `roots` for tropical polynomials
- [#4703](https://github.com/oscar-system/Oscar.jl/pull/4703) Add tropical prevarieties generated by intersecting tropical hypersurfaces
- [#4697](https://github.com/oscar-system/Oscar.jl/pull/4697) Remove broken tropical Groebner basis shortcut for binomial ideals
- [#4447](https://github.com/oscar-system/Oscar.jl/pull/4447) Add `positive_tropical_variety`
- [#4445](https://github.com/oscar-system/Oscar.jl/pull/4445) Add tropical linear spaces from graphs
- [#4061](https://github.com/oscar-system/Oscar.jl/pull/4061) Overhaul tropical varieties, add various new options

### Changes related to serializing data in the MRDI file format

- [#4331](https://github.com/oscar-system/Oscar.jl/pull/4331) Unify type encoding for similar types
- [#4162](https://github.com/oscar-system/Oscar.jl/pull/4162) Cleaner handling of type parameter serialization. This update forces entries of container types to share the same output of `Oscar.type_params` when serializing. Deserialization speed improvements.

### New features or extended functionality

- [#4797](https://github.com/oscar-system/Oscar.jl/pull/4797) Add experimental support for wreath Macdonald polynomials

### Only changes experimental parts of OSCAR

- [#4899](https://github.com/oscar-system/Oscar.jl/pull/4899) IntersectionTheory: Rename `abstract_projective_bundle` -> `projective_bundle` and `abstract_flag_bundle` -> `flag_bundle`
- [#4845](https://github.com/oscar-system/Oscar.jl/pull/4845) GroebnerWalk: Remove perturbed walk
- [#4780](https://github.com/oscar-system/Oscar.jl/pull/4780) Oscar Worker Pool and parallel functions functionality such as pmap.
- [#4772](https://github.com/oscar-system/Oscar.jl/pull/4772) Intersection theory: Introduce Gromov-Witten invariants
- [#4769](https://github.com/oscar-system/Oscar.jl/pull/4769) Intersection theory: extend documentation on Bott formula
- [#4764](https://github.com/oscar-system/Oscar.jl/pull/4764) Intersection Theory: Kontsevich spaces
- [#4100](https://github.com/oscar-system/Oscar.jl/pull/4100) Injective and irreducible resolutions of Q-graded modules
- [#2183](https://github.com/oscar-system/Oscar.jl/pull/2183) Add basics for quantum analogs

### Improvements or additions to documentation

- [#4758](https://github.com/oscar-system/Oscar.jl/pull/4758) Collapse docstrings in documentation to allow for easier navigation

### Changes related to the package AbstractAlgebra

- [#4894](https://github.com/oscar-system/Oscar.jl/pull/4895) Bump AbstractAlgebra to v0.45

### Changes related to the package AlgebraicSolving

- [#4894](https://github.com/oscar-system/Oscar.jl/pull/4894) Bump AlgebraicSolving to v0.9.0

### Changes related to the package Hecke

- [#4894](https://github.com/oscar-system/Oscar.jl/pull/4895) Bump Hecke to v0.36

### Changes related to the package Nemo

- [#4894](https://github.com/oscar-system/Oscar.jl/pull/4895) Bump Nemo to v0.50

### Changes related to the package Polymake

- [#4894](https://github.com/oscar-system/Oscar.jl/pull/4895) Bump Polymake.jl to v0.12

### Changes related to the package Singular

- [#4613](https://github.com/oscar-system/Oscar.jl/pull/4613) Optimize conversions from/to Singular

## [1.3.1](https://github.com/oscar-system/Oscar.jl/releases/tag/v1.3.1) - 2025-03-14

The following gives an overview of the changes compared to the previous release. This list is not
complete, many more internal or minor changes were made, but we tried to only list those changes
which we think might affect some users directly.

### Lie Theory

- [#4676](https://github.com/oscar-system/Oscar.jl/pull/4676) Fix `bracket` for `LieSubalgebra`s (the return value had the wrong type)
- [#4695](https://github.com/oscar-system/Oscar.jl/pull/4695) Fix `chevalley_basis` (sometimes returning wrong results)
- [#4695](https://github.com/oscar-system/Oscar.jl/pull/4695) Add `change_base_ring(::Field, ::LieAlgebra)`
- [#4695](https://github.com/oscar-system/Oscar.jl/pull/4695) Add `structure_constant_table(::LieAlgebra, basis)`

### Polyhedral Geometry

- [#4690](https://github.com/oscar-system/Oscar.jl/pull/4690) SubdivisionOfPoints: reject duplicate points

### Other fixed bugs

- [#4696](https://github.com/oscar-system/Oscar.jl/pull/4696) Restore Oscar<->Singular conversion for multivariate rational function fields

### Improvements or additions to documentation

- [#4701](https://github.com/oscar-system/Oscar.jl/pull/4701) Include linear solving doc page in the Oscar docs

## [1.3.0](https://github.com/oscar-system/Oscar.jl/releases/tag/v1.3.0) - 2025-02-28

The following gives an overview of the changes compared to the previous release. This list is not
complete, many more internal or minor changes were made, but we tried to only list those changes
which we think might affect some users directly.

### Highlights

- [#4399](https://github.com/oscar-system/Oscar.jl/pull/4399) Graduate Weyl groups and root systems from experimental to officially supported
- [#4294](https://github.com/oscar-system/Oscar.jl/pull/4294) Graduate elliptic surfaces from experimental

### Renamings

- [#4526](https://github.com/oscar-system/Oscar.jl/pull/4526) Rename `blowup` and `blow_up_points`  to ` blow_up` and `blow_up_points`, and other intersection theory changes
- [#4433](https://github.com/oscar-system/Oscar.jl/pull/4433) Rename `new_ray` to `exceptional_ray`
- [#4304](https://github.com/oscar-system/Oscar.jl/pull/4304) Change `hom` to `map` in in intersection theory
- [#4302](https://github.com/oscar-system/Oscar.jl/pull/4302) Rename `acting_domain` to `acting_group` and some improvements for group cosets
- [#4278](https://github.com/oscar-system/Oscar.jl/pull/4278) Rename `minimal_generating_set` for groups to `minimal_size_generating_set`

### Algebraic Geometry

- [#4615](https://github.com/oscar-system/Oscar.jl/pull/4615) Speed up computations for EllipticSurface by reduction to positive characteristic.
- [#4560](https://github.com/oscar-system/Oscar.jl/pull/4560) Add a database of Enriques surfaces
- [#4540](https://github.com/oscar-system/Oscar.jl/pull/4540) Support for computing automorphism groups of Enriques surfaces
- [#4534](https://github.com/oscar-system/Oscar.jl/pull/4534) Add `stabilizer_in_orthogonal_group` to compute stabilisers of timelike vectors
- [#4485](https://github.com/oscar-system/Oscar.jl/pull/4485) Make `saturation` for principal ideals faster by delegating to `remove`
- [#4393](https://github.com/oscar-system/Oscar.jl/pull/4393) Introduce `-inf` as potential output for the dimension of an ideal
- [#4352](https://github.com/oscar-system/Oscar.jl/pull/4352) Add documentation for coherent sheaves
- [#4345](https://github.com/oscar-system/Oscar.jl/pull/4345) Add a prototype for a moving lemma in concrete intersection theory
- [#4327](https://github.com/oscar-system/Oscar.jl/pull/4327) Add Eagon-Northcott complexes
- [#4272](https://github.com/oscar-system/Oscar.jl/pull/4272) Add Hasse-Schmidt derivatives
- [#4256](https://github.com/oscar-system/Oscar.jl/pull/4256) Print Betti tables more nicely
- [#4183](https://github.com/oscar-system/Oscar.jl/pull/4183) Add tweaks for enabling a 2-neighbor-step in characteristic 0

### Combinatorics

- [#4633](https://github.com/oscar-system/Oscar.jl/pull/4633) Fix `vertices(G)` to return all vertices, including isolated ones
- [#4499](https://github.com/oscar-system/Oscar.jl/pull/4499) Add `is_bipartite(::Graph)`
- [#4459](https://github.com/oscar-system/Oscar.jl/pull/4459) Make vertex labels in graph visualization start at 1
- [#4450](https://github.com/oscar-system/Oscar.jl/pull/4450) Add `in(::Int, ::Edge)`, `signed_incidence_matrix(::Graph{Undirected})`, `connectivity(::Graph{Undirected})`
- [#4441](https://github.com/oscar-system/Oscar.jl/pull/4441) Fix `degree(::Graph)` documentation
- [#4299](https://github.com/oscar-system/Oscar.jl/pull/4299) Fix leading zero bug in `matroid_hex`
- [#4270](https://github.com/oscar-system/Oscar.jl/pull/4270) Fix `BoundError` in `weights(hook_lengths(...))`

### Commutative Algebra

- [#4601](https://github.com/oscar-system/Oscar.jl/pull/4601) Fix Krull dimension of polynomial rings over number fields
- [#4596](https://github.com/oscar-system/Oscar.jl/pull/4596) Support computing Gröbner basis for lattice ideals using 4ti2, and improve `eliminate` to avoid recomputing a Gröbner basis
- [#4497](https://github.com/oscar-system/Oscar.jl/pull/4497) Allow variables of same positive degree in `monomials_of_degree`
- [#4492](https://github.com/oscar-system/Oscar.jl/pull/4492) Add links to commutative algebra tutorials
- [#4488](https://github.com/oscar-system/Oscar.jl/pull/4488) Overhaul preprocessing for `radical`, `primary_decomposition`, and friends over number fields
- [#4468](https://github.com/oscar-system/Oscar.jl/pull/4468) Graduate `present_finite_extension_ring` from experimental to supported (and fix a bug in it)
- [#4456](https://github.com/oscar-system/Oscar.jl/pull/4456) Generalize `monomials_of_degree` to allow graded rings with grading group Z
- [#4400](https://github.com/oscar-system/Oscar.jl/pull/4400) Speed up evaluation of maps of `MPolyRing`s which take variables to variables 
- [#4379](https://github.com/oscar-system/Oscar.jl/pull/4379) Fix bug in comparison function in module orderings
- [#4346](https://github.com/oscar-system/Oscar.jl/pull/4346) Fix `characteristic` method for localized rings
- [#4248](https://github.com/oscar-system/Oscar.jl/pull/4248) Add Cartan eilenberg resolutions of complexes
- [#4235](https://github.com/oscar-system/Oscar.jl/pull/4235) Add `monomial_basis` for `MPolyQuoLocRing`

### F-Theory Tools

- [#4636](https://github.com/oscar-system/Oscar.jl/pull/4636) Add attribute for more detailed info on tunable sections
- [#4597](https://github.com/oscar-system/Oscar.jl/pull/4597) Update `tune` function and related definitions/functions
- [#4562](https://github.com/oscar-system/Oscar.jl/pull/4562) Align names of properties and attributes among families of and individual G4-fluxes
- [#4503](https://github.com/oscar-system/Oscar.jl/pull/4503) Fix and extend strict_transform
- [#4500](https://github.com/oscar-system/Oscar.jl/pull/4500) Extend support for G4-flux families and individual G4-fluxes
- [#4491](https://github.com/oscar-system/Oscar.jl/pull/4491) Save updated intersection numbers after costly computation
- [#4487](https://github.com/oscar-system/Oscar.jl/pull/4487) Compute D3-tadpole for family of G4-fluxes
- [#4466](https://github.com/oscar-system/Oscar.jl/pull/4466) Introduce families of G4-fluxes
- [#4446](https://github.com/oscar-system/Oscar.jl/pull/4446) Add support for all vertical, well-quantized G4s that do not break the non-abelian gauge group
- [#4423](https://github.com/oscar-system/Oscar.jl/pull/4423) Use Zenodo data as artifact
- [#4422](https://github.com/oscar-system/Oscar.jl/pull/4422) Print known properties of G4-fluxes
- [#4343](https://github.com/oscar-system/Oscar.jl/pull/4343) Extend support for zero section and zero section class
- [#4286](https://github.com/oscar-system/Oscar.jl/pull/4286) Support computation of all well-quantized and vertical G4-fluxes
- [#4268](https://github.com/oscar-system/Oscar.jl/pull/4268) Implement ambient space models for G4-flux candidates
- [#4243](https://github.com/oscar-system/Oscar.jl/pull/4243) Support computation of basis of H^(2,2) and H^4

### Groups

- [#4608](https://github.com/oscar-system/Oscar.jl/pull/4608) Add `is_primitive` for $G$-sets
- [#4594](https://github.com/oscar-system/Oscar.jl/pull/4594) Add $G$-set docstrings to main documentation
- [#4582](https://github.com/oscar-system/Oscar.jl/pull/4582) Cache parents for permutations in `perm`, `cperm` and `@perm`
- [#4490](https://github.com/oscar-system/Oscar.jl/pull/4490) Make the `@perm` macro more powerful by supporting more argument variants and increased consistency
- [#4469](https://github.com/oscar-system/Oscar.jl/pull/4469) Add dedicated `stabilizer` methods for matrix groups for improved performance
- [#4465](https://github.com/oscar-system/Oscar.jl/pull/4465) Add `orbit_representatives_and_stabilizers` for 0-dim. subspaces
- [#4416](https://github.com/oscar-system/Oscar.jl/pull/4416) Add `is_atlas_character_table`
- [#4411](https://github.com/oscar-system/Oscar.jl/pull/4411) Support `character_field` for a vector of characters
- [#4409](https://github.com/oscar-system/Oscar.jl/pull/4409) Add the library of groups with at most 14 conjugacy classes
- [#4401](https://github.com/oscar-system/Oscar.jl/pull/4401) Add `isomorphic_subgroups`
- [#4378](https://github.com/oscar-system/Oscar.jl/pull/4378) Fix wrong result of `isomorphism(FPGroup, G, on_gens = true)` for trivial `G` with more than 0 generators
- [#4361](https://github.com/oscar-system/Oscar.jl/pull/4361) Speed up `cperm`
- [#4359](https://github.com/oscar-system/Oscar.jl/pull/4359) Add `cycle_length` for `PermGroupElem`
- [#4357](https://github.com/oscar-system/Oscar.jl/pull/4357) Improve `isomorphism` from pc-groups and fp-groups to `FinGenAbGroup`
- [#4337](https://github.com/oscar-system/Oscar.jl/pull/4337) Speed up `orbit_representatives_and_stabilizers` for not too large examples
- [#4319](https://github.com/oscar-system/Oscar.jl/pull/4319) Make `isomorphism(PcGroup, A)` for infinite abelian `A` work
- [#4311](https://github.com/oscar-system/Oscar.jl/pull/4311) Add action on matrices in row reduced echelon form
- [#4307](https://github.com/oscar-system/Oscar.jl/pull/4307) Speed up orbits of permutation groups on integers
- [#4298](https://github.com/oscar-system/Oscar.jl/pull/4298) Add `is_left` and `is_right` for `SubgroupTransversal`
- [#4281](https://github.com/oscar-system/Oscar.jl/pull/4281) Speed up `stabilizer` methods for the natural action of permutation groups, and the actions on sets, tuples, and vectors of pos. integers
- [#4265](https://github.com/oscar-system/Oscar.jl/pull/4265) Add `stabilizer` for $G$-sets with action on `Set`, `Vector` or `Tuple` objects
- [#4259](https://github.com/oscar-system/Oscar.jl/pull/4259) Add argument check to $G$-set constructor
- [#4202](https://github.com/oscar-system/Oscar.jl/pull/4202) Add `letters` function for `PcGroupElem`

### Lie Theory

- [#4676](https://github.com/oscar-system/Oscar.jl/pull/4676) Fix `bracket` for `LieSubalgebra`s (it sometimes returned the wrong result)
- [#4621](https://github.com/oscar-system/Oscar.jl/pull/4621) Add `(dual_)geometric_representation(::WeylGroup)`, `is_finite_order(::WeylGroupElem)` and `order(::WeylGroupElem)`
- [#4589](https://github.com/oscar-system/Oscar.jl/pull/4589) Implement `syllables`, `letters` and their inverse for `WeylGroup` and `WeylGroupElem`
- [#4536](https://github.com/oscar-system/Oscar.jl/pull/4536) Add `weyl_group(::Matrix{<:IntegerUnion})` convenience constructor
- [#4519](https://github.com/oscar-system/Oscar.jl/pull/4519) Add `map_word` and `parabolic_subgroup` for Weyl groups
- [#4478](https://github.com/oscar-system/Oscar.jl/pull/4478) Add `isomorphism(PermGroup, ::WeylGroup)` for missing irreducible types
- [#4414](https://github.com/oscar-system/Oscar.jl/pull/4414) Add `WeightLattice` serialization
- [#4374](https://github.com/oscar-system/Oscar.jl/pull/4374) Change Weyl group action to a right action
- [#4339](https://github.com/oscar-system/Oscar.jl/pull/4339) Add `demazure_character`
- [#4264](https://github.com/oscar-system/Oscar.jl/pull/4264) Add `isomorphism(::Type{PermGroup}, ::WeylGroup)` for some cases

### Number Theory

- [#4585](https://github.com/oscar-system/Oscar.jl/pull/4585) Improve usability of `abelian_closure(QQ)`, e.g. for elements representing real values allow comparisons and conversion to `Float64`
- [#4429](https://github.com/oscar-system/Oscar.jl/pull/4429) Fix `R()` when `R` is a `BoundedRing`
- [#4424](https://github.com/oscar-system/Oscar.jl/pull/4424) Fix `is_irreducible` for number field order elements
- [#4396](https://github.com/oscar-system/Oscar.jl/pull/4396) Fix regression in `galois_group`
- [#4370](https://github.com/oscar-system/Oscar.jl/pull/4370) Experimental: Add Clifford algebras and Clifford orders

### Polyhedral Geometry

- [#4482](https://github.com/oscar-system/Oscar.jl/pull/4482) Add `integer_hull` and `gomory_chvatal_closure`
- [#4454](https://github.com/oscar-system/Oscar.jl/pull/4454) Fix blowups along rays and singular cones
- [#4451](https://github.com/oscar-system/Oscar.jl/pull/4451) Add `tutte_lifting(::Graph)`
- [#4392](https://github.com/oscar-system/Oscar.jl/pull/4392) Allow checking containment of points in hyperplanes and halfspaces via `in`
- [#4260](https://github.com/oscar-system/Oscar.jl/pull/4260) Experimental: Add algebraic shifting

### Toric Geometry

- [#4529](https://github.com/oscar-system/Oscar.jl/pull/4529) Add lattice of one-parameter subgroups
- [#4509](https://github.com/oscar-system/Oscar.jl/pull/4509) Implement equality for toric varieties
- [#4508](https://github.com/oscar-system/Oscar.jl/pull/4508) Change diagonal to barycenter in toric blowup docs
- [#4505](https://github.com/oscar-system/Oscar.jl/pull/4505) Simplify code for toric blowup along a cone
- [#4484](https://github.com/oscar-system/Oscar.jl/pull/4484) Improve performance of strict transform under toric blowups
- [#4437](https://github.com/oscar-system/Oscar.jl/pull/4437) Support blow up along minimal supercone coordinates
- [#4266](https://github.com/oscar-system/Oscar.jl/pull/4266) Serialize cohomology classes
- [#4154](https://github.com/oscar-system/Oscar.jl/pull/4154) Support computing strict transform of an ideal in Cox ring

### Improvements or additions to documentation

- [#4549](https://github.com/oscar-system/Oscar.jl/pull/4549) Fix experimental not showing up in search
- [#4480](https://github.com/oscar-system/Oscar.jl/pull/4480) Add more links to tutorials
- [#4458](https://github.com/oscar-system/Oscar.jl/pull/4458) Add `weights` method for graded rings and improve some documentation

### Changes related to the package AbstractAlgebra

- [#4405](https://github.com/oscar-system/Oscar.jl/pull/4405) Update AbstractAlgebra to 0.44.0

### Changes related to the package GAP

- [#4421](https://github.com/oscar-system/Oscar.jl/pull/4421) Update GAP.jl to 0.13

### Changes related to the package Hecke

- [#4405](https://github.com/oscar-system/Oscar.jl/pull/4405) Update Hecke to 0.35.0

### Changes related to the package Nemo

- [#4645](https://github.com/oscar-system/Oscar.jl/pull/4645) Update Nemo to 0.49

### Changes related to the package Singular

- [#4548](https://github.com/oscar-system/Oscar.jl/pull/4548) Update Singular.jl to 0.25.0


## [1.2.2] - 2024-12-13

### Fixed

- Fixed a problem with `galois_group` [PR #4396](https://github.com/oscar-system/Oscar.jl/pull/4396)

- Fixed zero-dimensional cone in `cones` in PolyhedralGeometry. [PR
  #4336](https://github.com/oscar-system/Oscar.jl/pull/4336)

- Fixed up generic characteristic method for localized rings. [issue
  #4324](https://github.com/oscar-system/Oscar.jl/issues/4324) [PR
  #4346](https://github.com/oscar-system/Oscar.jl/pull/4346)

- Added hash method for `RayVector`. [partial fix for issue
  #2222](https://github.com/oscar-system/Oscar.jl/issues/2222) [PR
  #4354](https://github.com/oscar-system/Oscar.jl/pull/4354)
