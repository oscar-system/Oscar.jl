# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
tries to adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.1](https://github.com/oscar-system/Oscar.jl/releases/tag/v1.3.1) - 2025-03-14

The following gives an overview of the changes compared to the previous release. This list is not
complete, many more internal or minor changes were made, but we tried to only list those changes
which we think might affect some users directly.

### Lie Theory

- [#4676](https://github.com/oscar-system/Oscar.jl/pull/4676) Fix `bracket` for `LieSubalgebra`s (the return value had the wrong type)

### Polyhedral Geometry

- [#4690](https://github.com/oscar-system/Oscar.jl/pull/4690) SubdivisionOfPoints: reject duplicate points

### Other fixed bugs

- [#4696](https://github.com/oscar-system/Oscar.jl/pull/4696) Restore Oscar<->Singular conversion for multivariate rational function fields

### Improvements or additions to documentation

- [#4701](https://github.com/oscar-system/Oscar.jl/pull/4701) Include linear solving doc page in the Oscar docs

### **TODO** release notes: to be added

If there are any PRs listed below, check their title and labels.
When done, change their label to "release notes: use title".

- [#4695](https://github.com/oscar-system/Oscar.jl/pull/4695) Fix Chevalley basis; add more functionality

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
