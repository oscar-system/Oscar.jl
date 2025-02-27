# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
tries to adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.3](https://github.com/oscar-system/Oscar.jl/releases/tag/v1.2.3) - 2025-02-27

The following gives an overview of the changes compared to the previous release. This list is not
complete, many more internal or minor changes were made, but we tried to only list those changes
which we think might affect some users directly.

### Items being renamed

- [#4433](https://github.com/oscar-system/Oscar.jl/pull/4433) Rename `new_ray` to `exceptional_ray`

### Changes related to Algebraic Geometry

- [#4615](https://github.com/oscar-system/Oscar.jl/pull/4615) Speed up computations for EllipticSurface by reduction to positive characteristic.
- [#4560](https://github.com/oscar-system/Oscar.jl/pull/4560) Add a database of Enriques surfaces
- [#4485](https://github.com/oscar-system/Oscar.jl/pull/4485) Make `saturation` for principal ideals faster by delegating to `remove`

### Changes related to Combinatorics

- [#4633](https://github.com/oscar-system/Oscar.jl/pull/4633) Fix `vertices(G)` to return all vertices, including isolated ones
- [#4499](https://github.com/oscar-system/Oscar.jl/pull/4499) Add `is_bipartite(::Graph)`
- [#4459](https://github.com/oscar-system/Oscar.jl/pull/4459) Make vertex labels in graph visualization start at 1
- [#4450](https://github.com/oscar-system/Oscar.jl/pull/4450) Add `in(::Int, ::Edge)`, `signed_incidence_matrix(::Graph{Undirected})`, `connectivity(::Graph{Undirected})`
- [#4441](https://github.com/oscar-system/Oscar.jl/pull/4441) Fix `degree(::Graph)` documentation

### Changes related to Commutative Algebra

- [#4601](https://github.com/oscar-system/Oscar.jl/pull/4601) Fix Krull dimension of polynomial rings over number fields
- [#4596](https://github.com/oscar-system/Oscar.jl/pull/4596) Support computing Gröbner basis for lattice ideals using 4ti2, and improve `eliminate` to avoid recomputing a Gröbner basis

### Changes related to F-Theory Tools

- [#4647](https://github.com/oscar-system/Oscar.jl/pull/4647) Fix bug in G4-flux computation of big model
- [#4636](https://github.com/oscar-system/Oscar.jl/pull/4636) Added attribute for more detailed info on tunable sections
- [#4597](https://github.com/oscar-system/Oscar.jl/pull/4597) Update `tune` function and related definitions/functions
- [#4491](https://github.com/oscar-system/Oscar.jl/pull/4491) Save updated intersection numbers after costly computation
- [#4487](https://github.com/oscar-system/Oscar.jl/pull/4487) [FTheoryTools] Compute D3-tadpole for family of G4-fluxes
- [#4446](https://github.com/oscar-system/Oscar.jl/pull/4446) Add support for all vertical, well-quantized G4s that do not break the non-Abelian gauge group
- [#4422](https://github.com/oscar-system/Oscar.jl/pull/4422) Print known properties of G4-fluxes

### Changes related to Groups

- [#4608](https://github.com/oscar-system/Oscar.jl/pull/4608) Add `is_primitive` for $G$-sets
- [#4594](https://github.com/oscar-system/Oscar.jl/pull/4594) Add $G$-set docstrings to main documentation
- [#4582](https://github.com/oscar-system/Oscar.jl/pull/4582) Cache parents for permutations in `perm`, `cperm` and `@perm`
- [#4490](https://github.com/oscar-system/Oscar.jl/pull/4490) Make the `@perm` macro more powerful by supporting more argument variants and increased consistency
- [#4469](https://github.com/oscar-system/Oscar.jl/pull/4469) Add dedicated `stabilizer` methods for matrix groups for improved performance
- [#4465](https://github.com/oscar-system/Oscar.jl/pull/4465) Add `orbit_representatives_and_stabilizers` for 0-dim. subspaces
- [#4409](https://github.com/oscar-system/Oscar.jl/pull/4409) Add the library of groups with at most 14 conjugacy classes

### Changes related to Lie Theory

- [#4621](https://github.com/oscar-system/Oscar.jl/pull/4621) Add `(dual_)geometric_representation(::WeylGroup)`, `is_finite_order(::WeylGroupElem)` and `order(::WeylGroupElem)`
- [#4589](https://github.com/oscar-system/Oscar.jl/pull/4589) Implement `syllables`, `letters` and their inverse for `WeylGroup` and `WeylGroupElem`
- [#4536](https://github.com/oscar-system/Oscar.jl/pull/4536) Add `weyl_group(::Matrix{<:IntegerUnion})` convenience constructor
- [#4519](https://github.com/oscar-system/Oscar.jl/pull/4519) Add `map_word` and `parabolic_subgroup` for Weyl groups
- [#4478](https://github.com/oscar-system/Oscar.jl/pull/4478) Add `isomorphism(PermGroup, ::WeylGroup)` for missing irreducible types

### Changes related to Number Theory

- [#4429](https://github.com/oscar-system/Oscar.jl/pull/4429) Fix `R()` when `R` is a `BoundedRing`
- [#4424](https://github.com/oscar-system/Oscar.jl/pull/4424) Fix `is_irreducible` for number field order elements
- [#4370](https://github.com/oscar-system/Oscar.jl/pull/4370) Experimental: Add Clifford algebras and Clifford orders

### Changes related to Polyhedral Geometry

- [#4454](https://github.com/oscar-system/Oscar.jl/pull/4454) Fix blowups along rays and singular cones
- [#4451](https://github.com/oscar-system/Oscar.jl/pull/4451) Add `tutte_lifting(::Graph)`

### Changes related to Toric Schemes

- [#4505](https://github.com/oscar-system/Oscar.jl/pull/4505) Simplify code for toric blowup along a cone

### Changes related to Toric Varieties

- [#4508](https://github.com/oscar-system/Oscar.jl/pull/4508) Change diagonal to barycenter in toric blowup docs

### New features or extended functionality

- [#4497](https://github.com/oscar-system/Oscar.jl/pull/4497) Allow variables of same positive degree in `monomials_of_degree`
- [#4492](https://github.com/oscar-system/Oscar.jl/pull/4492) Add links to commutative algebra tutorials

### Changes related to the package AbstractAlgebra

- [#4405](https://github.com/oscar-system/Oscar.jl/pull/4405) Update AbstractAlgebra to 0.44.0, Nemo to 0.48.0, Hecke to 0.35.0

### Changes related to the package GAP

- [#4421](https://github.com/oscar-system/Oscar.jl/pull/4421) Update GAP.jl to 0.13

### Changes related to the package Nemo

- [#4645](https://github.com/oscar-system/Oscar.jl/pull/4645) Update Nemo to 0.49

### Changes related to the package Singular

- [#4548](https://github.com/oscar-system/Oscar.jl/pull/4548) Update Singular.jl to 0.25.0

### **TODO** release notes: to be added

If there are any PRs listed below, check their title and labels.
When done, change their label to "release notes: use title".

- [#4529](https://github.com/oscar-system/Oscar.jl/pull/4529) [ToricVarieties] Add lattice of one-parameter subgroups
- [#4522](https://github.com/oscar-system/Oscar.jl/pull/4522) Update docs subsection title
- [#4503](https://github.com/oscar-system/Oscar.jl/pull/4503) Fix FTheoryTools bug, extend strict_transform
- [#4500](https://github.com/oscar-system/Oscar.jl/pull/4500) [FTheoryTools] More improvements
- [#4486](https://github.com/oscar-system/Oscar.jl/pull/4486) [FTheoryTools] More improvements
- [#4484](https://github.com/oscar-system/Oscar.jl/pull/4484) Fast strict transform under toric blowups
- [#4482](https://github.com/oscar-system/Oscar.jl/pull/4482) new functions for integer programming
- [#4481](https://github.com/oscar-system/Oscar.jl/pull/4481) `QuadFormAndIsom`: fix Hermitian Miranda--Morrison (again...) + some cleanups
- [#4480](https://github.com/oscar-system/Oscar.jl/pull/4480) Add more links to tutorials
- [#4468](https://github.com/oscar-system/Oscar.jl/pull/4468) Fix bug, move present_finite_extension_ring out of experimental
- [#4466](https://github.com/oscar-system/Oscar.jl/pull/4466) [FTheoryTools] More improvements
- [#4458](https://github.com/oscar-system/Oscar.jl/pull/4458) some tweaks in docu, add method weights for graded rings
- [#4456](https://github.com/oscar-system/Oscar.jl/pull/4456) Generalize monomials_of_degree
- [#4437](https://github.com/oscar-system/Oscar.jl/pull/4437) Blow up along minimal supercone coordinates
- [#4430](https://github.com/oscar-system/Oscar.jl/pull/4430) Add `hash` methods for SLPolys
- [#4423](https://github.com/oscar-system/Oscar.jl/pull/4423) FTheoryTools: Use Zenodo data as artifact

### **TODO** Uncategorized PR

If there are any PRs listed below, either apply the same steps
as above, or change their label to "release notes: not needed".

- [#4587](https://github.com/oscar-system/Oscar.jl/pull/4587) Add missing division ops for QQab
- [#4577](https://github.com/oscar-system/Oscar.jl/pull/4577) Introduce global gap variable `Oscar_jl`
- [#4566](https://github.com/oscar-system/Oscar.jl/pull/4566) Internal variables names in toric blowups
- [#4565](https://github.com/oscar-system/Oscar.jl/pull/4565) Remove toric blowup arithmetic documentation
- [#4562](https://github.com/oscar-system/Oscar.jl/pull/4562) [FTheoryTools] Rename some functions
- [#4555](https://github.com/oscar-system/Oscar.jl/pull/4555) [FTheoryTools] More improvements
- [#4549](https://github.com/oscar-system/Oscar.jl/pull/4549) Fix experimental not showing up in search
- [#4540](https://github.com/oscar-system/Oscar.jl/pull/4540) Automorphism Groups of Enriques surfaces
- [#4534](https://github.com/oscar-system/Oscar.jl/pull/4534) Stabilisers of timelike vectors, some refactoring for K3Auto.jl
- [#4530](https://github.com/oscar-system/Oscar.jl/pull/4530) Simplify ToricBlowupMorphism
- [#4526](https://github.com/oscar-system/Oscar.jl/pull/4526) Intersection theory: fix issue #4512 and more
- [#4523](https://github.com/oscar-system/Oscar.jl/pull/4523) Remove toric blowups along ideals
- [#4509](https://github.com/oscar-system/Oscar.jl/pull/4509) Implement equality for toric varieties
- [#4507](https://github.com/oscar-system/Oscar.jl/pull/4507) Introduce `is_known` according to #4394.
- [#4488](https://github.com/oscar-system/Oscar.jl/pull/4488) Overhaul preprocessing for radical computations.
- [#4452](https://github.com/oscar-system/Oscar.jl/pull/4452) AlgebraicGeometry: added graph_curve(::Graph)
- [#4418](https://github.com/oscar-system/Oscar.jl/pull/4418) Update copyright and book publication year to 2025
- [#4417](https://github.com/oscar-system/Oscar.jl/pull/4417) Add conformance tests for `VarietyFunctionField`, add `base_ring_type` method
- [#4416](https://github.com/oscar-system/Oscar.jl/pull/4416) add `is_atlas_character_table`
- [#4415](https://github.com/oscar-system/Oscar.jl/pull/4415) fix two obvious duplicates in the list of groups with 14 classes
- [#4414](https://github.com/oscar-system/Oscar.jl/pull/4414) Move LieTheory serialization to src; add `WeightLattice` serialization
- [#4413](https://github.com/oscar-system/Oscar.jl/pull/4413) Caching of radicals
- [#4411](https://github.com/oscar-system/Oscar.jl/pull/4411) support `character_field` for a vector of characters
- [#4403](https://github.com/oscar-system/Oscar.jl/pull/4403) Blowup along a ray does not change orbifoldness
- [#4401](https://github.com/oscar-system/Oscar.jl/pull/4401) add `isomorphic_subgroups`
- [#4400](https://github.com/oscar-system/Oscar.jl/pull/4400) Further tuning on flattenings of rings
- [#4399](https://github.com/oscar-system/Oscar.jl/pull/4399) Move Weyl groups and root systems to `src/`
- [#4392](https://github.com/oscar-system/Oscar.jl/pull/4392) [PolyhedralGeometry] Check containment of points in hyperplanes and halfspaces
- [#4384](https://github.com/oscar-system/Oscar.jl/pull/4384) LieAlgebras: Adapt Demazure operator
- [#4327](https://github.com/oscar-system/Oscar.jl/pull/4327) Eagon northcott complex
- [#4154](https://github.com/oscar-system/Oscar.jl/pull/4154) Strict transform of an ideal in Cox ring


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
