# Changes — experimental/InjectiveResolutions

Summary of the changes to the `InjectiveResolutions` experimental package
since the version in OSCAR master (last upstream change: PR #5497), to
accompany the pull request. The package computes injective resolutions and
local cohomology of finitely generated modules over monoid algebras following
[HM05].

## New functionality

- **Affine semigroups.** New type `AffineSemigroup` with constructors
  `affine_semigroup` (from a matrix, a list of generators, or a monoid
  algebra), `semigroup_generators`, `ambient_dimension`, `is_pointed`, and
  cached polyhedral data (cone, bounding hyperplanes, zonotope). Monoid
  algebras store their semigroup; `monoid_algebra(Q::AffineSemigroup, k)`.
- **Non-normal monoid algebras.** `irreducible_resolution` and
  `injective_resolution` work for non-normal (unsaturated) `k[Q]`, using
  Algorithm 3.15 of [HM05] to compute the irreducible ideals in `k[Q]` from
  those in the saturation. Supporting functions: `saturation(kQ)`,
  `saturation_map`, `saturation_ideal`, `holes_module` (the module
  `k[Q_sat]/k[Q]`), `is_Q_graded`, semigroup membership tests.
  Local cohomology still requires a normal monoid algebra.
- **Corrected coefficient computation.** Algorithm 3.6 of [HM05] can produce
  a non-injective map into the irreducible hull when a socle element is
  supported on several generators. The coefficients are now chosen as a dual
  basis on each `ZF`-class of socle degrees (`_dual_basis_lambdas`), which
  makes the map injective; relevant relations are handled correctly.
- **Monomial matrices.** `monomial_matrix(i, res)` returns the differential
  `d^i` of an injective or irreducible resolution with its row and column
  labels (`MonomialMatrix`).
- **Bass numbers and minimality.** `graded_bass_numbers(M, p_F, i)` computes
  the graded Bass numbers at a face via `Ext`, and `is_minimal(res)` checks a
  computed injective resolution against them. `degrees_of_bass_numbers`
  returns the degrees of non-zero Bass numbers.
- **Injective hulls.** `injective_hull(M)` returns `E(M)` with the embedding.
- **Ideals.** `minimal_generating_set`, `number_of_generators`, `radical`,
  `intersect` of several ideals, `monoid_algebra_ideal` wrapper.
- **Getters** replacing field access: `injective_modules`, `cochain_maps`,
  `Q_graded_part`, `degree_shift`, `irreducible_sums`, `cochain_complex`,
  `is_exact`, `indecomposable_injectives`, `monoid_algebra`, `sectors`.

## Performance

- The shift needed to move all Bass numbers into `Q` is computed from a
  cheap over-approximation of their degrees read off a free resolution of
  the residue field (`compute_shift_bound`, default `shift = :bound`),
  instead of computing `Ext^j(k, M)` for every `j` (`shift = :helm_miller`).
  An LP/MILP-based minimal shift is available as `shift = :milp` and
  `shift = :milp_bound`.
- `degrees_of_bass_numbers` builds the residue field resolution once.
- Relevance checks for generators and relations in the coefficient
  computation precompute the polyhedra once per face; semigroup membership
  queries are cached per monoid algebra.
- Loop-invariant computations hoisted out of the main loops of
  `irreducible_hull` and `irreducible_resolution`.

## Bug fixes

- `local_cohomology_all`: a destructuring shadowed the loop counter, the
  `SectorPartitionLC` constructor was called with the wrong arguments, and the
  ideal-argument method compared a nonexistent field.
- `apply_gamma!`: crashed on an empty `hcat` when every summand was removed
  by `Γ_I`.
- `_compute_q_graded_part`: handles the empty case.
- `compute_shift`: no longer overshoots by one multiple of the ray sum.

## API cleanup

- Explicit imports instead of importing every OSCAR name; exports reduced to
  the public API; internal helpers unexported and named with a leading
  underscore.
- `check` keyword on `irreducible_resolution` and `injective_resolution`
  gates the internal exactness assertions.
- Printing follows the OSCAR conventions (capitalized first word, terse
  form for monoid algebras, ordinals in sector partition output).
- Removed `old_*` reference implementations and an unused combinatorial
  prototype.

## References

- [HM05] D. Helm, E. Miller, *Algorithms for graded injective resolutions and
  local cohomology over semigroup rings*. J. Symbolic Comput. 39 (2005),
  373–395.
