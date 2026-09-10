# Changes — experimental/InjectiveResolutions

Working summary of updates to the `InjectiveResolutions` experimental package, to accompany the upcoming pull request.

The package provides tools for computing injective resolutions and local cohomology of finitely generated modules over monoid algebras, implementing algorithms from [HM05]. This update focuses on performance and correctness fixes in `LocalCohomology.jl`.

---

## Performance

### `degrees_of_bass_numbers` — share the residue-field resolution
Previously, each cohomological degree `j ∈ 0:i` triggered a fresh `ext(k, M, j)`, which rebuilt the free resolution of the residue field `k = R_Q/m` from scratch. The resolution is now built **once** (to length `i+2`), `Hom(-, M)` is applied **once**, and each `Ext^j(k, M)` is obtained as a homology lookup on the resulting cochain complex.

### New: `degrees_of_bass_numbers_bound(M, i)`
A cheap over-approximation that avoids the `Hom(-, M)` and `homology` calls entirely:

$$
D := \{\deg(g) - \deg(e) : g \in \mathrm{gens}(M),\ e \in \mathrm{gens}(F_j),\ 0 \le j \le i+1\}.
$$

Since `Ext^j(k, M)` is a subquotient of `Hom(F_j, M) = ⊕_e M(\deg e)` and `supp(M) ⊆ ⋃_g (\deg(g) + Q)`, the Bass-number degrees lie in `D + Q`. Because `Q` is a semigroup, any shift `a` with `D + a ⊆ Q` automatically satisfies `(bass degrees) + a ⊆ Q`. The bound is only used as a sufficient set for `compute_shift_bound`; the exact `degrees_of_bass_numbers` is unchanged.

This is a large speedup on **quotient-ring** monoid algebras, where the residue-field resolution has fast-growing ranks and the `Hom(F_j, M) ≅ M^{n_j}` term dominates the runtime.

### New: `compute_shift_bound(M, i)`
Shift computation driven by `degrees_of_bass_numbers_bound`. Returns a valid (possibly larger) shift; used by `injective_resolution` by default.

### `compute_shift` — per-point translation
Rewritten so each Bass-number degree is translated by `c` only until it lands in `Q`, rather than re-testing the entire list on every iteration. Returns the minimal `j·c` (the previous version always over-shot by one `c`).

### `injective_resolution`
Now uses `compute_shift_bound`. The previous behaviour is preserved as `old_injective_resolution`.

### `_coefficients_normal` (split out of `coefficients`)
The per-face relevance checks for generators and relations involve polyhedral intersection of Minkowski sums. The polyhedra `\deg(g) + cone(Q)` (resp. `\deg(r) + cone(Q)`) depend only on the face, not on the socle basis element `b`, so they are now precomputed once per face and reused across all `b ∈ Bp`. Eliminates `|Bp| × (|gens(N)| + |relations(N)|)` redundant `convex_hull` / Minkowski-sum constructions per face.

### Misc. loop-invariant hoisting
- `irreducible_hull`: `(ideal(kQ, []) * Mi)[1]` hoisted out of the face loop.
- `irreducible_resolution`: `monomial_basis(R_Q, degree(Mi[ii]))[1]` hoisted out of the inner column loop.

---

## Bug fixes — `LocalCohomology.jl`

### `local_cohomology_all`: shadowed loop variable
The destructure `J, phi, psi, (j, k) = apply_gamma!(...)` was clobbering the outer loop counter `j`, silently shifting its meaning mid-iteration. Renamed the destructured names to `j0, k0`.

### `local_cohomology_all`: wrong `SectorPartitionLC` constructor call
The call `SectorPartitionLC(M, j, I.ideal, Hj, maps_needed(kQ, Hj))` passed 5 arguments to a 3-argument constructor, and used `I.ideal` (an `Ideal`) where a `MonoidAlgebraIdeal` is expected. Fixed: construct via the 3-arg constructor, then set `.sectors` and `.maps` separately.

### `apply_gamma!`: empty `hcat` crash
`transpose(hcat(rows_phi...))` (and the analogous `rows_psi` line) crashed when every indecomposable injective in `J^0` (resp. `J^1`) was filtered out by `Γ_I` — `hcat()` on an empty vector raised an error. Guarded both sites.

---

## API additions

Reference implementations kept for correctness/performance comparison alongside the new versions:

- `old_degrees_of_bass_numbers` — the original `ext`-loop implementation
- `old_compute_shift` — the original re-translate-all-degrees loop
- `old_injective_resolution` — uses the exact `compute_shift`

New public functions:

- `degrees_of_bass_numbers_bound`
- `compute_shift_bound`

---

## Other additions

<!-- Anna: fill in details of the items below that you added independently of the
     perf/bug-fix work. Each was already in your working tree by the time this
     document was started, so I left them as stubs rather than guessing. -->

- `saturation_ideal`, `saturation_map` — TODO
- `holes_module` — TODO
- `is_Q_graded` — TODO
- `injective_hull` — TODO
- `compute_Q_graded_part`: now handles the empty `indec_injectives` case (returns the zero subquotient of `graded_free_module(kQ, 0)`).

---

## Known issues / not addressed in this update

- `irreducible_resolution` produces a non-exact complex for multi-generator (non-cyclic) modules. This is a **pre-existing** failure (already documented with `#this test fails` in `test/injective_res.jl` before this update). Diagnosis ongoing; localized to the multi-generator branch of `coefficients` / map assembly (assertion `is_homogeneous(_c_b[i])` holds, so the all-ones-evaluation step is not the cause; shape of `_lambda` matches `ngens(Mi)`).
- The two `@assert is_injective(fi)` / `@assert is_welldefined(fi)` checks in `irreducible_resolution`'s main loop are still always-on. Each runs a kernel computation per iteration; gating them behind a debug flag is an outstanding minor speedup, left for a follow-up to keep this PR focused on already-validated changes.

---

## References

- [HM05] D. Helm, E. Miller, *Algorithms for graded injective resolutions and local cohomology over semigroup rings*. J. Symbolic Comput. 39 (2005), 373–395.
