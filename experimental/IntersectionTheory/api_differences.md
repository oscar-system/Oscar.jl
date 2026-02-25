# API Differences: Oscar IntersectionTheory vs jieaosong/IntersectionTheory

This document catalogs the differences between Oscar's `experimental/IntersectionTheory`
and the standalone [jieaosong/IntersectionTheory](https://github.com/jieaosong/IntersectionTheory)
package. Both originate from the same codebase but have diverged significantly.

---

## 1. Architecture

| Aspect | Oscar | Standalone |
|--------|-------|------------|
| **Package type** | Experimental module within Oscar.jl | Independent Julia package |
| **Ring system** | Oscar's `MPolyDecRing` / `MPolyQuoRing` (graded polynomial rings) | Custom `ChRing` / `ChRingElem` wrapping Singular's `PolyRing` |
| **Attributes** | Oscar's `@attributes` macro, `get/set_attribute!` | Same (migrated from `@declare_other` / `get/set_special`) |
| **Coefficient rings** | Oscar's ring types (`QQFieldElem`, `ZZRingElem`, etc.) | Singular/Nemo types (same names after upstream modernization) |
| **Dependencies** | Full Oscar stack (Nemo, Hecke, Singular, Polymake, GAP, ...) | Singular.jl, Nemo.jl, GAP.jl only |

---

## 2. Naming Conventions

Oscar follows the [Oscar naming conventions](https://docs.oscar-system.org/stable/DeveloperDocumentation/styleguide/)
with `snake_case` and more descriptive names. The standalone package uses shorter names.

### 2.1 Constructor Functions

| Oscar | Standalone | Notes |
|-------|------------|-------|
| `abstract_projective_space(n)` | `proj(n)` | |
| `abstract_grassmannian(k, n)` | `grassmannian(k, n)` | Oscar also has alias `grassmannian` |
| `abstract_flag_variety(dims...)` | `flag(dims...)` | Oscar also has alias `flag` |
| `abstract_point()` | `point()` | |
| `abstract_variety(n, R)` | `variety(n, R)` | |
| `abstract_variety(n, symbols, degs)` | `variety(n, symbols, degs)` | |
| `abstract_variety(n, bundles)` | `variety(n, bundles)` | |
| `abstract_variety(n)` | `variety(n)` | |
| `abstract_bundle(X, ch)` | `bundle(X, ch)` | |
| `abstract_bundle(X, r, c)` | `bundle(X, r, c)` | |
| `abstract_hirzebruch_surface(n)` | — | Not in standalone |
| `abstract_curve(g)` | `curve(g)` | **Oscar does not have this** |
| — | `K3(g)` | **Not in Oscar** |
| — | `quadric(n)` | **Not in Oscar** |
| — | `cayley_plane()` | **Not in Oscar** |
| — | `cayley_grassmannian()` | **Not in Oscar** |
| `projective_bundle(F)` | `proj(F)` | Standalone reuses `proj` |
| `flag_bundle(F, ranks)` | `flag(ranks, F)` | |

### 2.2 Getter Functions

| Oscar | Standalone | Notes |
|-------|------------|-------|
| `tautological_bundles(X)` | `bundles(X)` | |
| `structure_map(X)` | `X.struct_map` | Oscar has both getter and field `X.structure_map` |
| `polarization(X)` | `X.O1` | Oscar has getter function; standalone uses field directly |
| `point_class(X)` | `X.point` | Oscar has getter function |
| `chow_ring(X)` | `X.ring` | Oscar has getter function; both have `X.ring` field |
| `tangent_bundle(X)` | `tangent_bundle(X)` | Same |
| `cotangent_bundle(X)` | `cotangent_bundle(X)` | Same |
| `canonical_bundle(X)` | `canonical_bundle(X)` | Same |
| `canonical_class(X)` | `canonical_class(X)` | Same |
| `base(X)` | `base_ring(X)` | Different function names |

### 2.3 Characteristic Classes

| Oscar | Standalone | Notes |
|-------|------------|-------|
| `chern_character(F)` | `ch(F)` | |
| `total_chern_class(F)` | `chern(F)` | |
| `chern_class(F, k)` | `chern(k, F)` | **Argument order differs** |
| `top_chern_class(F)` | `ctop(F)` | |
| `total_segre_class(F)` | `segre(F)` | |
| `segre_class(F, k)` | `segre(k, F)` | **Argument order differs** |
| `todd_class(F)` | `todd(F)` | |
| `todd_class(X)` | `todd(X)` | Same pattern |
| `euler_number(X)` | `euler(X)` | |
| `euler_characteristic(F)` | `chi(F)` | |
| `euler_pairing(F, G)` | `chi(F, G)` | |
| `total_pontryagin_class(F)` | `pontryagin(F)` | |
| `pontryagin_class(F, k)` | `pontryagin(k, F)` | **Argument order differs** |

### 2.4 Bundle Operations

| Oscar | Standalone | Notes |
|-------|------------|-------|
| `symmetric_power(F, k)` | `symmetric_power(k, F)` | **Argument order differs** |
| `exterior_power(F, k)` | `exterior_power(k, F)` | **Argument order differs** |
| `schur_functor(F, λ)` | `schur_functor(λ, F)` | **Argument order differs** |
| `dual(F)` | `dual(F)` | Same |
| `det(F)` | `det(F)` | Same |

### 2.5 Morphism Functions

| Oscar | Standalone | Notes |
|-------|------------|-------|
| `map(X, Y, pullbacks)` | `hom(X, Y, pullbacks)` | Different name |
| `pullback(f, x)` | `pullback(f, x)` | Same |
| `pushforward(f, x)` | `pushforward(f, x)` | Same (Oscar has namespace issues, needs `IntersectionTheory.pushforward`) |
| `identity_map(X)` | `identity_hom(X)` | Different name |
| `compose(f, g)` | `f * g` | Oscar has function; standalone uses `*` operator |

### 2.6 Blow-up Functions

| Oscar | Standalone | Notes |
|-------|------------|-------|
| `blow_up(i)` → `(Bl, E, j)` | `blowup(i)` → `(Bl, PN)` | **Different return value!** |
| `blow_up_points(X, n)` | `blowup_points(n, X)` | Different name and argument order |
| `extend_inclusion(i)` | `_inclusion(i)` | Oscar exports it; standalone is internal |

**Important:** The blow-up return value differs significantly:
- **Oscar:** Returns `(Bl, E, j)` where `E` is the exceptional divisor variety,
  `j` is the inclusion `E → Bl`.
- **Standalone:** Returns `(Bl, PN)` where `PN` is the exceptional divisor variety.
  The maps `j: PN → Bl` and `g: PN → center` are stored as
  `get_attribute(PN, :projections)`.

### 2.7 Other Functions

| Oscar | Standalone | Notes |
|-------|------------|-------|
| `zero_locus_section(F)` | `section_zero_locus(F)` | Different name |
| `complete_intersection(X, degs)` | `complete_intersection(X, degs)` | Same |
| `degeneracy_locus(k, F, G)` | `degeneracy_locus(k, F, G)` | Same |
| `schubert_class(G, λ)` | `schubert_class(G, λ)` | Same |
| `schubert_classes(m, G)` | `schubert_classes(m, G)` | Same |
| `trim!(X)` | `trim!(X)` | Same |
| `pair_trim!(X)` | `pair_trim!(X)` | Same |
| `simplify(x)` | `simplify(x)` | Same |

---

## 3. Type Names

| Oscar | Standalone |
|-------|------------|
| `AbstractVariety` | `AbsVariety` |
| `AbstractBundle` | `AbsBundle` |
| `AbstractVarietyMap` | `AbsVarietyHom` |
| `MPolyDecRingOrQuo` | `ChRing` |
| `MPolyDecRingOrQuoElem` | `ChRingElem` |
| `TnVariety` | `TnVariety` |
| `TnBundle` | `TnBundle` |
| `TnRep` | `TnRep` |

---

## 4. Struct Field Names

### AbstractVariety / AbsVariety

| Oscar field | Standalone field | Type (Oscar) | Type (Standalone) |
|-------------|------------------|--------------|-------------------|
| `dim` | `dim` | `Int` | `Int` |
| `ring` | `ring` | `MPolyDecRingOrQuo` | `ChRing` |
| `base` | `base` | `Ring` | `Ring` |
| `point` | `point` | `MPolyDecRingOrQuoElem` | `ChRingElem` |
| `O1` | `O1` | `MPolyDecRingOrQuoElem` | `ChRingElem` |
| `T` | `T` | `AbstractBundle` | `AbsBundle` |
| `bundles` | `bundles` | `Vector{AbstractBundle}` | `Vector{AbsBundle}` |
| `structure_map` | `struct_map` | `AbstractVarietyMap` | `AbsVarietyHom` |

### AbstractBundle / AbsBundle

| Oscar field | Standalone field |
|-------------|------------------|
| `parent` | `parent` |
| `rank` | `rank` |
| `ch` | `ch` |
| `chern` | `chern` |

### AbstractVarietyMap / AbsVarietyHom

| Oscar field | Standalone field |
|-------------|------------------|
| `domain` | `domain` |
| `codomain` | `codomain` |
| `dim` | `dim` |
| `pullback` | `pullback` (type differs: `AffAlgHom` vs `ChAlgHom`) |
| `pushforward` | `pushforward` (type differs: `MapFromFunc` vs `FunctionalMap`) |
| `O1` | `O1` |
| `T` | `T` |

---

## 5. Features Only in Oscar

1. **`abstract_hirzebruch_surface(n)`** — Convenience constructor for Hirzebruch surfaces.
2. **`kontsevich_moduli_space(d, n)`** — Moduli of stable maps $\overline{M}_{0,n}(\mathbb{P}^r, d)$.
3. **`gromov_witten_invariant(d, degs)`** — Gromov-Witten invariants via localization.
4. **`instanton_number(d, degs)`** — Instanton numbers via mirror symmetry.
5. **`lines_on_hypersurface` / `linear_subspaces_on_hypersurface`** — Counting linear subspaces on hypersurfaces via Bott localization.
6. **`hilbert_polynomial(F)` / `hilbert_polynomial(X)`** — Hilbert polynomial computation.
7. **`a_hat_genus(k, X)` / `a_hat_genus(X)`** — $\hat{A}$-genus (also in standalone under same name).
8. **`l_genus(k, X)` / `l_genus(X)`** — $L$-genus.
9. **`signature(X)`** — Hirzebruch signature.
10. **`chern_number(X, λ)` / `chern_numbers(X)`** — Chern number computation.
11. **`betti_numbers(X)` / `basis(X)`** — Explicit Chow ring basis.
12. **`intersection_matrix(X)` / `dual_basis(X)`** — Intersection-theoretic data.
13. **`product(X, Y)`** — Explicit product of varieties.
14. **`graph(f)`** — Graph of a morphism.
15. **`identity_map(X)`** — Identity morphism.
16. **`compose(f, g)`** — Explicit composition.

---

## 6. Features Only in Standalone

1. **`curve(g)`** — Construct a curve of genus $g$ (accepts `RingElement` or `String`).
2. **`K3(g)`** — Construct a K3 surface of genus $g$.
3. **`quadric(n)` / `quadric(P)`** — Quadric hypersurface with spinor bundles.
4. **`cayley_plane()`** — Cayley plane $\mathbf{OP}^2$.
5. **`cayley_grassmannian()`** — Cayley Grassmannian $\mathbf{CG}$.
6. **`Moduli.jl`** — Moduli of matrices and twisted cubics (`matrix_moduli`, `twisted_cubics`).
7. **`Weyl.jl`** — Weyl group actions on flag varieties (`homogeneous_variety`).
8. **`Cobord.jl`** — Cobordism ring and related computations.
9. **Threaded `_integral`** — Multi-threaded Bott formula evaluation using `Threads.@spawn`.
10. **`add_rels!(R, rels)`** — Add relations to an existing `ChRing`.
11. **`parameters(vars...)`** — Create symbolic parameters from a function field.
12. **`libgober_wood_polynomial`** — Libgober-Wood polynomial computation.
13. **`milnor(X)`** — Milnor number.
14. **`homogeneous_variety(G)`** — Construct homogeneous varieties of exceptional Lie types.

Note: Items 6–8 are **commented out** in Oscar's includes (`Moduli.jl`, `Weyl.jl`),
and `Cobord.jl` was never ported.

---

## 7. Argument Order Summary

The most impactful difference between the two implementations is **argument order**.
Oscar puts the "container" (bundle/variety) first; the standalone package puts
numerical/partition arguments first. This follows different conventions:
Oscar follows Julia/Oscar style (object first), while the standalone package
follows the mathematical notation (operation parameter first).

| Operation | Oscar | Standalone |
|-----------|-------|------------|
| Chern class | `chern_class(F, k)` | `chern(k, F)` |
| Segre class | `segre_class(F, k)` | `segre(k, F)` |
| Pontryagin class | `pontryagin_class(F, k)` | `pontryagin(k, F)` |
| Symmetric power | `symmetric_power(F, k)` | `symmetric_power(k, F)` |
| Exterior power | `exterior_power(F, k)` | `exterior_power(k, F)` |
| Schur functor | `schur_functor(F, λ)` | `schur_functor(λ, F)` |
| Schubert class | `schubert_class(G, λ)` | `schubert_class(G, λ)` |

Note: `schubert_class` has the **same** argument order in both.

---

## 8. Internal Differences

### Ring Infrastructure

- **Oscar** uses `MPolyDecRing` (graded multivariate polynomial ring from Oscar/AbstractAlgebra)
  with `MPolyQuoRing` for quotients. The truncation/grading is managed by Oscar's ring system.
- **Standalone** uses a custom `ChRing` type wrapping Singular's `PolyRing` with manual
  weight tracking (`w::Vector{Int}`) and an optional `sideal` for the ideal.
  Truncation is stored as an attribute (`:truncate`).

### Morphism Construction

- **Oscar** uses `hom(R, S, images)` from Oscar to create `AffAlgHom`.
- **Standalone** uses `Singular.AlgebraHomomorphism` directly.

### Bott Localization

- **Oscar** has a sequential implementation of `integral(c::TnBundleChern)`.
  No batch `_integral` function, no threading.
- **Standalone** has both `integral(c::TnBundleChern)` (sequential) AND
  a batch `_integral(F, pp)` function with threaded computation using
  `Threads.@spawn` and pre-allocated arrays for maximum performance.

### Blow-up

- **Oscar** constructs the blow-up ring using `graded_polynomial_ring` and
  `present_finite_extension_ring` from Oscar's algebra infrastructure.
  Returns `(Bl, E, j)` with the inclusion map.
- **Standalone** constructs the blow-up ring using Singular's `PolynomialRing`
  and `_pushfwd` for the finite extension. Returns `(Bl, PN)` with
  projections stored as attributes.
