```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Riquier rankings in OSCAR](@id riquier_rankings_in_OSCAR)

In OSCAR, we define rankings, i.e., total orderings of ``X``, by combining the natural *less-than* relation on ``\underline{m}``
with a customizable total ordering on ``\mathbb{N}_0^n``. The latter is constructed as a matrix ordering; see
[`index_ordering_matrix`](@ref).

The way of combining of these two total orderings to obtain a total ordering of ``X`` is specified by an ordered partition of
``\underline{m}``, i.e. by grouping the elements of ``\underline{m}`` into blocks. The first block is considered largest and so on.
See [`partition`](@ref partition(r::ActionPolyRingRanking)).

Consider two elements ``x_1 = (i_1, J_1), x_2 = (i_2, J_2) \in X``. Then ``x_1 < x_2`` if and only if:
1. The block of ``i_1`` is smaller than the one of ``i_2``.
2. ``i_1`` and ``i_2`` are in the same block and ``J_1 < J_2`` with respect to the total ordering on ``\mathbb{N}_0^n``.
3. ``i_1`` and ``i_2`` are in the same block, ``J_1 = J_2`` and ``i_1 < i_2``.

---

Each action polynomial ring has an internal field `ranking`, which can be modified using the [`set_ranking!`](@ref)-method.

```@docs
set_ranking!
```

!!! note
    [`set_ranking!`](@ref) is also called upon construction with [`difference_polynomial_ring`](@ref) or
    [`differential_polynomial_ring`](@ref). These constructors allow for the same keywords.

---

```@docs
parent(r::ActionPolyRingRanking)
ranking(r::DifferencePolyRing)
riquier_matrix(r::ActionPolyRingRanking)
index_ordering_matrix(r::ActionPolyRingRanking)
partition(r::ActionPolyRingRanking)
```
