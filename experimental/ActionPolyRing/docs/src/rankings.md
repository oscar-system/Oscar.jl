```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Rankings](@id actionpolyranking)

Let ``A`` be an action polynomial ring with ``m`` elementary symbols and ``n`` commuting action maps.  
Define ``\underline{m} := \{1, \ldots, m\}``, with ``m \in \mathbb{N}``.  

Many algorithms require the jet variables of ``A`` to be totally ordered. This is achieved using *rankings*.  
Analogous to how monomial orderings in a multivariate polynomial ring in ``n`` variables correspond to total
orderings of ``\mathbb{N}_0^n``, rankings of jet variables correspond to total orderings of ``X := \underline{m} \times \mathbb{N}_0^n``.  

Specifically, we identify a jet variable ``(u_i)_J \in A`` with the pair ``(i, J) \in X``.

## [Riquier rankings](@id actionpolyriquierranking)

The rankings we use are called *Riquier rankings*. By definition, these are rankings of ``X`` that extend to a ranking of
``\{1\} \times \mathbb{N}^{m+n}``.  

Equivalently, there exists a positive integer ``s`` and an ``s \times (m+n)`` real matrix ``M`` such that the total ordering
of the jet variables defined by ``X`` coincides with the ordering obtained from the [matrix ordering](@ref "Matrix Orderings")
on ``\mathbb{N}_0^{m+n}`` defined by ``M``.

For this construction, we identify a jet variable ``(u_i)_J \in A`` with ``(e_i, J) \in \mathbb{N}_0^{m+n}``,
where ``e_i`` is the ``i``-th unit row and restrict ourselves to integer matrices ``M``. In this context, we call ``M``
a *Riquier matrix*.

!!! note
    Not all Riquier rankings are obtained from integral Riquier matrices. However, this is the case if we only require a total
    ordering of a finite subset of ``X``. Thus, only considering integer matrices is sufficient for practical use.

### Constructing Riquier rankings

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
