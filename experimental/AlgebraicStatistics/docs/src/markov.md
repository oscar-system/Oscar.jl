```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Discrete random variables

The joint probability distribution of random variables ``X_1, \ldots, X_n``
is given by a tensor of order ``n``. If the random variable ``X_i`` takes
``d_i`` states, the tensor is of format ``d_1 \times \cdots \times d_n``
and consists of non-negative real numbers ``p_{x_1 \cdots x_n}``, for all
choices ``x_i \in [d_i]``, which sum to ``1``. The functions below deal
with the ambient polynomial ring in these ``p`` variables, special forms
in them like marginals, and conditional independence ideals.

```@docs
markov_ring(rvs::Pair{<:VarName, <:AbstractArray}...; unknown::VarName="p", K::Field=QQ, cached=false)
tensor_ring(rvs::Pair{<:VarName, <:AbstractArray}...; unknown::VarName="p", K::Field=QQ, cached=false)
ring(R::MarkovRing)
random_variables(R::MarkovRing)
ci_statements(R::MarkovRing)
unknowns(R::MarkovRing)
gens(R::MarkovRing)
find_random_variables(R::MarkovRing, K)
find_state(R::MarkovRing, K, x)
state_space(R::MarkovRing, K=random_variables(R))
marginal(R::MarkovRing, K, x)
ci_ideal(R::MarkovRing, stmts)
```
