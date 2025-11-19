# -*- Markov rings for discrete random variables -*-

struct MarkovRing
  ring::MPolyRing
  gens::GenDict{<:Tuple}
  state_spaces::Vector{AbstractArray}
end

@doc raw"""
    markov_ring(F::Field, states::Int...; unknown::VarName="p", cached=false)::MarkovRing
    markov_ring(states::Int...; unknown::VarName="p", cached=false)::MarkovRing

The polynomial ring whose unknowns are the entries of a probability tensor.
The argument `states` specifies the state space sizes of the random variables,
i.e., the dimensions of the tensor. The polynomial ring being constructed
will have one variable for each element in the cartesian product of
`1:s for s in states`. It is a multivariate polynomial ring whose variables are
named `p[...]` and whose coefficient field `F` is by default `QQ`.

If `cached` is `true`, the internally generated polynomial ring will be cached.

The name `tensor_ring` is an alias for the constructor `markov_ring` because
that is really what a `MarkovRing` is: the coordinate ring of tensors of a
fixed format. The name `MarkovRing` is kept for compatibility in terminology
with the Macaulay2 package `GraphicalModels`.

## Examples

```jldoctest
julia> R = markov_ring(2,2,2,2)
Markov ring over rational field for 4 random variables and states (2, 2, 2, 2)
```
"""
function markov_ring(F::Field, states::Int...; unknown::VarName="p", cached=false)
  @req length(states) >= 1 "need at least one random variable"
  @req all(s -> s >= 1, states) "all state spaces must have at least one element"
  random_variables = 1:length(states);
  state_spaces = [1:s for s in states];
  varindices = collect(Iterators.product(state_spaces...))
  varnames = [["$(unknown)[$(join(i, ", "))]" for i in varindices]...]
  R, p = polynomial_ring(F, varnames; cached=cached)
  d = GenDict{typeof(first(varindices))}(Dict([varindices[i] => p[i] for i in 1:length(varindices)]))
  return MarkovRing(R, d, state_spaces)
end

markov_ring(states::Int...; unknown::VarName="p", cached=false) =
  markov_ring(QQ, states...; unknown=unknown, cached=cached)

const tensor_ring = markov_ring

function Base.show(io::IO, R::MarkovRing)
  io = pretty(io)
  print(io, "Markov ring over ", Lowercase(), base_ring(R.ring), " for $(length(R.state_spaces)) random variables and states ",
    "(" * join([string(length(x)) for x in R.state_spaces], ", ") * ")",
  )
end

@doc raw"""
    ring(R::MarkovRing)

Return the multivariate polynomial ring inside `R`.

## Examples

```jldoctest
julia> R = markov_ring(2,2,2,2)
Markov ring over rational field for 4 random variables and states (2, 2, 2, 2)

julia> ring(R)
Multivariate polynomial ring in 16 variables p[1, 1, 1, 1], p[2, 1, 1, 1], p[1, 2, 1, 1], p[2, 2, 1, 1], ..., p[2, 2, 2, 2]
  over rational field
```
"""
function ring(R::MarkovRing)
  return R.ring
end

@doc raw"""
    gens(R::MarkovRing)

Return a dictionary for indexing the generators of `R` by their states.

## Examples

```jldoctest
julia> R = markov_ring(2,2,2,2)
Markov ring over rational field for 4 random variables and states (2, 2, 2, 2)

julia> gens(R)[1, 1, 2, 1]
p[1, 1, 2, 1]
```
"""
function gens(R::MarkovRing)
  return R.gens
end

@doc raw"""
    random_variables(R::MarkovRing)

Return the vector of `[1, ..., n]` where `n` is the number of random variables in the MarkovRing
(equivalently `n` is the order of the tensor).

## Examples

```jldoctest
julia> R = markov_ring(2, 2, 2, 2)
Markov ring over rational field for 4 random variables and states (2, 2, 2, 2)

julia> random_variables(R)
4-element Vector{Int64}:
 1
 2
 3
 4
```
"""
function random_variables(R::MarkovRing)
    return collect(1:length(R.state_spaces))
end

@doc raw"""
    ci_statements(R::MarkovRing)

Return all the `CIStmt` objects which can be formed on the `random_variables(R)`.

```jldoctest
julia> R = markov_ring(2, 2, 2, 2)
Markov ring over rational field for 4 random variables and states (2, 2, 2, 2)

julia> ci_statements(R)
24-element Vector{CIStmt}:
 [1 _||_ 2 | {}]
 [1 _||_ 2 | 3]
 [1 _||_ 2 | 4]
 [1 _||_ 2 | {3, 4}]
 [1 _||_ 3 | {}]
 [1 _||_ 3 | 2]
 [1 _||_ 3 | 4]
 [1 _||_ 3 | {2, 4}]
 [1 _||_ 4 | {}]
 [1 _||_ 4 | 2]
 ⋮
 [2 _||_ 3 | {1, 4}]
 [2 _||_ 4 | {}]
 [2 _||_ 4 | 1]
 [2 _||_ 4 | 3]
 [2 _||_ 4 | {1, 3}]
 [3 _||_ 4 | {}]
 [3 _||_ 4 | 1]
 [3 _||_ 4 | 2]
 [3 _||_ 4 | {1, 2}]
```
"""
ci_statements(R::MarkovRing) = ci_statements(random_variables(R))

@doc raw"""
    state_space(R::MarkovRing, K=random_variables(R))

Return all states that the random subvector indexed by `K` can attain
in the ring `R`. The result is an `Iterators.product` iterator unless
`K` has only one element in which case it is a vector.

## Examples

```jldoctest
julia> R = markov_ring(2, 3, 2, 4)
Markov ring over rational field for 4 random variables and states (2, 3, 2, 4)

julia> collect(state_space(R, [1, 4]))
2×4 Matrix{Tuple{Int64, Int64}}:
 (1, 1)  (1, 2)  (1, 3)  (1, 4)
 (2, 1)  (2, 2)  (2, 3)  (2, 4)
```
"""
function state_space(R::MarkovRing, K=random_variables(R))
  return length(K) == 1 ? R.state_spaces[K[1]] : Iterators.product([R.state_spaces[i] for i in K]...)
end

# -*- CI equations for MarkovRing -*-

# p is a permutation of 1:n and x is a vector of length n.
# This method returns the components of x permuted by p,
# e.g. apply_permutation([3,2,1], [0,1,2]) == [2,1,0].
function apply_permutation(p, x)
  px = Vector{Int64}(undef, length(p))
  for i in p
    px[p[i]] = x[i]
  end
  return px
end

@doc raw"""
    marginal(R::MarkovRing, K, x)

Return a marginal as a sum of unknowns from `R`. The argument `K` lists
random variables which are fixed to the event `x`; all other random
variables in `R` are summed over their respective state spaces.

## Examples

```jldoctest
julia> R = markov_ring(2, 3, 2, 4)
Markov ring over rational field for 4 random variables and states (2, 3, 2, 4)

julia> marginal(R, [1, 4], [2, 1]) # sum all p[i,j,k,k] where i=2 and l=1
p[2, 1, 1, 1] + p[2, 2, 1, 1] + p[2, 3, 1, 1] + p[2, 1, 2, 1] + p[2, 2, 2, 1] + p[2, 3, 2, 1]
```
"""
function marginal(R::MarkovRing, K, x)
  p = R.gens
  N = random_variables(R)
  M = setdiff(N, K)
  summands = Vector()
  for y in state_space(R, M)
    idx = apply_permutation(
      [K..., M...],
      [x..., y...]
    )
    push!(summands, p[idx...])
  end
  return length(summands) == 1 ? summands[1] : +(summands...)
end

@doc raw"""
    ci_ideal(R::MarkovRing, stmts)::MPolyIdeal

Return the ideal for the conditional independence statements given by `stmts`.

## Examples

```jldoctest
julia> R = markov_ring(2, 2, 2)
Markov ring over rational field for 3 random variables and states (2, 2, 2)

julia> ci_ideal(R, [CI"1,3|2", CI"2,3|1"])
Ideal generated by
  p[1, 1, 1]*p[2, 1, 2] - p[2, 1, 1]*p[1, 1, 2]
  p[1, 2, 1]*p[2, 2, 2] - p[2, 2, 1]*p[1, 2, 2]
  p[1, 1, 1]*p[1, 2, 2] - p[1, 2, 1]*p[1, 1, 2]
  p[2, 1, 1]*p[2, 2, 2] - p[2, 2, 1]*p[2, 1, 2]
```
"""
function ci_ideal(R::MarkovRing, stmts)
  eqs = Vector()
  for stmt in stmts
    # The proper CI equations are the familiar 2x2 determinants in marginals of the probability tensor.
    QI = state_space(R, stmt.I)
    QJ = state_space(R, stmt.J)
    IJK = [stmt.I..., stmt.J..., stmt.K...]
    for z in state_space(R, stmt.K), x in subsets(collect(QI), 2), y in subsets(collect(QJ), 2)
      p11 = marginal(R, IJK, [x[1], y[1], z...]);
      p12 = marginal(R, IJK, [x[1], y[2], z...]);
      p21 = marginal(R, IJK, [x[2], y[1], z...]);
      p22 = marginal(R, IJK, [x[2], y[2], z...]);
      push!(eqs, p11*p22 - p12*p21)
    end
  end
  return ideal(eqs)
end
