@doc raw"""
    torus_group(k::Field, m::Int)

Return the torus $(K^{\ast})^m$ defined over `k`, where $K$ is an algebraic closure of `k`.

!!! note
    In the context of computing invariant rings, there is no need to deal with the group structure of a torus: The torus $(K^{\ast})^m$ is specified by just giving $k$ and $m$.

# Examples
```jldoctest
julia> T = torus_group(QQ,2)
Torus of rank 2
  defined over QQ
```
"""
torus_group(k::Field, m::Int) = TorusGroup(k, m)

@doc raw"""
    rank(T::TorusGroup)

Return the rank of `T`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> rank(T)
2
```
"""
rank(G::TorusGroup) = G.rank

@doc raw"""
    base_ring(T::TorusGroup)

Return the field over which `T` is defined.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> base_ring(T)
Rational field
```
"""
base_ring(G::TorusGroup) = G.field

function Base.show(io::IO, ::MIME"text/plain", G::TorusGroup)
  io = pretty(io)
  println(io, "Torus of rank ", rank(G))
  print(terse(io), Indent(), "defined over ", Lowercase(), base_ring(G))
  print(io, Dedent())
end

function Base.show(io::IO, G::TorusGroup)
  if is_terse(io)
    print(io, "Torus")
  else
    io = pretty(io)
    print(io, "Torus of rank ", rank(G))
    print(terse(io), " defined over ", Lowercase(), base_ring(G))
  end
end

@doc raw"""
    representation_from_weights(T::TorusGroup, W::Union{ZZMatrix, Matrix{<:IntegerUnion}, Vector{<:IntegerUnion}})

Return the diagonal action of `T` with weights given by `W`.

The matrix `W` has one row for each coordinate of the vector space on which `T` acts, and
$m$ = `rank(T)` columns: if the $i$-th row is $(w_1, \dots, w_m)$, then
$t = (t_1, \dots, t_m) \in T$ multiplies the $i$-th coordinate by $t_1^{w_1} \cdots t_m^{w_m}$.
If $m = 1$, then `W` may also be a vector.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1])
Representation
  of torus of rank 2 defined over QQ
  with weights Vector{ZZRingElem}[[-1, 1], [-1, 1], [2, -2], [0, -1]]
```
"""
function representation_from_weights(
  G::TorusGroup, W::Union{ZZMatrix,Matrix{<:IntegerUnion},Vector{<:IntegerUnion}}
)
  n = rank(G)
  V = weights_from_matrix(n, W)
  return RepresentationTorusGroup(G, V)
end

function weights_from_matrix(
  n::Int, W::Union{ZZMatrix,Matrix{<:IntegerUnion},Vector{<:IntegerUnion}}
)
  V = Vector{Vector{ZZRingElem}}()
  if W isa Vector
    @req n == 1 "a vector of weights requires a torus of rank 1"
    for i in 1:length(W)
      push!(V, [ZZRingElem(W[i])])
    end
  else
    @req n == ncols(W) "the weight matrix must have rank(T) = $n columns and one row per coordinate"
    for i in 1:nrows(W)
      push!(V, [ZZRingElem(W[i, j]) for j in 1:ncols(W)])
    end
  end
  return V
end

@doc raw"""
    weights(r::RepresentationTorusGroup)

Return the weights of `r`, one for each coordinate of the vector space on which the torus acts diagonally.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> weights(r)
4-element Vector{Vector{ZZRingElem}}:
 [-1, 1]
 [-1, 1]
 [2, -2]
 [0, -1]
```
"""
weights(R::RepresentationTorusGroup) = R.weights

@doc raw"""
    group(r::RepresentationTorusGroup)

Return the torus group represented by `r`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> group(r)
Torus of rank 2
  defined over QQ
```
"""
group(R::RepresentationTorusGroup) = R.group

vector_space_dim(R::RepresentationTorusGroup) = length(weights(R))
degree(R::RepresentationTorusGroup) = vector_space_dim(R)

function Base.show(io::IO, ::MIME"text/plain", R::RepresentationTorusGroup)
  io = pretty(io)
  println(io, "Representation")
  println(io, Indent(), "of ", Lowercase(), group(R))
  print(io, "with weights ", weights(R))
  print(io, Dedent())
end

function Base.show(io::IO, R::RepresentationTorusGroup)
  if is_terse(io)
    print(io, "Representation of torus group")
  else
    print(io, "Representation of degree ", degree(R), " of torus group")
  end
end

@doc raw"""
    invariant_ring(r::RepresentationTorusGroup)

Return the invariant ring of the torus group represented by `r`.

!!! note
    The creation of invariant rings is lazy in the sense that no explicit computations are done until specifically invoked (for example, by the `fundamental_invariants` function).

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r)
Invariant ring
  of graded multivariate polynomial ring in 4 variables over QQ
  under group action of torus of rank 2 defined over QQ
```
"""
invariant_ring(R::RepresentationTorusGroup) = TorGroupInvarRing(R)

@doc raw"""
    fundamental_invariants(RT::TorGroupInvarRing)

Return a system of fundamental invariants for `RT`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r);

julia> fundamental_invariants(RT)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 X[1]^2*X[3]
 X[1]*X[2]*X[3]
 X[2]^2*X[3]
```
"""
function fundamental_invariants(z::TorGroupInvarRing)
  if !isdefined(z, :fundamental)
    R = z.representation
    z.fundamental = torus_invariants_fast(weights(R), polynomial_ring(z))
  end
  return copy(z.fundamental)
end

@doc raw"""
    reynolds_operator(RT::TorGroupInvarRing, f::MPolyRingElem)

Return the image of `f` under the Reynolds operator corresponding to `RT`.

This is the projection onto the invariant monomials: it keeps the terms of `f` whose
monomial $x^a$ has weight $\sum_i a_i w_i = 0$, where $w_i$ is the weight of $x_i$.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r);

julia> X = gens(polynomial_ring(RT));

julia> reynolds_operator(RT, (X[1] + X[2] + X[3] + X[4])^3)
3*X[1]^2*X[3] + 6*X[1]*X[2]*X[3] + 3*X[2]^2*X[3]
```
"""
function reynolds_operator(R::TorGroupInvarRing, elem::MPolyRingElem)
  W = weights(representation(R))
  @req ngens(parent(elem)) == length(W) "the number of variables must be the number of weights"
  has_weight_zero(a) = all(
    j -> is_zero(sum(a[i] * W[i][j] for i in eachindex(W); init=zero(ZZ))),
    1:rank(group(R)),
  )
  ctx = MPolyBuildCtx(parent(elem))
  for (c, a) in zip(coefficients(elem), exponents(elem))
    has_weight_zero(a) && push_term!(ctx, c, a)
  end
  return finish(ctx)
end

@doc raw"""
    null_cone_ideal(RT::TorGroupInvarRing)

Return Hilbert's ideal, the ideal of `polynomial_ring(RT)` generated by the invariants
of positive degree.

Its zero set is the null cone, the set of all $v \in V$ such that the closure of the
orbit $T v$ contains $0$. The generators are the fundamental invariants of `RT`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r);

julia> null_cone_ideal(RT)
Ideal generated by
  X[1]^2*X[3]
  X[1]*X[2]*X[3]
  X[2]^2*X[3]
```
"""
@attr ideal_type(PolyRingT) function null_cone_ideal(
  R::TorGroupInvarRing{FldT,PolyRingElemT,PolyRingT}
) where {FldT,PolyRingElemT,PolyRingT}
  # the ring is passed as there are no generators if only the constants are invariant
  return ideal(polynomial_ring(R), fundamental_invariants(R))
end

# Algorithm 4.3.1 from [DK15]. Computes Torus invariants without Reynolds operator.
function torus_invariants_fast(W::Vector{Vector{ZZRingElem}}, R::MPolyRing)
  # no check that length(W[i]) for all i is the same
  length(W) == ngens(R) || error(
    "number of weights must be equal to the number of generators of the polynomial ring"
  )
  n = length(W)
  r = length(W[1])
  #step 2
  if length(W[1]) == 1
    M = zero_matrix(ZZ, n, 1)
    for i in 1:n
      M[i, 1] = W[i][1]
    end
    C1 = lattice_points(convex_hull(M))
  else
    M = zero_matrix(ZZ, 2 * n, r)
    for i in 1:n
      M[i, 1:r] = 2 * r * W[i]
      M[n + i, 1:r] = -2 * r * W[i]
    end
    C1 = lattice_points(convex_hull(M))
  end

  #get a Vector{Vector{ZZRingElem}} from Vector{PontVector{ZZRingElem}}
  C = map(Vector{ZZRingElem}, C1)

  # 0 is not in C if r = 1 and all weights have the same sign
  index_0 = findfirst(is_zero, C)
  isnothing(index_0) && return elem_type(R)[]

  #step 3
  S = Vector{Vector{elem_type(R)}}()
  U = Vector{Vector{elem_type(R)}}()
  for point in C
    vars = [gen(R, i) for i in 1:n if point == W[i]]
    push!(S, vars)
    push!(U, copy(vars))
  end
  #step 4
  count = 0
  while true
    for j in 1:length(U)
      if length(U[j]) != 0
        m = U[j][1]
        w = C[j] #weight_of_monomial(m, W)
        #step 5 - 7
        for i in 1:n
          u = m * gen(R, i)
          v = w + W[i]
          if v in C
            index = findfirst(==(v), C)
            c = true
            for elem in S[index]
              if is_divisible_by(u, elem)
                c = false
                break
              end
            end
            if c == true
              push!(S[index], u)
              push!(U[index], u)
            end
          end
        end
        deleteat!(U[j], findall(==(m), U[j]))
      else
        count += 1
      end
    end
    if count == length(U)
      return S[index_0]
    else
      count = 0
    end
  end
end

@doc raw"""
    hilbert_series(RT::TorGroupInvarRing)

Return the Hilbert series of `RT` as `hilbert_series` does for the graded affine algebra
returned by `affine_algebra(RT)`: a pair $(p, q)$ of univariate polynomials such that
$p/q$ is the Hilbert series, where $q = (1-t^{d_1}) \cdots (1-t^{d_s})$ and
$d_1, \dots, d_s$ are the degrees of the fundamental invariants.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1]);

julia> RT = invariant_ring(r);

julia> hilbert_series(RT)
(-t^6 + 1, (-t^3 + 1)^3)
```
"""
hilbert_series(R::TorGroupInvarRing) = hilbert_series(affine_algebra(R)[1])
