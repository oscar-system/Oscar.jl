@doc raw"""
    torus_group(K::Field, m::Int)

Return the torus $(K^{\ast})^m$.

!!! note
    In the context of computing invariant rings, there is no need to deal with the group structure of a torus: The torus $(K^{\ast})^m$ is specified by just giving $K$ and $m$.
# Examples
```jldoctest
julia> T = torus_group(QQ,2)
Torus of rank 2
  over QQ
```
"""
torus_group(F::Field, n::Int) = TorusGroup(F, n)

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
  print(terse(io), Indent(), "over ", Lowercase(), base_ring(G))
  print(io, Dedent())
end

function Base.show(io::IO, G::TorusGroup)
  if is_terse(io)
    print(io, "Torus")
  else
    io = pretty(io)
    print(io, "Torus of rank ", rank(G))
    print(terse(io), " over ", Lowercase(), base_ring(G))
  end
end

@doc raw"""
    representation_from_weights(T::TorusGroup, W::Union{ZZMatrix, Matrix{<:Integer}, Vector{<:Int}})

Return the diagonal action of `T` with weights given by `W`.

# Examples
```jldoctest
julia> T = torus_group(QQ,2);

julia> r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1])
Representation
  of torus of rank 2 over QQ
  with weights Vector{ZZRingElem}[[-1, 1], [-1, 1], [2, -2], [0, -1]]
```
"""
function representation_from_weights(
  G::TorusGroup, W::Union{ZZMatrix,Matrix{<:Integer},Vector{<:Int}}
)
  n = rank(G)
  V = weights_from_matrix(n, W)
  return RepresentationTorusGroup(G, V)
end

function weights_from_matrix(n::Int, W::Union{ZZMatrix,Matrix{<:Integer},Vector{<:Int}})
  V = Vector{Vector{ZZRingElem}}()
  if W isa Vector
    n == 1 || error("Incompatible weights")
    for i in 1:length(W)
      push!(V, [ZZRingElem(W[i])])
    end
  else
    n == ncols(W) || error("Incompatible weights")
    #assume columns = G.group[2]
    for i in 1:nrows(W)
      push!(V, [ZZRingElem(W[i, j]) for j in 1:ncols(W)])
    end
  end
  return V
end

@doc raw"""
    weights(r::RepresentationTorusGroup)

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
  over QQ
```
"""
group(R::RepresentationTorusGroup) = R.group

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
    # I don't know what else to print here
    print(io, "Representation of torus group")
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
  under group action of torus of rank 2 over QQ
```
"""
invariant_ring(R::RepresentationTorusGroup) = TorGroupInvarRing(R)

polynomial_ring(R::TorGroupInvarRing) = R.poly_ring
group(R::TorGroupInvarRing) = R.group
representation(R::TorGroupInvarRing) = R.representation

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

function Base.show(io::IO, ::MIME"text/plain", R::TorGroupInvarRing)
  io = pretty(io)
  println(io, "Invariant ring")
  println(io, Indent(), "of ", Lowercase(), polynomial_ring(R))
  print(io, "under group action of ", Lowercase(), group(R))
  print(io, Dedent())
end

function Base.show(io::IO, R::TorGroupInvarRing)
  if is_terse(io)
    print(io, "Invariant ring")
  else
    io = pretty(io)
    print(io, "Invariant ring of ")
    print(terse(io), Lowercase(), group(R))
  end
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
    affine_algebra(RT::TorGroupInvarRing)

Return the invariant ring `RT` as an affine algebra (this amounts to compute the algebra syzygies among the fundamental invariants of `RT`).

In addition, if `A` is this algebra, and `R` is the polynomial ring of which `RT` is a subalgebra,
return the inclusion homomorphism  `A` $\hookrightarrow$ `R` whose image is `RT`.

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

julia> affine_algebra(RT)
(Quotient of multivariate polynomial ring by ideal (-t[1]*t[3] + t[2]^2), Hom: quotient of multivariate polynomial ring -> graded multivariate polynomial ring)
```
"""
function affine_algebra(R::TorGroupInvarRing)
  if !isdefined(R, :presentation)
    V = fundamental_invariants(R)
    s = length(V)
    weights_ = zeros(Int, s)
    for i in 1:s
      weights_[i] = total_degree(V[i])
    end
    S, _ = graded_polynomial_ring(
      base_ring(group(representation(R))), :t => 1:s; weights=weights_, cached=false
    )
    R_ = polynomial_ring(R)
    StoR = hom(S, R_, V)
    I = kernel(StoR)
    Q, StoQ = quo(S, I)
    QtoR = hom(Q, R_, V)
    R.presentation = QtoR
  end
  return domain(R.presentation), R.presentation
end
