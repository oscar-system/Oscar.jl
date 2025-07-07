
@attributes mutable struct MPolyDecRing{T, S} <: AbstractAlgebra.MPolyRing{T}
  R::S
  D::FinGenAbGroup
  d::Vector{FinGenAbGroupElem}
  lt::Any
  hilbert_series_parent::Generic.LaurentPolyWrapRing{ZZRingElem, ZZPolyRing}
  multi_hilbert_series_parent::Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}

  function MPolyDecRing(R::S, d::Vector{FinGenAbGroupElem}) where {S}
    @req !(R isa MPolyDecRing) "cannot graded polynomial ring which is already decorated"
    @assert length(d) == ngens(R)
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    return r
  end
  function MPolyDecRing(R::S, d::Vector{FinGenAbGroupElem}, lt) where {S}
    @req !(R isa MPolyDecRing) "cannot filter polynomial ring which is already decorated"
    @assert length(d) == ngens(R)
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    r.lt = lt
    return r
  end
end

generator_degrees(S::MPolyDecRing) = S.d

@doc raw"""
    grading_group(R::MPolyDecRing)

If `R` is, say, `G`-graded, then return `G`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> G = grading_group(R)
Z

julia> H = free_abelian_group(1)
Z

julia> G == H
false

```
"""
grading_group(R::MPolyDecRing) = R.D

function is_graded(R::MPolyDecRing)
   return !isdefined(R, :lt)
end

_grading(R::MPolyDecRing) = R.d

@doc raw"""
    is_graded(R::MPolyRing)

Return `true` if `R` is graded, `false` otherwise.
"""
is_graded(R::MPolyRing) = false

is_filtered(R::MPolyDecRing) = isdefined(R, :lt)
is_filtered(::MPolyRing) = false

function show(io::IO, ::MIME"text/plain", W::MPolyDecRing)
  AbstractAlgebra.@show_name(io, W)
  AbstractAlgebra.@show_special(io, W)
  io = pretty(io)
  R = forget_decoration(W)
  print(io, R)
  if is_trivial(grading_group(W))
    print(io, " graded by the trivial group")
  else
    if is_filtered(W)
      println(io, " filtrated by")
    else
      println(io, " graded by")
    end
    g = gens(R)
    print(io, Indent())
    for i = 1:ngens(R)
      if i == ngens(R)
         print(io, "$(g[i]) -> $(W.d[i].coeff)")
      else
         println(io, "$(g[i]) -> $(W.d[i].coeff)")
      end
    end
  end
  print(io, Dedent())
#  println(IOContext(io, :compact => true, ), W.d)
end

function Base.show(io::IO, W::MPolyDecRing)
  Hecke.@show_name(io, W)
  Hecke.@show_special(io, W)
  io = pretty(io)
  if is_terse(io)
    if is_filtered(W)
      print(io, "Filtered multivariate polynomial ring")
    else
      print(io, "Graded multivariate polynomial ring")
    end
  else
    R = forget_decoration(W)
    if is_filtered(W)
      print(io, "Filtered ", Lowercase(), R)
    else
      print(io, "Graded ", Lowercase(), R)
    end
  end
end


function decorate(R::MPolyRing)
  A = abelian_group([0])
  S = MPolyDecRing(R, [1*A[1] for i = 1: ngens(R)], (x,y) -> x[1] < y[1])
  return S, map(R, gens(R))
end

@doc raw"""
    grade(R::MPolyRing, W::AbstractVector{<:IntegerUnion})

Given a vector `W` of `ngens(R)` integers, create a free abelian group of type `FinGenAbGroup`
given by one free generator, and convert the entries of `W` to elements of that group. Then
create a $\mathbb Z$-graded ring by assigning the group elements as weights to the variables
of `R`, and return the new ring, together with the vector of variables.

    grade(R::MPolyRing)

As above, where the grading is the standard $\mathbb Z$-grading on `R`.

# Examples
```jldoctest grade-ex
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W);

julia> S
Multivariate polynomial ring in 3 variables over QQ graded by
  x -> [1]
  y -> [2]
  z -> [3]

julia> T, (x, y, z) = grade(R);

julia> T
Multivariate polynomial ring in 3 variables over QQ graded by
  x -> [1]
  y -> [1]
  z -> [1]
```

Grading an already graded polynomial ring is not supported.
```jldoctest grade-ex
julia> grade(S)
ERROR: ArgumentError: cannot graded polynomial ring which is already decorated
```
"""
function grade(R::MPolyRing, W::AbstractVector{<:IntegerUnion})
  @assert length(W) == ngens(R)
  A = abelian_group([0])
  set_attribute!(A, :show_elem => show_special_elem_grad)
  S = MPolyDecRing(R, [i*A[1] for i = W])
  return S, map(S, gens(R))
end

function grade(R::MPolyRing)
  A = abelian_group([0])
  S = MPolyDecRing(R, [1*A[1] for i = 1: ngens(R)])
  return S, map(S, gens(R))
end

@doc raw"""
    grade(R::MPolyRing, W::AbstractVector{<:AbstractVector{<:IntegerUnion}})

Given a vector `W` of `ngens(R)` integer vectors of the same size `m`, say, create a free
abelian group of type `FinGenAbGroup` given by `m` free generators, and convert the vectors in
`W` to elements of that group. Then create a $\mathbb Z^m$-graded ring by assigning the group
elements as weights to the variables of `R`, and return the new ring, together with the vector
of variables.

    grade(R::MPolyRing, W::Union{ZZMatrix, AbstractMatrix{<:IntegerUnion}})

As above, converting the columns of `W`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> W = [1 1 0 0 0; 0 0 1 1 1]
2×5 Matrix{Int64}:
 1  1  0  0  0
 0  0  1  1  1

julia> S, x = grade(R, W);

julia> S
Multivariate polynomial ring in 5 variables over QQ graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  x[3] -> [0 1]
  x[4] -> [0 1]
  x[5] -> [0 1]

```
"""
function grade(R::MPolyRing, W::AbstractVector{<:AbstractVector{<:IntegerUnion}})
  @assert length(W) == ngens(R)
  n = length(W[1])
  @assert all(x->length(x) == n, W)

  A = abelian_group(zeros(Int, n))
  set_attribute!(A, :show_elem => show_special_elem_grad)
  S = MPolyDecRing(R, [A(w) for w = W])
  return S, map(S, gens(R))
end

function grade(R::MPolyRing, W::Union{ZZMatrix, AbstractMatrix{<:IntegerUnion}})
  @assert size(W, 2) == ngens(R)
  A = abelian_group(zeros(Int, size(W, 1)))
  set_attribute!(A, :show_elem => show_special_elem_grad)
  ###S = MPolyDecRing(R, [A(view(W, :, i)) for i = 1:size(W, 2)])
  S = MPolyDecRing(R, [A(W[:, i]) for i = 1:size(W, 2)])
  return S, map(S, gens(R))
end


@doc raw"""
    weights(R::MPolyDecRing)

Given a graded multivariate polynomial ring `R`, return the weights (degrees) of the variables of `R`.

     weights(::Type{Vector{Int}}, R::MPolyDecRing)

Given a $\mathbb Z^m$-graded multivariate polynomial ring `R`, return the weights (degrees) of the variables of `R`, converted to vectors of integer numbers.

    weights(::Type{Int}, R::MPolyDecRing)

Given a $\mathbb Z$-graded multivariate polynomial ring `R`, return the weights (degrees) of `R`, converted to integer numbers.

# Examples
```jldoctest
julia> G = abelian_group([0, 0, 2, 2])
Finitely generated abelian group
  with 4 generators and 4 relations and relation matrix
  [0   0   0   0]
  [0   0   0   0]
  [0   0   2   0]
  [0   0   0   2]

julia> W = [G[1]+G[3]+G[4], G[2]+G[4], G[1]+G[3], G[2], G[1]+G[2]];

julia> R, x = graded_polynomial_ring(QQ, :x => 1:5; weights = W);

julia> weights(R)
5-element Vector{FinGenAbGroupElem}:
 [1, 0, 1, 1]
 [0, 1, 0, 1]
 [1, 0, 1, 0]
 [0, 1, 0, 0]
 [1, 1, 0, 0]

julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]];

julia> R, x = graded_polynomial_ring(QQ, :x => 1:4, W);

julia> weights(R)
4-element Vector{FinGenAbGroupElem}:
 [1 0]
 [0 1]
 [1 0]
 [4 1]

julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> weights(R)
3-element Vector{FinGenAbGroupElem}:
 [1]
 [1]
 [1]

```
"""
function weights(R::MPolyDecRing)
 @assert is_graded(R)
 return [degree(x) for x in gens(R)]
end

function weights(::Type{Vector{Int}}, R::MPolyDecRing)
 @assert is_zm_graded(R)
 G = grading_group(R)
 ws = Vector{Int}[]
 for i = 1:ngens(R)
   wi = [Int(degree(gens(R)[i])[1]) for i=1:ngens(G)]
   push!(ws, wi)
 end
 return ws
end

function weights(::Type{Int}, R::MPolyDecRing)
 @assert is_z_graded(R)
   return [Int(degree(x)[1]) for x in gens(R)]
end

@doc raw"""
    is_standard_graded(R::MPolyDecRing)

Return `true` if `R` is standard $\mathbb Z$-graded, `false` otherwise.

# Examples
```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights = [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> is_standard_graded(S)
false
```
"""
function is_standard_graded(R::MPolyDecRing)
  is_z_graded(R) || return false
  A = grading_group(R)
  W = R.d
  for i = 1:length(W)
     if W[i] != A[1]
        return false
     end
  end
  return true
end

is_standard_graded(::MPolyRing) = false

@doc raw"""
    is_z_graded(R::MPolyDecRing)

Return `true` if `R` is $\mathbb Z$-graded, `false` otherwise.

!!! note
    Writing `G = grading_group(R)`, we say that `R` is $\mathbb Z$-graded if
    `G` is free abelian of rank `1`, and `ngens(G) == 1`.

# Examples
```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights = [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> is_z_graded(S)
true
```
"""
function is_z_graded(R::MPolyDecRing)
  is_graded(R) || return false
  A = grading_group(R)
  return ngens(A) == 1 && torsion_free_rank(A) == 1 && is_free(A)
end

is_z_graded(::MPolyRing) = false

@doc raw"""
    is_zm_graded(R::MPolyDecRing)

Return `true` if `R` is $\mathbb Z^m$-graded for some $m$, `false` otherwise.

!!! note
    Writing `G = grading_group(R)`, we say that `R` is $\mathbb Z^m$-graded
    `G` is free abelian of rank `m`, and `ngens(G) == m`.

# Examples
```jldoctest
julia> G = abelian_group([0, 0, 2, 2])
Finitely generated abelian group
  with 4 generators and 4 relations and relation matrix
  [0   0   0   0]
  [0   0   0   0]
  [0   0   2   0]
  [0   0   0   2]

julia> W = [G[1]+G[3]+G[4], G[2]+G[4], G[1]+G[3], G[2], G[1]+G[2]];

julia> S, x = graded_polynomial_ring(QQ, :x => 1:5; weights = W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> is_zm_graded(S)
false

```
```jldoctest
julia> G = abelian_group(ZZMatrix([1 -1]))
Finitely generated abelian group
  with 2 generators and 1 relation and relation matrix
  [1   -1]

julia> g = gen(G, 1)
Abelian group element [0, 1]

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z], W);

julia> R
Multivariate polynomial ring in 4 variables over QQ graded by
  w -> [0 1]
  x -> [0 1]
  y -> [0 1]
  z -> [0 1]

julia> is_free(G)
true

julia> is_zm_graded(R)
false

```
"""
function is_zm_graded(R::MPolyDecRing)
  is_graded(R) || return false
  A = grading_group(R)
  return is_free(A) && ngens(A) == torsion_free_rank(A)
end

is_zm_graded(::MPolyRing) = false

@doc raw"""
    is_positively_graded(R::MPolyDecRing)

Return `true` if `R` is positively graded, `false` otherwise.

!!! note
    We say that `R` is *positively graded* by a finitely generated abelian group $G$ if the coefficient ring of `R` is a field,
    $G$ is free, and each graded part $R_g$, $g\in G$, has finite dimension.

# Examples
```jldoctest
julia> S, (t, x, y) = graded_polynomial_ring(QQ, [:t, :x, :y]; weights = [-1, 1, 1])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[t, x, y])

julia> grading_group(S)
Z

julia> is_positively_graded(S)
false

```

```jldoctest
julia> G = abelian_group([0, 2])
Finitely generated abelian group
  with 2 generators and 2 relations and relation matrix
  [0   0]
  [0   2]

julia> is_free(G)
false

julia> W = [gen(G, 1)+gen(G, 2), gen(G, 1)]
2-element Vector{FinGenAbGroupElem}:
 [1, 1]
 [1, 0]

julia> S, (x, y) = graded_polynomial_ring(QQ, [:x, :y]; weights = W)
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> is_positively_graded(S)
false

```
"""
@attr Bool function is_positively_graded(R::MPolyDecRing)
  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  is_graded(R) || return false
  G = grading_group(R)
  is_free(G) || return false
  if ngens(G) == torsion_free_rank(G)
    W = reduce(vcat, [x.coeff for x = R.d])
    if is_positive_grading_matrix(transpose(W))
       return true
    end
  end
  try
    homogeneous_component(R, zero(G))
  catch e
    if e isa ArgumentError && e.msg == "Polyhedron not bounded"
      return false
    else
      rethrow(e)
    end
  end
  return true
end

is_positively_graded(::MPolyRing) = false

@doc raw"""
    graded_polynomial_ring(C::Ring, args...; weights = nothing, kwargs...)

Create a multivariate [`polynomial_ring`](@ref polynomial_ring(R, [:x])) with
coefficient ring `C` and variables as described by `args...` (using the exact
same syntax as for `polynomial_ring`), and [`grade`](@ref) this ring
according to the data provided by the keyword argument `weights`.
Return the graded ring as an object of type `MPolyDecRing`, together with the variables.

!!! note
    If no `weights` are entered, the returned ring is standard $\mathbb Z$-graded, i.e. all variables are graded with weight `1`.
 
!!! note
    `kwargs` allows one to set the same keywords as for `polynomial_ring`. 

# Examples
```jldoctest
julia> G = abelian_group([0, 0, 2, 2])
Finitely generated abelian group
  with 4 generators and 4 relations and relation matrix
  [0   0   0   0]
  [0   0   0   0]
  [0   0   2   0]
  [0   0   0   2]

julia> W1 = [G[1]+G[3]+G[4], G[2]+G[4], G[1]+G[3], G[2], G[1]+G[2]];
 
julia> R1, x, y = graded_polynomial_ring(QQ, :x => 1:2, :y => 1:3, W1);

julia> R1
Multivariate polynomial ring in 5 variables over QQ graded by
  x[1] -> [1 0 1 1]
  x[2] -> [0 1 0 1]
  y[1] -> [1 0 1 0]
  y[2] -> [0 1 0 0]
  y[3] -> [1 1 0 0]

julia> W2 = [[1, 0], [0, 1], [1, 0], [4, 1]];

julia> R2, x = graded_polynomial_ring(QQ, 4, :x; weights = W2);

julia> R2
Multivariate polynomial ring in 4 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [0 1]
  x3 -> [1 0]
  x4 -> [4 1]

julia> R3, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights = [1, 2, 3]);

julia> R3
Multivariate polynomial ring in 3 variables over QQ graded by
  x -> [1]
  y -> [2]
  z -> [3]

julia> R4, x = graded_polynomial_ring(QQ, :x => 1:3);

julia> R4
Multivariate polynomial ring in 3 variables over QQ graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

julia> R5, x, y = graded_polynomial_ring(QQ, :x => 1:3, :y => (1:2, 1:2); weights = 1:7);

julia> R5
Multivariate polynomial ring in 7 variables over QQ graded by
  x[1] -> [1]
  x[2] -> [2]
  x[3] -> [3]
  y[1, 1] -> [4]
  y[2, 1] -> [5]
  y[1, 2] -> [6]
  y[2, 2] -> [7]

```
"""
function graded_polynomial_ring(C::Ring, args...; weights = nothing, kwargs...)
  if weights === nothing
    # no weights kwarg given: for backwards compatibility also check if
    # the last regular argument might be a weight vector and if so, use it.
    if args[end] isa Union{Vector{<:IntegerUnion}, Vector{<:Vector{<:IntegerUnion}}, Matrix{<:IntegerUnion}, ZZMatrix, Vector{FinGenAbGroupElem}}
      weights = args[end]
      args = args[1:end-1]
    end
  end

  # pass all arguments (except possibly the last one if it contains weights)
  # on to polynomial_ring
  R, v... = polynomial_ring(C, args...; kwargs...)

  # if no weights were given as last argument or via the `weights` kwarg,
  # then we now use a default value where we assign weight 1 to each variable
  if weights === nothing
    weights = ones(Int, nvars(R))
  end

  # TODO: MPolyDecRing resp. `grade` should support `cached` kwarg
  S, = grade(R, weights)

  # return a result of the same "shape" as that returned by polynomial_ring
  return S, _map_recursive(S, v)...
end

# helper for graded_polynomial_ring
_map_recursive(S::NCRing, x::NCRingElem) = S(x)
_map_recursive(S::NCRing, a::AbstractArray) = map(x -> _map_recursive(S, x), a)
_map_recursive(S::NCRing, a::Tuple) = map(x -> _map_recursive(S, x), a)

filtrate(R::MPolyRing) = decorate(R)

function show_special_elem_grad(io::IO, a::FinGenAbGroupElem)
  if get(io, :compact, false)
    print(io, a.coeff)
  else
    print(io, "$(a.coeff)")
  end
end

function filtrate(R::MPolyRing, v::Vector{Int})
  A = abelian_group([0])
  set_attribute!(A, :show_elem => show_special_elem_grad)
  S = MPolyDecRing(R, [i*A[1] for i = v], (x,y) -> x[1] < y[1])
  return S, map(S, gens(R))
end

function filtrate(R::MPolyRing, v::Vector{FinGenAbGroupElem}, lt)
  S = MPolyDecRing(R, v, lt)
  return S, map(S, gens(R))
end

@doc raw"""
    grade(R::MPolyRing, W::Vector{FinGenAbGroupElem})

Given a vector `W` of `ngens(R)` elements of a finitely generated abelian group `G`, say, create 
a `G`-graded ring by assigning the entries of `W` as weights to the variables of `R`. Return
the new ring as an object of type `MPolyDecRing`, together with the vector of variables.

# Examples
```jldoctest
julia> R, (t, x, y) = polynomial_ring(QQ, [:t, :x, :y])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[t, x, y])

julia> typeof(R)
QQMPolyRing

julia>  typeof(x)
QQMPolyRingElem

julia> G = abelian_group([0])
Z

julia> g = gen(G, 1)
Abelian group element [1]

julia> S, (t, x, y) = grade(R, [-g, g, g])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[t, x, y])

julia> S
Multivariate polynomial ring in 3 variables over QQ graded by
  t -> [-1]
  x -> [1]
  y -> [1]

julia> typeof(S)
MPolyDecRing{QQFieldElem, QQMPolyRing}

julia> S isa MPolyRing
true

julia> typeof(x)
MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}

```

```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0])
Z^2

julia> g = gens(G)
2-element Vector{FinGenAbGroupElem}:
 [1, 0]
 [0, 1]

julia> W = [g[1], g[1], g[2], g[2], g[2]];

julia> S, _ = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> S
Multivariate polynomial ring in 5 variables over QQ graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  x[3] -> [0 1]
  x[4] -> [0 1]
  x[5] -> [0 1]

julia> typeof(x[1])
QQMPolyRingElem

julia> x = map(S, x)
5-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x[1]
 x[2]
 x[3]
 x[4]
 x[5]

julia> typeof(x[1])
MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}

```
```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
Finitely generated abelian group
  with 4 generators and 4 relations and relation matrix
  [0   0   0   0]
  [0   0   0   0]
  [0   0   2   0]
  [0   0   0   2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]]
5-element Vector{FinGenAbGroupElem}:
 [1, 0, 1, 1]
 [0, 1, 0, 1]
 [1, 0, 1, 0]
 [0, 1, 0, 0]
 [1, 1, 0, 0]

julia> S, x = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> S
Multivariate polynomial ring in 5 variables over QQ graded by
  x[1] -> [1 0 1 1]
  x[2] -> [0 1 0 1]
  x[3] -> [1 0 1 0]
  x[4] -> [0 1 0 0]
  x[5] -> [1 1 0 0]

```
"""
function grade(R::MPolyRing, v::AbstractVector{FinGenAbGroupElem})
  S = MPolyDecRing(R, v)
  return S, map(S, gens(R))
end

mutable struct MPolyDecRingElem{T, S} <: MPolyRingElem{T}
  f::S
  parent
  function MPolyDecRingElem(f::S, p) where {S}
    r = new{elem_type(base_ring(f)), S}(f, p)
#    if is_graded(p) && length(r) > 1
#      if !is_homogeneous(r)
#        error("element not homogeneous")
#      end
#both wrong and undesired.
#    end
    return r
  end
end

function Base.deepcopy_internal(f::MPolyDecRingElem{T, S}, dict::IdDict) where {T, S}
  return MPolyDecRingElem(Base.deepcopy_internal(forget_decoration(f), dict), f.parent)
end

function show(io::IO, w::MPolyDecRingElem)
  show(io, forget_decoration(w))
end

parent(a::MPolyDecRingElem{T, S}) where {T, S} = a.parent::MPolyDecRing{T, parent_type(S)}

symbols(R::MPolyDecRing) = symbols(forget_decoration(R))
number_of_variables(R::MPolyDecRing) = number_of_variables(forget_decoration(R))

elem_type(::Type{MPolyDecRing{T, S}}) where {T, S} = MPolyDecRingElem{T, elem_type(S)}
parent_type(::Type{MPolyDecRingElem{T, S}}) where {T, S} = MPolyDecRing{T, parent_type(S)}

(W::MPolyDecRing)() = MPolyDecRingElem(forget_decoration(W)(), W)
(W::MPolyDecRing)(i::Int) = MPolyDecRingElem(forget_decoration(W)(i), W)
(W::MPolyDecRing)(f::Singular.spoly) = MPolyDecRingElem(forget_decoration(W)(f), W)

### Coercion of elements of the underlying polynomial ring
# into the graded ring.
function (W::MPolyDecRing{S, T})(f::U) where {S, T, U<:MPolyRingElem}
  if parent_type(U) === T
    @assert forget_decoration(W) === parent(f)
    return MPolyDecRingElem(f, W)
  else
    return W(forget_decoration(W)(f))
  end
end

function (W::MPolyDecRing)(f)
  return W(forget_decoration(W)(f))
end


function (W::MPolyDecRing{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T
  return W(forget_decoration(W)(c, e))
end

(W::MPolyDecRing)(g::MPolyDecRingElem; check::Bool=true) = MPolyDecRingElem(forget_decoration(W)(forget_decoration(g)), W)
one(W::MPolyDecRing) = MPolyDecRingElem(one(forget_decoration(W)), W)
zero(W::MPolyDecRing) = MPolyDecRingElem(zero(forget_decoration(W)), W)

################################################################################
#
#  Binary operations
#
################################################################################

for T in [:(+), :(-), :(*)]
  @eval ($T)(a::MPolyDecRingElem,
             b::MPolyDecRingElem) = MPolyDecRingElem($T(forget_decoration(a), forget_decoration(b)), parent(a))
end

divexact(a::MPolyDecRingElem, b::MPolyDecRingElem; check::Bool=true) = MPolyDecRingElem(divexact(forget_decoration(a), forget_decoration(b); check=check), parent(a))


################################################################################
#
#  Unitary operations
#
################################################################################

-(a::MPolyDecRingElem)   = MPolyDecRingElem(-forget_decoration(a), parent(a))

################################################################################
#
#  Binary ad hoc operations
#
################################################################################

divexact(a::MPolyDecRingElem, b::RingElem; check::Bool=true) = MPolyDecRingElem(divexact(forget_decoration(a), b; check=check), parent(a))

divexact(a::MPolyDecRingElem, b::Integer; check::Bool=true) = MPolyDecRingElem(divexact(forget_decoration(a), b; check=check), parent(a))

divexact(a::MPolyDecRingElem, b::Rational; check::Bool=true) = MPolyDecRingElem(divexact(forget_decoration(a), b; check=check), parent(a))

for T in [:(-), :(+)]
  @eval ($T)(a::MPolyDecRingElem,
             b::RingElem) = MPolyDecRingElem($(T)(forget_decoration(a), b), parent(a))

  @eval ($T)(a::MPolyDecRingElem,
             b::Integer) = MPolyDecRingElem($(T)(forget_decoration(a), b), parent(a))

  @eval ($T)(a::MPolyDecRingElem,
             b::Rational) = MPolyDecRingElem($(T)(forget_decoration(a), b), parent(a))

  @eval ($T)(a::RingElem,
             b::MPolyDecRingElem) = MPolyDecRingElem($(T)(a, forget_decoration(b)), b.parent)

  @eval ($T)(a::Integer,
             b::MPolyDecRingElem) = MPolyDecRingElem($(T)(a, forget_decoration(b)), b.parent)

  @eval ($T)(a::Rational,
             b::MPolyDecRingElem) = MPolyDecRingElem($(T)(a, forget_decoration(b)), b.parent)
end

function *(a::MPolyDecRingElem{T, S}, b::T) where {T <: RingElem, S}
  return MPolyDecRingElem(forget_decoration(a) * b, parent(a))
end

function *(a::T, b::MPolyDecRingElem{T, S}) where {T <: RingElem, S}
  return b * a
end

################################################################################
#
#  Factoring, division, ...
#
################################################################################

function factor(x::MPolyDecRingElem)
  R = parent(x)
  F = factor(forget_decoration(x))
  D = Dict{elem_type(R), Int}(R(i) => e for (i, e) in F)
  return Fac(R(unit(F)), D)
end

function factor_squarefree(x::MPolyDecRingElem)
  R = parent(x)
  F = factor_squarefree(forget_decoration(x))
  D = Dict{elem_type(R), Int}(R(i) => e for (i, e) in F)
  return Fac(R(unit(F)), D)
end

function gcd(x::MPolyDecRingElem, y::MPolyDecRingElem)
  R = parent(x)
  return R(gcd(forget_decoration(x), forget_decoration(y)))
end

function div(x::MPolyDecRingElem, y::MPolyDecRingElem)
  R = parent(x)
  return R(div(forget_decoration(x), forget_decoration(y)))
end

function divrem(x::MPolyDecRingElem, y::MPolyDecRingElem)
  R = parent(x)
  q, r = divrem(forget_decoration(x), forget_decoration(y))
  return R(q), R(r)
end

################################################################################
#
#  Equality
#
################################################################################

==(a::MPolyDecRingElem, b::MPolyDecRingElem) = forget_decoration(a) == forget_decoration(b)

^(a::MPolyDecRingElem, i::Int) = MPolyDecRingElem(forget_decoration(a)^i, parent(a))

length(a::MPolyDecRingElem) = length(forget_decoration(a))

@doc raw"""
    total_degree(f::MPolyDecRingElem)

Return the total degree of `f`.

Given a set of variables ``x = \{x_1, \ldots, x_n\}``, the *total degree* of a monomial ``x^\alpha=x_1^{\alpha_1}\cdots x_n^{\alpha_n}\in\text{Mon}_n(x)`` is the sum of the ``\alpha_i``. The *total degree* of a polynomial `f`  is the maximum of the total degrees of its monomials.

!!! note
    The notion of total degree does not dependent on weights given to the variables.
"""
total_degree(a::MPolyDecRingElem) = total_degree(forget_decoration(a))

AbstractAlgebra.monomial(a::MPolyDecRingElem, i::Int) = parent(a)(AbstractAlgebra.monomial(forget_decoration(a), i))
AbstractAlgebra.coeff(a::MPolyDecRingElem, i::Int) = AbstractAlgebra.coeff(forget_decoration(a), i)
AbstractAlgebra.term(a::MPolyDecRingElem, i::Int) = parent(a)(AbstractAlgebra.term(forget_decoration(a), i))
AbstractAlgebra.exponent_vector(a::MPolyDecRingElem, i::Int) = AbstractAlgebra.exponent_vector(forget_decoration(a), i)
AbstractAlgebra.exponent_vector(a::MPolyDecRingElem, i::Int, ::Type{T}) where T = AbstractAlgebra.exponent_vector(forget_decoration(a), i, T)
AbstractAlgebra.exponent(a::MPolyDecRingElem, i::Int, j::Int) = AbstractAlgebra.exponent(forget_decoration(a), i, j)
AbstractAlgebra.exponent(a::MPolyDecRingElem, i::Int, j::Int, ::Type{T}) where T = AbstractAlgebra.exponent(forget_decoration(a), i, j, T)

function has_weighted_ordering(R::MPolyDecRing)
  grading_to_ordering = false
  w_ord = degrevlex(R) # dummy, not used
  # This is not meant to be exhaustive, there a probably more gradings which one
  # can meaningfully translate into a monomial ordering
  # However, we want to stick to global orderings.
  if is_z_graded(R)
    w = Int[ R.d[i].coeff[1] for i = 1:ngens(R) ]
    if all(isone, w)
      w_ord = degrevlex(R)
      grading_to_ordering = true
    elseif all(>(0), w)
      w_ord = wdegrevlex(R, w)
      grading_to_ordering = true
    end
  end
  return grading_to_ordering, w_ord
end

@attr MonomialOrdering{T} function default_ordering(R::T) where {T<:MPolyDecRing}
  fl, w_ord = has_weighted_ordering(R)
  fl && return w_ord
  return degrevlex(R)
end

function singular_poly_ring(R::MPolyDecRing; keep_ordering::Bool = false)
  if !keep_ordering
    return singular_poly_ring(forget_decoration(R), default_ordering(R))
  end
  return singular_poly_ring(forget_decoration(R); keep_ordering)
end

MPolyCoeffs(f::MPolyDecRingElem) = MPolyCoeffs(forget_decoration(f))
MPolyExponentVectors(f::MPolyDecRingElem) = MPolyExponentVectors(forget_decoration(f))

function push_term!(M::MPolyBuildCtx{<:MPolyDecRingElem{T, S}}, c::T, expv::Vector{Int}) where {T <: RingElement, S}
  if iszero(c)
    return M
  end
  len = length(M.poly.f) + 1
  set_exponent_vector!(M.poly.f, len, expv)
  setcoeff!(M.poly.f, len, c)
  return M
end

function set_exponent_vector!(f::MPolyDecRingElem, i::Int, exps::Vector{Int})
  f.f = set_exponent_vector!(forget_decoration(f), i, exps)
  return f
end

function finish(M::MPolyBuildCtx{<:MPolyDecRingElem})
  f = sort_terms!(M.poly.f)
  f = combine_like_terms!(M.poly.f)
  return parent(M.poly)(f)
end

function jacobian_matrix(f::MPolyDecRingElem)
  R = parent(f)
  n = nvars(R)
  return matrix(R, n, 1, [derivative(f, i) for i=1:n])
end

function jacobian_ideal(f::MPolyDecRingElem)
  R = parent(f)
  n = nvars(R)
  return ideal(R, [derivative(f, i) for i=1:n])
end

function jacobian_matrix(g::Vector{<:MPolyDecRingElem})
  R = parent(g[1])
  n = nvars(R)
  @assert all(x->parent(x) === R, g)
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

@doc raw"""
    degree(f::MPolyDecRingElem)

Given a homogeneous element `f` of a graded multivariate ring, return the degree of `f`.

    degree(::Type{Vector{Int}}, f::MPolyDecRingElem)

Given a homogeneous element `f` of a $\mathbb Z^m$-graded multivariate ring, return the degree of `f`, converted to a vector of integer numbers.

    degree(::Type{Int}, f::MPolyDecRingElem)

Given a homogeneous element `f` of a $\mathbb Z$-graded multivariate ring, return the degree of `f`, converted to an integer number.

# Examples
```jldoctest
julia> G = abelian_group([0, 0, 2, 2])
Finitely generated abelian group
  with 4 generators and 4 relations and relation matrix
  [0   0   0   0]
  [0   0   0   0]
  [0   0   2   0]
  [0   0   0   2]

julia> W = [G[1]+G[3]+G[4], G[2]+G[4], G[1]+G[3], G[2], G[1]+G[2]];

julia> S, x = graded_polynomial_ring(QQ, :x => 1:5; weights = W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[2]^2+2*x[4]^2
x[2]^2 + 2*x[4]^2

julia> degree(f)
Abelian group element [0, 2, 0, 0]

julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 0]
 [4, 1]

julia> R, x = graded_polynomial_ring(QQ, :x => 1:4, W)
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^4*x[2]+x[4]
x[1]^4*x[2] + x[4]

julia> degree(f)
[4 1]

julia> degree(Vector{Int}, f)
2-element Vector{Int64}:
 4
 1

julia>  R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^6+y^3+z^2
x^6 + y^3 + z^2

julia> degree(f)
[6]

julia> typeof(degree(f))
FinGenAbGroupElem

julia> degree(Int, f)
6

julia> typeof(degree(Int, f))
Int64
```
"""
function degree(a::MPolyDecRingElem; check::Bool=true)
  !check && !is_filtered(parent(a)) && return _degree_fast(a)
  # TODO: Also provide a fast track for the filtered case.
  @req !iszero(a) "Element must be non-zero"
  W = parent(a)
  w = W.D[0]
  first = true
  d = W.d
  for c = MPolyExponentVectors(forget_decoration(a))
    u = W.D[0]
    for i=1:length(c)
      u += c[i]*d[i]
    end
    if first
      first = false
      w = u
    elseif is_filtered(W)
      w = W.lt(w, u) ? u : w
    else
      w == u || error("element not homogeneous")
    end
  end
  return w
end

function _degree_fast(a::MPolyDecRingElem)
  f = forget_grading(a)
  w = parent(a).d
  z = zero(grading_group(parent(a)))
  is_zero(f) && return z
  for (c, e) in zip(coefficients(f), exponents(f))
    !iszero(c) && return sum(b*w[i] for (i, b) in enumerate(e); init=z)
  end
end

function degree(::Type{Int}, a::MPolyDecRingElem; check::Bool=true)
  @assert is_z_graded(parent(a))
  return Int(degree(a; check)[1])
end

function degree(::Type{Vector{Int}}, a::MPolyDecRingElem); check::Bool=true
  @assert is_zm_graded(parent(a))
  d = degree(a; check)
  return Int[d[i] for i=1:ngens(parent(d))]
end

@doc raw"""
    is_homogeneous(f::MPolyDecRingElem)

Given an element `f` of a graded multivariate ring, return `true` if `f` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^2+y*z
x^2 + y*z

julia> is_homogeneous(f)
false

julia> W = [1 2 1 0; 3 4 0 1]
2×4 Matrix{Int64}:
 1  2  1  0
 3  4  0  1

julia> S, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z], W)
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w, x, y, z])

julia> F = w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3
w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3

julia> is_homogeneous(F)
true
```
"""
function is_homogeneous(F::MPolyDecRingElem)
  D = parent(F).D
  d = parent(F).d
  S = nothing
  u = zero(D)
  for c = MPolyExponentVectors(forget_decoration(F))
    u = zero!(u)
    for i=1:length(c)
      u = addmul_delayed_reduction!(u, d[i], c[i])
    end
    u = reduce!(u)
    if S === nothing
      S = deepcopy(u)
    elseif S != u
      return false
    end
  end
  return true
end

@doc raw"""
    homogeneous_components(f::MPolyDecRingElem{T, S}) where {T, S}

Given an element `f` of a graded multivariate ring, return the homogeneous components of `f`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^2+y+z
x^2 + y + z

julia> homogeneous_components(f)
Dict{FinGenAbGroupElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  [2] => x^2 + y
  [3] => z

julia> G = abelian_group([0, 0, 2, 2])
Finitely generated abelian group
  with 4 generators and 4 relations and relation matrix
  [0   0   0   0]
  [0   0   0   0]
  [0   0   2   0]
  [0   0   0   2]

julia> W = [G[1]+G[3]+G[4], G[2]+G[4], G[1]+G[3], G[2], G[1]+G[2]];

julia> S, x = graded_polynomial_ring(QQ, :x => 1:5; weights = W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[1]^2+x[3]^2+x[5]^2
x[1]^2 + x[3]^2 + x[5]^2

julia> homogeneous_components(f)
Dict{FinGenAbGroupElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  [2, 2, 0, 0] => x[5]^2
  [2, 0, 0, 0] => x[1]^2 + x[3]^2
```
"""
function homogeneous_components(a::MPolyDecRingElem{T, S}) where {T, S}
  D = parent(a).D
  d = parent(a).d
  h = Dict{elem_type(D), typeof(a)}()
  W = parent(a)
  R = forget_decoration(W)
  # First assemble the homogeneous components into the build contexts.
  # Afterwards compute the polynomials.
  hh = Dict{elem_type(D), MPolyBuildCtx{S, DataType}}()
  dmat = reduce(vcat, [d[i].coeff for i in 1:length(d)])
  tmat = zero_matrix(ZZ, 1, nvars(R))
  res_mat = zero_matrix(ZZ, 1, ncols(dmat))
  for (c, e) = Base.Iterators.zip(AbstractAlgebra.coefficients(forget_decoration(a)), AbstractAlgebra.exponent_vectors(forget_decoration(a)))
    # this is non-allocating
    for i in 1:length(e)
      tmat[1, i] = e[i]
    end
    mul!(res_mat, tmat, dmat)
    u = FinGenAbGroupElem(D, res_mat)
    if haskey(hh, u)
      ctx = hh[u]
      push_term!(ctx, c, e)
    else
      # We put u in the dictionary
      # Make a fresh res_mat, which can be used the for the next u
      res_mat = deepcopy(res_mat)
      ctx = MPolyBuildCtx(R)
      push_term!(ctx, c, e)
      hh[u] = ctx
    end
  end
  hhh = Dict{elem_type(D), typeof(a)}()
  for (u, C) in hh
    hhh[u] = W(finish(C))
  end

  return hhh
end

@doc raw"""
    homogeneous_component(f::MPolyDecRingElem, g::FinGenAbGroupElem)

Given an element `f` of a graded multivariate ring, and given an element
`g` of the grading group of that ring, return the homogeneous component of `f` of degree `g`.

    homogeneous_component(f::MPolyDecRingElem, g::Vector{<:IntegerUnion})

Given an element `f` of a $\mathbb  Z^m$-graded multivariate ring `R`, say, and given
a vector `g` of $m$ integers, convert `g` into an element of the grading group of `R`,
and return the homogeneous component of `f` whose degree is that element.

    homogeneous_component(f::MPolyDecRingElem, g::IntegerUnion)

Given an element `f` of a $\mathbb  Z$-graded multivariate ring `R`, say, and given
an integer `g`, convert `g` into an element of the grading group of `R`, and return the
homogeneous component of `f` whose degree is that element.

# Examples
```jldoctest
julia> G = abelian_group([0, 0, 2, 2])
Finitely generated abelian group
  with 4 generators and 4 relations and relation matrix
  [0   0   0   0]
  [0   0   0   0]
  [0   0   2   0]
  [0   0   0   2]

julia> W = [G[1]+G[3]+G[4], G[2]+G[4], G[1]+G[3], G[2], G[1]+G[2]];

julia> S, x = graded_polynomial_ring(QQ, :x => 1:5; weights = W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[1]^2+x[3]^2+x[5]^2
x[1]^2 + x[3]^2 + x[5]^2

julia> homogeneous_component(f, 2*G[1])
x[1]^2 + x[3]^2

```

```jldoctest
julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 0]
 [4, 1]

julia> R, x = graded_polynomial_ring(QQ, :x => 1:4; weights = W)
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^2*x[2]+x[4]
x[1]^2*x[2] + x[4]

julia> homogeneous_component(f, [2, 1])
x[1]^2*x[2]

julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights = [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^2+y+z
x^2 + y + z

julia> homogeneous_component(f, 1)
0

julia> homogeneous_component(f, 2)
x^2 + y

julia> homogeneous_component(f, 3)
z
```
"""
function homogeneous_component(a::MPolyDecRingElem, g::FinGenAbGroupElem)
  R = forget_decoration(parent(a))
  r = R(0)
  d = parent(a).d
  for (c, m) = Base.Iterators.zip(MPolyCoeffs(forget_decoration(a)), Generic.MPolyMonomials(forget_decoration(a)))
    e = exponent_vector(m, 1)
    u = parent(a).D[0]
    for i=1:length(e)
      u += e[i]*d[i]
    end
    if u == g
      r += c*m
    end
  end
  return parent(a)(r)
end

function homogeneous_component(a::MPolyDecRingElem, g::IntegerUnion)
  @assert is_z_graded(parent(a))
  return homogeneous_component(a, grading_group(parent(a))([g]))
end

function homogeneous_component(a::MPolyDecRingElem, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(parent(a))
  return homogeneous_component(a, grading_group(parent(a))(g))
end

base_ring(W::MPolyDecRing) = base_ring(forget_decoration(W))
base_ring_type(::Type{MPolyDecRing{T, S}}) where {T, S} = base_ring_type(S)
number_of_generators(W::MPolyDecRing) = number_of_generators(forget_decoration(W))
gens(W::MPolyDecRing) = map(W, gens(forget_decoration(W)))
gen(W::MPolyDecRing, i::Int) = W(gen(forget_decoration(W), i))
is_gen(a::MPolyDecRingElem) = is_gen(forget_grading(a))

function show_homo_comp(io::IO, M)
  (W, d) = get_attribute(M, :data)
  n = AbstractAlgebra.get_name(W)
  io = pretty(io)
  if n !== nothing
    print(io, LowercaseOff(), "$(n)_$(d.coeff) of dim $(dim(M))")
  else
    print(io, "homogeneous component of ", Lowercase(), W, " of degree ")
    print(IOContext(io, :compact => true), d)
  end
end

@doc raw"""
    monomial_basis(R::MPolyDecRing, g::FinGenAbGroupElem)

Given a polynomial ring `R` over a field which is graded by a free
group of type `FinGenAbGroup`, and given an element `g` of that group,
return the monomials of degree `g` in `R`.

    monomial_basis(R::MPolyDecRing, W::Vector{<:IntegerUnion})

Given a $\mathbb  Z^m$-graded polynomial ring `R` over a field and
a vector `W` of $m$ integers, convert `W` into an element `g` of the grading
group of `R` and proceed as above.

    monomial_basis(R::MPolyDecRing, d::IntegerUnion)

Given a $\mathbb  Z$-graded polynomial ring `R` over a field and
an integer `d`, convert `d` into an element `g` of the grading
group of `R` and proceed as above.

!!! note
    If the component of the given degree is not finite dimensional, an error message will be thrown.

# Examples
```jldoctest
julia> T, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> G = grading_group(T)
Z

julia> L = monomial_basis(T, 2)
6-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 z^2
 y*z
 y^2
 x*z
 x*y
 x^2
```
"""
function monomial_basis(W::MPolyDecRing, d::FinGenAbGroupElem)
  #TODO: lazy: ie. no enumeration of points
  #      apparently it is possible to get the number of points faster than the points
  #TODO: in the presence of torsion, this is wrong. The component
  #      would be a module over the deg-0-sub ring.
  @req coefficient_ring(W) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  D = W.D
  is_free(D) || error("Grading group must be free")
  h = hom(free_abelian_group(ngens(W)), D, W.d)
  fl, p = has_preimage_with_preimage(h, d)
  R = base_ring(W)
  B = elem_type(W)[]
  if fl
     k, im = kernel(h)
     #need the positive elements in there...
     #Ax = b, Cx >= 0
     C = identity_matrix(ZZ, ngens(W))
     A = reduce(vcat, [x.coeff for x = W.d])
     k = solve_mixed(transpose(A), transpose(d.coeff), C)
     for ee = 1:nrows(k)
       e = k[ee, :]
       a = MPolyBuildCtx(forget_decoration(W))
       push_term!(a, R(1), [Int(e[i]) for i in 1:length(e)])
       push!(B, W(finish(a)))
     end
  end
  return B
end


function monomial_basis(R::MPolyDecRing, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(R)
  return monomial_basis(R, grading_group(R)(g))
end

function monomial_basis(R::MPolyDecRing, g::IntegerUnion)
  @assert is_z_graded(R)
  return monomial_basis(R, grading_group(R)([g]))
end

@doc raw"""
    homogeneous_component(R::MPolyDecRing, g::FinGenAbGroupElem)

Given a polynomial ring `R` over a field which is graded by a free
group, and given an element `g` of that group,
return the homogeneous component of `R` of degree `g` as a standard
vector space. Additionally, return the map which sends an element
of that vector space to the corresponding monomial in `R`.

    homogeneous_component(R::MPolyDecRing, W::Vector{<:IntegerUnion})

Given a $\mathbb  Z^m$-graded polynomial ring `R` over a field, and given
a vector `W` of $m$ integers, convert `W` into an element `g` of the grading
group of `R` and proceed as above.

    homogeneous_component(R::MPolyDecRing, d::IntegerUnion)

Given a $\mathbb  Z$-graded polynomial ring `R` over a field, and given
an integer `d`, convert `d` into an element `g` of the grading group of `R`
proceed as above.

!!! note
    If the component is not finite dimensional, an error will be thrown.

# Examples
```jldoctest
julia> W = [1 1 0 0 0; 0 0 1 1 1]
2×5 Matrix{Int64}:
 1  1  0  0  0
 0  0  1  1  1

julia> S, _ = graded_polynomial_ring(QQ, :x => 1:2, :y => 1:3; weights = W);

julia> G = grading_group(S)
Z^2

julia> L = homogeneous_component(S, [1, 1]);

julia> L[1]
S_[1 1] of dim 6

julia> FG = gens(L[1]);

julia> EMB = L[2]
Map defined by a julia-function with inverse
  from S_[1 1] of dim 6
  to graded multivariate polynomial ring in 5 variables over QQ

julia> for i in 1:length(FG) println(EMB(FG[i])) end
x[2]*y[3]
x[2]*y[2]
x[2]*y[1]
x[1]*y[3]
x[1]*y[2]
x[1]*y[1]
```
"""
function homogeneous_component(W::MPolyDecRing, d::FinGenAbGroupElem)
  #TODO: lazy: ie. no enumeration of points
  #      apparently it is possible to get the number of points faster than the points
  #TODO: in the presence of torsion, this is wrong. The component
  #      would be a module over the deg-0-sub ring.
  R = base_ring(W)
  B = monomial_basis(W, d)
  M, h = vector_space(R, B, target = W)
  set_attribute!(M, :show => show_homo_comp, :data => (W, d))
  add_relshp(M, W, x -> sum(x[i] * B[i] for i=1:length(B)))
#  add_relshp(W, M, g)
  return M, h
end

function homogeneous_component(R::MPolyDecRing, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(R)
  return homogeneous_component(R, grading_group(R)(g))
end

function homogeneous_component(R::MPolyDecRing, g::IntegerUnion)
  @assert is_z_graded(R)
  return homogeneous_component(R, grading_group(R)([g]))
end


@doc raw"""
    vector_space(K::Field, polys::Vector{T}; target = nothing) where {T <: MPolyRingElem}

Return a `K`-vector space `V` and an isomorphism from `V` to the `K`-vector
space spanned by the polynomials in `polys`.

Note that all polynomials must have the same parent `R`, and `K` must be equal
to `base_ring(R)`. The optional keyword argument `target` can be used to
specify `R` when `polys` is empty.

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> polys = [x + y, x - y, 2x + 3y];

julia> V, VtoPoly = vector_space(QQ, polys)
(Vector space of dimension 2 over QQ, Map: V -> R)

julia> VtoPoly.(basis(V))
2-element Vector{QQMPolyRingElem}:
 x
 y
```
"""
function vector_space(K::Field, polys::Vector{T}; target = nothing) where {T <: MPolyRingElem}
  local R
  if length(polys) == 0
    R = target
    @assert R !== nothing
  else
    R = parent(polys[1])
    @assert target === nothing || target === R
  end
  @assert base_ring(R) == K
  expvec = Dict{Vector{Int}, Int}()
  expvec_idx = Vector{Vector{Int}}()
  M = sparse_matrix(K)

  # take polynomials and turn them into sparse row vectors
  for f in polys
    pos = Vector{Int}()
    val = Vector{elem_type(K)}()
    for (c, e) = Base.Iterators.zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
      i = get!(expvec, e) do
        push!(expvec_idx, e)
        length(expvec_idx)
      end
      push!(pos, i)
      push!(val, c)
    end
    push!(M, sparse_row(K, pos, val))
  end

  # row reduce
  d = rref!(M)

  # turn the reduced sparse rows back into polynomials
  b = Vector{elem_type(R)}(undef, d)
  for i in 1:d
    s = MPolyBuildCtx(R)
    for (k,v) in M[i]
      push_term!(s, v, expvec_idx[k])
    end
    b[i] = finish(s)
  end

  # create a standard vector space of the right dimension
  F = free_module(K, length(b); cached = false)

  # helper computing images
  function img(x)
    sum([x[i] * b[i] for i in 1:length(b) if !is_zero_entry(x.v, 1, i)]; init = zero(R))
  end

  # helper computing preimages: takes a polynomial and maps it to a vector in F
  function pre(x::T)
    @assert parent(x) == R
    pos = Vector{Int}()
    val = Vector{elem_type(K)}()
    for (c, e) = Base.Iterators.zip(AbstractAlgebra.coefficients(x), AbstractAlgebra.exponent_vectors(x))
      i = get(expvec, e) do
        error("not in image")
      end
      push!(pos, i)
      push!(val, c)
    end
    v = sparse_row(K, pos, val)
    fl, a = can_solve_with_solution(M, v)
    if !fl
      error("not in image")
    end
    return F(dense_row(a, length(b)))
  end

  h = MapFromFunc(F, R, img, pre)
  return F, h
end

###########################################
# needs re-thought
function (W::MPolyDecRing)(m::Generic.FreeModuleElem)
  h = has_relshp(parent(m), W)
  if h !== nothing
    return h(m)
  end
  error("no coercion possible")
end

#########################################
function add_relshp(R, S, h)
  #this assumes that h is essentially a canonical map from R -> S
  D = get_attribute!(Dict{Any, Any}, R, :relshp)::Dict{Any, Any}
  if haskey(D, S)
    error("try to add double")
  end
  D[S] = h
end

function has_relshp(R, S)
  r = get_attribute(R, :relshp)
  if r === nothing
    return r
  end
  if haskey(r, S)
    return r[S]
  end
  #now the hard bit: traverse the graph, not falling into cycles.
end

############################################################################
############################################################################
############################################################################

###############################################################################
### z-graded Hilbert series stuff using Singular for finding the Hilbert series
###############################################################################

mutable struct HilbertData
  coeffs::Vector{BigInt}
  weights::Vector{Int32}
  I::MPolyIdeal
  function HilbertData(I::MPolyIdeal)
    R = base_ring(I)
    @req is_z_graded(R) "The base ring must be ZZ-graded"

    W = R.d
    W = [Int(W[i][1]) for i = 1:ngens(R)]

    @req minimum(W) > 0 "The weights must be positive"
    @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
    @req all(is_homogeneous, gens(I)) "The generators of the ideal must be homogeneous"

    S = singular_groebner_generators(I, false, true)
    cf = Singular.hilbert_series_data(S, W)
    return new(cf, W, I)
  end
  function HilbertData(B::IdealGens)
    return HilbertData(Oscar.MPolyIdeal(B))
  end
end

function hilbert_series(H::HilbertData)
  Zt, t = ZZ[:t]
  den = prod([1-t^w for w in H.weights])
  h = Zt(map(ZZRingElem, H.coeffs[1:end]))
  return h, den
end

function hilbert_series_reduced(H::HilbertData)
  h, den = hilbert_series(H)
  _, h, den = gcd_with_cofactors(h, den)
  return den(0)*h, den(0)*den
end

#Decker-Lossen, p23/24
function hilbert_polynomial(H::HilbertData)

  @req all(isone, H.weights) "All weights must be 1"

  q, dn = hilbert_series_reduced(H)
  a = QQFieldElem[]
  nf = QQFieldElem(1)
  d = degree(dn)-1
  for i=1:d+1
    push!(a, q(1)//nf)
    nf *= i
    q = derivative(q)
  end
  Qt, t = QQ[:t]
  d==-1 && return zero(Qt)
  bin = one(Qt)
  b = QQPolyRingElem[]
  for i=0:d
    push!(b, (-1)^(d-i)*a[d-i+1]*bin)
    bin *= (t+i+1)*QQFieldElem(1, i+1)
  end
  return sum(b)
end

function Oscar.degree(H::HilbertData)

  @req all(isone, H.weights) "All weights must be 1"

  P = hilbert_polynomial(H)
  if iszero(P)
     q, _ = hilbert_series_reduced(H)
     return q(1)
  end
  deg = leading_coefficient(P)*factorial(ZZ(degree(P)))
  @assert isone(denominator(deg))
  return numerator(deg)
end

function _rational_function_to_power_series(P::QQRelPowerSeriesRing, n, d)
  Qt, t = QQ[:t]
  nn = map_coefficients(QQ, n; parent = Qt)
  dd = map_coefficients(QQ, d; parent = Qt)
  gg, ee, _ = gcdx(dd, t^max_precision(P))
  @assert isone(gg)
  nn = Hecke.mullow(nn, ee, max_precision(P))
  c = collect(coefficients(nn))
  return P(map(QQFieldElem, c), length(c), max_precision(P), 0)
end

function _rational_function_to_power_series(P::QQRelPowerSeriesRing, f)
  return _rational_function_to_power_series(P, numerator(f), denominator(f))
end

@doc raw"""
    expand(f::FracFieldElem{QQPolyRingElem}, d::Int) -> RelPowerSeries

Given a rational function $f$ over the rationals, expand $f$ as a power series
up to terms of degree $d$.

```jldoctest
julia> Qx, x = QQ[:x];

julia> expand(1//(1 - x^2), 5)
1 + t^2 + t^4 + O(t^6)
```
"""
function expand(f::Generic.FracFieldElem{QQPolyRingElem}, d::Int)
  T, t = power_series_ring(QQ, d+1, "t")
  return _rational_function_to_power_series(T, f)
end

function (P::QQRelPowerSeriesRing)(H::HilbertData)
  n, d = hilbert_series_reduced(H)
  return _rational_function_to_power_series(P, n, d)
end

function hilbert_series_expanded(H::HilbertData, d::Int)
   T, t = power_series_ring(QQ, d+1, "t")
   return T(H)
end

function hilbert_function(H::HilbertData, d::Int)
   d < 0 && QQ(0)
   HS = hilbert_series_expanded(H, d)
   return coeff(HS, d)
end

function Base.show(io::IO, h::HilbertData)
  print(io, "Hilbert Series for $(h.I), data: $(h.coeffs), weights: $(h.weights)")  ###new
end

############################################################################
### Homogenization and Dehomogenization
############################################################################

## 2024-02-06  New UI to homogenization & dehomogenization
## used partly the code from 2023-07 (see below)


struct Homogenizer
  P::MPolyRing              # original poly ring (?not graded?)
  P_homog::MPolyDecRing     # graded poly ring: same vars as P, plus g extra homogenizing vars
  VarMap::Vector{Int}       # var[k] in P embeds to var[VarMap[k]] in P_homog
  HVars::Vector{Int}        # the homogenizing vars in p_homog are var[k] with k in HVars
  W::ZZMatrix               # weight matrix for vars in P: k-th column is degree of var[k]
  FastMethodOK::Bool        # true if the weights are "strongly positive", so we can use the fast method

  # STANDARD GRADED case
  function Homogenizer(P::MPolyRing{T}, h::VarName; pos::Int=1+ngens(P))  where T
    n = nvars(P)
    @req pos in 1:n+1 "Index out of range"
    vars = copy(symbols(P))
    insert!(vars, pos, Symbol(h))
    P_homog,_ = graded_polynomial_ring(base_ring(P), vars)
    VarMap = vcat(1:pos-1, pos+1:n+1)
    HVars = [pos]
    W = ZZMatrix(1,n, [1  for i in 1:n])
    return new(P, P_homog, VarMap, HVars, W, true)
  end #function

  # Grading specified by COLUMNS of W
  function Homogenizer(P::MPolyRing{T}, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, h::VarName; pos::Int=1+ngens(P))  where T
    n = nvars(P)
    @req size(W,2) == n  "Weight matrix not compatible with polynomial ring"
    @req pos in 1:n+1  "Index out of range"
    vars = copy(symbols(P))
    grading_dimension = size(W,1)
    if grading_dimension == 1
      insert!(vars, pos, Symbol(h))
    else
      for i = 1:grading_dimension
        insert!(vars, pos-1+i, Symbol("$h[$i]"))  # using names like h[1], h[2], etc
      end
    end
    G = abelian_group(zeros(Int, grading_dimension))
    W = ZZMatrix(W)
    WH = hcat(Matrix(W[:, 1:(pos-1)]),
              Matrix(identity_matrix(ZZ, grading_dimension)),
              Matrix(W[:, pos:n]))
    P_homog,_ = graded_polynomial_ring(base_ring(P), vars; weights = [G(WH[:, i]) for i = 1:size(WH, 2)], cached = false)

    VarMap = vcat(1:pos-1, pos+grading_dimension:n+grading_dimension)
    HVars = collect(pos:pos+grading_dimension-1)
    return new(P, P_homog, VarMap, HVars, W, is_positive_grading_matrix(W))
  end #function

end # struct


# Duplicate constructor but with lower-case initial letter: standard graded case
@doc raw"""
    homogenizer(P::MPolyRing, h::VarName;  pos::Int=1+ngens(P))

Create a "homogenizing operator" assuming a standard grading; `h` is the name of the homogenizing variable; `pos` indicates where to put the homogenizing variable in the list of generators of the graded polynomial ring (default is after all the other variables).

# Examples
```jldoctest
julia> P, (x,y) = polynomial_ring(QQ, [:x, :y]);

julia> H = homogenizer(P, "h");

julia> F = H(x^2+y)
x^2 + y*h

julia> parent(F)
Multivariate polynomial ring in 3 variables over QQ graded by
  x -> [1]
  y -> [1]
  h -> [1]

julia> V = H.([x^2+y, x+y^2]);

julia> parent(V[1]) == parent(V[2])
true

julia> H(ideal([x^2+y]))
Ideal generated by
  x^2 + y*h
```
"""
function homogenizer(P::MPolyRing{T}, h::VarName; pos::Int=1+ngens(P))  where T
  return Homogenizer(P, h; pos=pos)
end

# Duplicate constructor but with lower-case initial letter: W-graded case
@doc raw"""
    homogenizer(P::MPolyRing, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, h::VarName;  pos::Int=1+ngens(P))

Create a "homogenizing operator" using the grading specified by the columns of `W`; `h` is the prefix for the homogenizing variables; `pos` indicates where to put the homogenizing variables in the list of generators of the graded polynomial ring (default is after all the other variables).

# Examples
```jldoctest
julia> P, (x,y) = polynomial_ring(QQ, [:x, :y]);

julia> W = ZZMatrix(2,2, [2,3,5,7]);

julia> H = homogenizer(P, W, "h");

julia> F = H(x^2+y)
x^2 + y*h[1]*h[2]^3

julia> parent(F)
Multivariate polynomial ring in 4 variables over QQ graded by
  x -> [2 5]
  y -> [3 7]
  h[1] -> [1 0]
  h[2] -> [0 1]

julia> V = H.([x^2+y, x+y^2]);

julia> parent(V[1]) == parent(V[2])
true

julia> H(ideal([x^2+y]))
Ideal generated by
  x^2 + y*h[1]*h[2]^3
```
"""
function homogenizer(P::MPolyRing{T}, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, h::VarName; pos::Int=1+ngens(P))  where T
  return Homogenizer(P, W, h; pos=pos)
end


# Returns a list of the homogenizing variables (elements of P_homog)
function homogenizing_variables(H::Homogenizer)
  HVarList = [gen(H.P_homog,j)  for j in H.HVars]
  return HVarList;
end


# Make a Dehomogenizer from a Homogenizer
# (same fields as Homogenizer see comments/doc for struct Homogenizer, above).
# I have chosen to copy all fields so that it would, if necessary, be possible
# to re-create the Homogenizer from its Dehomogenizer.
struct Dehomogenizer
  P::MPolyRing
  P_homog::MPolyDecRing
  VarMap::Vector{Int}
  HVars::Vector{Int}
  W::ZZMatrix
  FastMethodOK::Bool

  function Dehomogenizer(H::Homogenizer)
    return new(H.P, H.P_homog, H.VarMap, H.HVars, H.W)
  end #function
end # struct


# Duplicate constructor but with lower-case initial letter
@doc raw"""
    dehomogenizer(H::Homogenizer)

Create a "dehomogenizing operator" from a `Homogenizer`; it is effectively a polynomial ring homomorphism mapping all homogenizing variables to 1.  A `Dehomogenizer` is a post-inverse for the `Homogenizer` it was created from.

# Examples
```jldoctest
julia> P, (x,y) = polynomial_ring(QQ, [:x, :y]);

julia> H = homogenizer(P, "h");

julia> DH = dehomogenizer(H);

julia> F = H(x^2+y)
x^2 + y*h

julia> DH(F)
x^2 + y

julia> V = [x^2+y, x*y+y^2]; HV = H.(V);

julia> parent(DH(HV[1])) == P  &&  parent(DH(HV[2])) == P
true

julia> DH(H(ideal(V)))
Ideal generated by
  x*y + y^2
  x^2 + y
  y^3 + y^2
```
"""
function dehomogenizer(H::Homogenizer)
  return Dehomogenizer(H)
end


# Use a homogenizer as a function/operator to actually homogenize
function (H::Homogenizer)(f::MPolyRingElem)
  @req parent(f) == H.P  "Not in correct polynomial ring"
  return Oscar._homogenization(f, H.P_homog, H.HVars[1])
end

function (H::Homogenizer)(I::MPolyIdeal; extra_gens_flag::Bool = false)
  return _homogenization(I, H; extra_gens_flag)
end


# Use a dehomogenizer as a function/operator to actually dehomogenize
function (DH::Dehomogenizer)(f::MPolyDecRingElem)
  @req parent(f) == DH.P_homog   "Not in correct graded polynomial ring"
  return Oscar._dehomogenization(f, DH.P, DH.HVars[1], length(DH.HVars))
end

function (DH::Dehomogenizer)(I::MPolyIdeal{MPolyDecRingElem{T1,T2}})  where T1  where T2
  @req base_ring(I) == DH.P_homog  "Dehomogenizer not compatible with ideal"
  return ideal(DH.(gens(I)));
end


# ==================================================================
# 2023-06-30 START: New impl of homogenization
# By John Abbott, based on K+R Book "Computational Commutative Algebra"
# This code delegates to _homogenization_via_saturation for non-positive gradings
# ==================================================================

## ------------------------------------------------------------------
## INTERNAL FUNCTION:  _homogenization_via_saturation
## Homogenization of an ideal -- works for all gradings (incl. non-positive),
## but is slow.  If grading is positive over ZZ^1: use Singular (and wdegrevlex);
## if grading is positive over ZZ^m (with m > 1) use new "fast/clever" code.
## Author: John Abbott  2023-07-14
## ------------------------------------------------------------------

# Used only in _homogenization_via_saturation (immediately below).
# This function returns gens(I) and possibly some more "small" elements of I
# (obtained from 1 or more groebner bases of I).  If extra_gens_flag is true,
# we compute one more groebner bases of I in the hope that they contain some small elements.
function _gens_for_homog_via_sat(I::MPolyIdeal{T}, extra_gens_flag::Bool) where {T <: MPolyRingElem}
  G = gens(I)
  if isempty(G)  throw("Ideal must have at least 1 generator"); end;
  OrigR = parent(G[1])      # Since I is not zero, it has at least 1 gen.
  # Next few lines: we adjoin some more small, redundant gens.
  # !!HEURISTIC!!  We use 2*AveNumTerms as size limit for redundant gens we shall adjoin.
  AveNumTerms = sum([length(g)  for g in G])/length(G) ## floating-point!
  if extra_gens_flag
    # compute DegRevLex GB (?and maybe DegLex GB?)
    assume_gb_is_cached = groebner_basis(I; ordering=degrevlex(OrigR))
    #??Good idea??        assume_gb_is_cached = groebner_basis(I; ordering=deglex(OrigR))
  end
  for GB in values(I.gb)
    extra_gens = filter(f -> (length(f) < 2*AveNumTerms), elements(GB))
    G = vcat(G,  extra_gens)
  end
  return G
end

# Fully general, but typically much slower than specialized method below (for positive gradings)
# Optional kwarg:
#   extra_gens_flag: if false, inhibits computing grobner basis speculatively in _gens_for_homog_via_sat
function _homogenization_via_saturation(I::MPolyIdeal{T},  H::Homogenizer; extra_gens_flag::Bool = false) where {T <: MPolyRingElem}
  P = base_ring(I)
  @req P == H.P  "Homogenizer not for this ring"
  if is_zero(I)  # special handling for ideal(0)
    return ideal([zero(H.P_homog)]) # return zero ideal in !!P_homog!!
  end
  Hgens = H.(_gens_for_homog_via_sat(I, extra_gens_flag))
  if length(Hgens) == 1  # short-cut for principal ideal
    return ideal(Hgens)
  end
  DoSatByProduct = false  # true means saturate by product of homogenizing vars;
                          # false means saturate successively by each homogenizing var (faster in most cases?)
  HVars = homogenizing_variables(H)
  if DoSatByProduct
    Ih = saturation(ideal(Hgens), ideal(prod(HVars)))
  else # not DoSatByProduct, so do a cascade of saturations
    Ih = ideal(Hgens)
    for hvar in HVars
      Ih = saturation(Ih, ideal([hvar]))
      # IMPROVE???  if length(gens(Ih)) == 1  return Ih;  end
    end
  end
  return Ih
end



# Saturate a polynomial by a variable -- equiv to  f/gcd(f,h^deg(f))
function _sat_poly_by_var(f::MPolyRingElem{T}, h::MPolyRingElem{T})  where { T <: RingElement }
  # REQUIRE h is a variable i.e. monic poly of deg 1 with just 1 term
  @req parent(f) == parent(h)  "Polynomial and variable must be in same ring"
  @req (#=is_monic(h) &&=# length(h) == 1 && total_degree(h) == 1)  "Arg 2 must be a variable"
  if is_zero(f)
    return f;
  end
  # Below i is index of the variable h
  i = findfirst(>(0), exponent_vector(h,1))
  n = length(f)
  EV = AbstractAlgebra.exponent_vectors(f)  # list or iterator
  multiplicity = exponent_vector(f,1)[i]
  # Loop below computes min exponent of h in the terms of f
  for expv in EV
    multiplicity = min(multiplicity, expv[i])
    if is_zero(multiplicity)
      return f
    end
  end
  return f/h^multiplicity;
end

# This function should be in a file of utilities (JAA thinks)
function kronecker_delta(i::Int, j::Int)
  return (i == j) ? 1 : 0
end


# Check that W defines a positive grading (defn 4.2.4 from R+K vol 2): namely
#     each col contains a non-zero entry, and
#     first non-zero going down the col is positive.
function is_positive_grading_matrix(W::Union{ZZMatrix, Matrix{<:IntegerUnion}})
  nrows = size(W,1)
  ncols = size(W,2)
  for c in 1:ncols
    IsGoodCol = false
    for r in 1:nrows
      if is_zero(W[r,c])
        continue
      end
      if W[r,c] < 0
        return false
      end
      IsGoodCol = true
      break
    end
    if !IsGoodCol
      return false
    end
  end
  return true
end


# --------------------------------------------
# This homogenization impl seems to be usefully faster with positive ZZ^m-grading for m > 1

# This is the MAIN WORKER FUNCTION (for homogenizing ideals).  It has an optional kwarg
# (*) extra_gens_flag (relevant only for non-positive gradings); default is false, but
#     if true it triggers computation of a "needless" groebner basis in function
#     _gens_for_homog_via_sat which can overall speed up homogenization.
function _homogenization(I::MPolyIdeal{T}, H::Homogenizer; extra_gens_flag::Bool = false)  where {T <: MPolyRingElem}
  ## ASSUME NumCols(W) = nvars(P)
  P = base_ring(I) # initial polynomial ring
  @req P == H.P  "Homogenizer not for the ring of this ideal"
  # Handle zero ideal as special case: (quick and easy)
  if  is_zero(I)
    return ideal(H.P_homog, [])  # ideal(0) in correct ring
  end
  if length(gens(I)) == 1
    return ideal([H(gen(I,1))])
  end
  # The fast method below is valid only for positive gradings;
  # so if not positive graded, delegate to _homogenization_via_saturation (general but slower)
  if !H.FastMethodOK
    return _homogenization_via_saturation(I, H; extra_gens_flag=extra_gens_flag)
  end
  # >>> POSITIVELY GRADED <<<  case henceforth
  if is_z_graded(H.P_homog)
    # Case: grading is over ZZ^1, so delegate to faster special-case function.
    # For correctness see Kreuzer+Robbiano (vol.2) Tutorial 53 part (d).
    if all(x -> (degree(x)[1] == 1), gens(H.P_homog))
      W_ordering = degrevlex(P)
    else
      W_ordering = wdegrevlex(P, [Int(H.W[1,j])  for j in 1:nvars(P)]) # or get columns from H.W ???
    end
    GB = groebner_basis(I; ordering=W_ordering)
    return ideal(H.(elements(GB)))
  end
  # Case: positive not ZZ^1-graded
  # Handle one ideal as special case (NB is_one might be costly!)
  if  is_one(I)
    return ideal([one(H.P_homog)])  # ideal(1) in correct ring
  end
  # Case: positive grading over ZZ^k with k > 1
  # >>> JUSTIFICATION <<<  from book Kreuzer+Robbiano "Computational Commutative Algebra"
  # Combination of Cor 4.3.8(a) [main algorithm],  Cor 4.4.9 [saturation by indet], and
  # Tutorial 37(g) [compute saturation by product via "cascade"]
  G = gens(I);
  if length(G) == 1  # short-cut for principal ideal
    return ideal(H.(G));
  end
  grading_dimension = nrows(H.W)   ###ngens(H.P_homog) - ngens(H.P);
  H2 = homogenizer(P, H.W, "hh"; pos = 1+ngens(P)) # hh vars come last -- different from homogenizer H
  G = H2.(G)
  Ph = parent(G[1])  # equiv H2.P_homog  # Ph is graded ring, "homog-extn" of P.
  N = ngens(Ph)
  num_x = ngens(P)
  # Build ordering matrix: weights matrix followed by identity mat, underneath is a invlex matrix
  Id = reduce(hcat, [[kronecker_delta(i,j)  for i in 1:grading_dimension]  for j in 1:grading_dimension])
  RevLexMat = reduce(hcat, [[-kronecker_delta(i+j, 1+N)  for i in 1:N]  for j in 1:N])
  W_int = reduce(hcat, [[Int(H.W[i,j])  for i in 1:grading_dimension]  for j in 1:nvars(P)])
  M = hcat(W_int, Id)
  for j in 1:grading_dimension
    M_complete = vcat(M, RevLexMat)
    JohnsOrdering = matrix_ordering(Ph, M_complete)
    Ih = ideal(G)
    G = groebner_basis(Ih; ordering=JohnsOrdering)  # option complete_reduction seemed to make no measurable difference
    G = elements(G)
    h_j = gen(Ph, ngens(Ph)+1-j)
    G = [_sat_poly_by_var(g, h_j)  for g in G]
    # Loop below: rotate cols in RevLexMat which correspond to the "hh" variables:
    for i in 1:grading_dimension
      col1 = (i+j > 1+grading_dimension) ? (N+2+grading_dimension-i-j) : (N+2-i-j)
      col2 = (col1 == N-grading_dimension+1) ? N : col1-1
      RevLexMat[i, col1] = 0
      RevLexMat[i, col2] = -1
    end
  end
  # Now map result into correct ring (?? suffice to do H.(DH2.(G)) ??)
  P_homog = H.P_homog
  pos = H.HVars[1]
    images = [zero(P_homog)  for i in 1:N]
    for i in 1:pos-1
      images[i] = gen(P_homog,i)
    end
    for i in pos:num_x
      images[i] = gen(P_homog,i+grading_dimension)
    end
    for i in 1:grading_dimension
      images[num_x+i] = gen(P_homog, pos+i-1)
    end
    reorder_variables = hom(Ph, P_homog, images)
  G = reorder_variables.(G)
  return ideal(G) # now in correct ring
end

# ============================================
# 2024-02-06 END: New UI for homogenization
# ============================================


# This function is old (not sure who wrote it)
function _homogenization(f::MPolyRingElem, S::MPolyDecRing, start_pos::Int)
    # ASSUME start_pos is in 1:1+ngens(parent(f))
   len_gg  = ngens(S.D)
      #= get exponent vectors and enlarge them for S =#
   exps = Vector{Int}[]
   for e in AbstractAlgebra.exponent_vectors(f)
       for j in 1:len_gg
           insert!(e, start_pos+j-1, 0)
       end
       push!(exps, e)
   end

   #= compute degrees in S =#
   if is_z_graded(S) == true
       degs    = Int[]
       gens_gg = [Int(S.d[i][1])  for i in 1:ngens(S)]
   else
       degs    = Vector{Int}[]
       gens_gg = [[Int(S.d[i][j]) for j in 1:ngens(S.D)] for i in 1:ngens(S)]
   end

   for e in exps
       push!(degs, sum((e.*gens_gg)))
   end
   # get top degrees
   tdeg  = vcat(maximum(hcat(degs...); dims=2)...)
       #= finally generate homogeneous version of f in S, called F =#
   F  = MPolyBuildCtx(S)

   for (i, (c, e)) = enumerate(zip(AbstractAlgebra.coefficients(f), exps))
       for j in 1:len_gg
           e[start_pos+j-1] = tdeg[j]-degs[i][j]
       end
       push_term!(F, c, e)
   end

   return finish(F)
end

# ============================================
# 2023-06-30 END: New impl of homogenization
# ============================================



# This is also old code: if it ain't broke, don't fix it!

function _dehomogenization(F::MPolyDecRingElem, R::MPolyRing, pos::Int, m::Int)
  B = MPolyBuildCtx(R)
  for (c,e) = zip(AbstractAlgebra.coefficients(F), AbstractAlgebra.exponent_vectors(F))
    deleteat!(e, pos:pos+m-1)
    push_term!(B, c, e)
  end
  return finish(B)
end



################################################################################
#
#  Evaluation
#
################################################################################

(f::MPolyDecRingElem)(x...) = evaluate(f, collect(x))

################################################################################
#
#  Promote rule
#
################################################################################

function AbstractAlgebra.promote_rule(::Type{MPolyDecRingElem{S, T}}, ::Type{U}) where {S, T, U}
  if AbstractAlgebra.promote_rule(T, U) === T
    return MPolyDecRingElem{S, T}
  else
    return Union{}
  end
end

function AbstractAlgebra.promote_rule(::Type{MPolyDecRingElem{S, T}}, ::Type{MPolyDecRingElem{S, T}}) where {S, T}
  return MPolyDecRingElem{S, T}
end

################################################################################
#
#  Random homogeneous polynomials
#
################################################################################

#create homogeneous polynomials in graded rings
#TODO: make this work for non-standard gradings
@doc raw"""
    rand(S::MPolyDecRing, term_range, deg_range, v...)

Create a random homogeneous polynomial with a random number of
terms (`rand(term_range)`) and of random degree (`rand(deg_range)`)
and random coefficients via `v...`.
"""
function rand(S::MPolyDecRing, term_range, deg_range, v...)
  f = zero(forget_decoration(S))
  d = rand(deg_range)
  for i=1:rand(term_range)
    t = forget_decoration(S)(rand(base_ring(S), v...))
    if iszero(t)
      continue
    end
    for j=1:ngens(forget_decoration(S))-1
      t *= gen(forget_decoration(S), j)^rand(0:(d-total_degree(t)))
    end
    #the last exponent is deterministic...
    t *= gen(forget_decoration(S), ngens(forget_decoration(S)))^(d-total_degree(t))
    f += t
  end
  return S(f)
end

################################################################################
#
#  Forgetting the decoration / grading
#
################################################################################

@doc raw"""
    forget_decoration(R::MPolyDecRing)

Return the underlying undecorated ring.
"""
forget_decoration(R::MPolyDecRing) = R.R

@doc raw"""
    forget_grading(R::MPolyDecRing)

Return the underlying ungraded undecorated ring.
"""
forget_grading(R::MPolyDecRing) = forget_decoration(R)

@doc raw"""
    forget_decoration(f::MPolyDecRingElem)

Return the element in the underlying undecorated ring.
"""
forget_decoration(f::MPolyDecRingElem) = f.f

@doc raw"""
    forget_grading(f::MPolyDecRingElem)

Return the element in the underlying ungraded undecorated ring.
"""
forget_grading(f::MPolyDecRingElem) = forget_decoration(f)

@doc raw"""
    forget_decoration(I::MPolyIdeal{<:MPolyDecRingElem})

Return the ideal in the underlying undecorated ring.
"""
function forget_decoration(I::MPolyIdeal{<:MPolyDecRingElem})
  R = forget_decoration(base_ring(I))
  return ideal(R, map(forget_decoration, gens(I)))
end

@doc raw"""
    forget_grading(I::MPolyIdeal{<:MPolyDecRingElem})

Return the ideal in the underlying ungraded undecorated ring.
"""
forget_grading(I::MPolyIdeal{<:MPolyDecRingElem}) = forget_decoration(I)

Generic.internal_ordering(S::MPolyDecRing) = internal_ordering(S.R)

#############truncation#############

@doc raw"""
    truncate(I::MPolyIdeal, g::FinGenAbGroupElem)

Given a (homogeneous) ideal `I` in a $\mathbb Z$-graded multivariate polynomial ring
with positive weights, return the truncation of `I` at degree `g`.

    truncate(I::MPolyIdeal, d::Int)

Given an ideal `I` as above, and given an integer `d`, convert `d`
into an element `g` of the grading group of `base_ring(I)` and proceed as above.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x, y^4, z^6])
Ideal generated by
  x
  y^4
  z^6

julia> truncate(I, 3)
Ideal generated by
  x*z^2
  x*y*z
  x*y^2
  x^2*z
  x^2*y
  x^3
  y^4
  z^6
```

```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [3,2,1]);

julia> I = ideal(R, [x, y^4, z^6])
Ideal generated by
  x
  y^4
  z^6

julia> truncate(I, 3)
Ideal generated by
  x
  y^4
  z^6

julia> truncate(I, 4)
Ideal generated by
  x*z
  z^6
  y^4
```
"""
function truncate(I::MPolyIdeal, g::FinGenAbGroupElem)
  return truncate(I, Int(g[1]))
end

function  truncate(I::MPolyIdeal, d::Int)
  R = base_ring(I)
  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  @req is_z_graded(R) "The base ring must be ZZ-graded"
  W = R.d
  W = [Int(W[i][1]) for i = 1:ngens(R)]
  @req minimum(W) > 0 "The weights must be positive"
  if is_zero(I)
     return I
  end
  dmin = minimum(degree(Int, x) for x in gens(I))
  if  d <= dmin
     return I
  end
  V = sort(gens(I); by=a -> degree(Int, a))
  RES = elem_type(R)[]
  s = dmin
  B = monomial_basis(R, d-s)
  for i = 1:length(V)
       if degree(Int, V[i]) < d
          if degree(Int, V[i]) > s
             s = degree(Int, V[i])
             B = monomial_basis(R, d-s)
          end
          append!(RES, B .*V[i])
       else
           push!(RES, V[i])
       end
  end
  return ideal(R, RES)
end

################################################################################
#
# Minimal generating set
#
################################################################################

@doc raw"""
    minimal_generating_set(I::MPolyIdeal{<:MPolyDecRingElem})

Given a (homogeneous) ideal `I` in a graded multivariate polynomial ring
over a field, return an array containing a minimal set of generators
of `I`. If `I` is the zero ideal, an empty list is returned.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> V = [x, z^2, x^3+y^3, y^4, y*z^5];

julia> I = ideal(R, V)
Ideal generated by
  x
  z^2
  x^3 + y^3
  y^4
  y*z^5

julia> minimal_generating_set(I)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x
 z^2
 y^3

julia> I = ideal(R, zero(R))
Ideal generated by
  0

julia> minimal_generating_set(I)
MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[]
```
"""
function minimal_generating_set(I::MPolyIdeal{<:MPolyDecRingElem})
  # This only works / makes sense for homogeneous ideals. So far ideals in an
  # MPolyDecRing are forced to be homogeneous though.

  R = base_ring(I)

  @assert is_graded(R)

  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"

  if !isempty(I.gb)
    # make sure to not recompute a GB from scratch on the singular
    # side if we have one
    G = first(values(I.gb))
    _, sing_min = Singular.mstd(singular_generators(G, G.ord))
    return filter(!iszero, (R).(gens(sing_min)))
  else
    sing_gb, sing_min = Singular.mstd(singular_generators(I))
    ring = base_ring(I)
    computed_gb = IdealGens(ring, sing_gb, true)
    I.gb[computed_gb.ord] = computed_gb
    return filter(!iszero, (R).(gens(sing_min)))
  end
end

##################regularity#######################

@doc raw"""
    cm_regularity(I::MPolyIdeal)

Given a (homogeneous) ideal `I` in a standard $\mathbb Z$-graded multivariate polynomial ring
with coefficients in a field, return the Castelnuovo-Mumford regularity of I.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> I = ideal(R, [y^2*z − x^2*w, z^4 − x*w^3]);

julia> cm_regularity(I)
6

julia> minimal_betti_table(I)
degree: 0  1
------------
     3: 1  -
     4: 1  -
     5: -  -
     6: -  1
------------
 total: 2  1
```
"""
function cm_regularity(I::MPolyIdeal)
   @req is_standard_graded(base_ring(I)) "The base ring is not standard ZZ-graded"
   B = minimal_betti_table(I)
   S = as_dictionary(B)
   V = [x[2][1] - x[1] for x in keys(S)]
  return maximum(V)
end

#######################################################
#######################################################
@doc raw"""
    degree(I::MPolyIdeal)

Given a (homogeneous) ideal `I` in a standard $\mathbb Z$-graded multivariate polynomial ring,
return the degree of `I` (that is, the degree of the quotient of `base_ring(I)` modulo `I`).
Otherwise, return the degree of the homogenization of `I` with respect to the
standard $\mathbb Z$-grading.

!!! note
    Geometrically, the degree of a homogeneous ideal as above is the number
    of intersection points of its projective variety with a generic linear
    subspace of complementary dimension (counted with multiplicities).
    See also [MS21](@cite).

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> I = ideal(R, [y-x^2, x-z^3])
Ideal generated by
  -x^2 + y
  x - z^3

julia> degree(I)
6
```
"""
@attr ZZRingElem function degree(I::MPolyIdeal)
  P = base_ring(I)
  if !is_standard_graded(P)
    H = homogenizer(P, "_h")
    I = H(I)
  end
  A, _ = quo(base_ring(I), I)
  return degree(A)
end

