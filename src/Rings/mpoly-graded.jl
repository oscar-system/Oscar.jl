
@attributes mutable struct MPolyDecRing{T, S} <: AbstractAlgebra.MPolyRing{T}
  R::S
  D::GrpAbFinGen
  d::Vector{GrpAbFinGenElem}
  lt
  function MPolyDecRing(R::S, d::Vector{GrpAbFinGenElem}) where {S}
    @assert length(d) == ngens(R)
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    return r
  end
  function MPolyDecRing(R::S, d::Vector{GrpAbFinGenElem}, lt) where {S}
    @assert length(d) == ngens(R)
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    r.lt = lt
    return r
  end
end

@doc raw"""
    grading_group(R::MPolyDecRing)

If `R` is, say, `G`-graded, then return `G`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> grading_group(R)
GrpAb: Z
```
"""
grading_group(R::MPolyDecRing) = R.D

function is_graded(R::MPolyDecRing)
   return !isdefined(R, :lt)
end

@doc raw"""
    is_graded(R::MPolyRing)

Return `true` if `R` is graded, `false` otherwise.
"""
is_graded(R::MPolyRing) = false

is_filtered(R::MPolyDecRing) = isdefined(R, :lt)
is_filtered(::MPolyRing) = false

function show(io::IO, ::MIME"text/plain", W::MPolyDecRing)
  Hecke.@show_name(io, W)
  Hecke.@show_special(io, W)
  io = pretty(io)
  R = forget_decoration(W)
  print(io, R)
  if is_filtered(W)
    println(io, " filtrated by ")
  else
    println(io, " graded by ")
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
  print(io, Dedent())
#  println(IOContext(io, :compact => true, ), W.d)
end

function Base.show(io::IO, W::MPolyDecRing)
  Hecke.@show_name(io, W)
  Hecke.@show_special(io, W)
  io = pretty(io)
  if get(io, :supercompact, false)
    # no nested printing
    if is_filtered(W)
      print(io, "Filtered multivariate polynomial ring")
    else
      print(io, "Graded multivariate polynomial ring")
    end
  else
    # nested printing allowed, preferably supercompact
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
    grade(R::MPolyRing, W::Vector{<:IntegerUnion})

Given a vector `W` of `ngens(R)` integers, create a free abelian group of type `GrpAbFinGen` 
given by one free generator, and convert the entries of `W` to elements of that group. Then 
create a $\mathbb Z$-graded ring by assigning the group elements as weights to the variables 
of `R`, and return the new ring, together with the vector of variables.

    grade(R::MPolyRing)

As above, where the grading is the standard $\mathbb Z$-grading on `R`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> T, (x, y, z) = grade(R)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])
```
"""
function grade(R::MPolyRing, W::Vector{<:IntegerUnion})
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
    grade(R::MPolyRing, W::Vector{<:Vector{<:IntegerUnion}})
 
Given a vector `W` of `ngens(R)` integer vectors of the same size `m`, say, create a free 
abelian group of type `GrpAbFinGen` given by `m` free generators, and convert the vectors in
`W` to elements of that group. Then create a $\mathbb Z^m$-graded ring by assigning the group 
elements as weights to the variables of `R`, and return the new ring, together with the vector
of variables.

    grade(R::MPolyRing, W::Union{ZZMatrix, Matrix{<:IntegerUnion}})

As above, converting the columns of `W`.

# Examples
```jldoctest
julia> R, x, y = polynomial_ring(QQ, "x" => 1:2, "y" => 1:3)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2]], QQMPolyRingElem[y[1], y[2], y[3]])

julia> W = [1 1 0 0 0; 0 0 1 1 1]
2×5 Matrix{Int64}:
 1  1  0  0  0
 0  0  1  1  1

julia> grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], y[1], y[2], y[3]])
```
"""
function grade(R::MPolyRing, W::Vector{<:Vector{<:IntegerUnion}})
  @assert length(W) == ngens(R)
  n = length(W[1])
  @assert all(x->length(x) == n, W)

  A = abelian_group(zeros(Int, n))
  set_attribute!(A, :show_elem => show_special_elem_grad) 
  S = MPolyDecRing(R, [A(w) for w = W])
  return S, map(S, gens(R))
end

function grade(R::MPolyRing, W::Union{ZZMatrix, Matrix{<:IntegerUnion}})
  @assert size(W, 2) == ngens(R)
  A = abelian_group(zeros(Int, size(W, 1)))
  set_attribute!(A, :show_elem => show_special_elem_grad) 
  ###S = MPolyDecRing(R, [A(view(W, :, i)) for i = 1:size(W, 2)])
  S = MPolyDecRing(R, [A(W[:, i]) for i = 1:size(W, 2)])
  return S, map(S, gens(R))
end



@doc raw"""
    is_standard_graded(R::MPolyDecRing)

Return `true` if `R` is standard $\mathbb Z$-graded, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
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
    Writing `G = grading_group(R)`, we say that `R` is $\mathbb Z$-graded if `is_free(G) && ngens(G) == rank(G) == 1` evaluates to `true`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> is_z_graded(S)
true
```
"""
function is_z_graded(R::MPolyDecRing)
  is_graded(R) || return false
  A = grading_group(R)
  return ngens(A) == 1 && rank(A) == 1 && is_free(A)
end

is_z_graded(::MPolyRing) = false

@doc raw"""
    is_zm_graded(R::MPolyDecRing)

Return `true` if `R` is $\mathbb Z^m$-graded for some $m$, `false` otherwise.

!!! note
    Writing `G = grading_group(R)`, we say that `R` is $\mathbb Z^m$-graded if `is_free(G) && ngens(G) == rank(G) == m` evaluates to `true`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, "x" => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> is_zm_graded(S)
false

julia> G = abelian_group(ZZMatrix([1 -1]));

julia> g = gen(G, 1)
Element of G with components [0 1]

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], W);

julia> is_free(G)
true

julia> is_zm_graded(R)
false
```
"""
function is_zm_graded(R::MPolyDecRing)
  is_graded(R) || return false
  A = grading_group(R)
  return is_free(A) && ngens(A) == rank(A)
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
julia> R, (t, x, y) = polynomial_ring(QQ, ["t", "x", "y"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[t, x, y])

julia> G = abelian_group([0])
GrpAb: Z

julia> S, (t, x, y) = grade(R, [-1, 1, 1])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[t, x, y])

julia> is_positively_graded(S)
false

julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> G = abelian_group([0, 2])
(General) abelian group with relation matrix
[0 0; 0 2]

julia> W = [gen(G, 1)+gen(G, 2), gen(G, 1)]
2-element Vector{GrpAbFinGenElem}:
 Element of G with components [1 1]
 Element of G with components [1 0]

julia> S, (x, y) = grade(R, W)
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
  if ngens(G) == rank(G)
    W = vcat([x.coeff for x = R.d])
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
    graded_polynomial_ring(C::Ring, V, W; ordering=:lex)
    graded_polynomial_ring(C::Ring, n::Int, s::T, W; ordering=:lex) where T<:VarName

Create a multivariate [`polynomial_ring`](@ref polynomial_ring(R, [:x])) with
coefficient ring `C` and variables which print according to the variable names
in `V` (respectively as "\$(s)1" up to "\$s\$n"), and [`grade`](@ref) this ring
according to the data provided by `W`. Return the graded ring as an object of
type `MPolyDecRing`, together with the vector of variables.

    graded_polynomial_ring(C::Ring, V; ordering=:lex)

As above, where the grading is the standard $\mathbb Z$-grading.

# Examples
```jldoctest
julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 0]
 [4, 1]

julia> R, x = graded_polynomial_ring(QQ, 4, :x, W)
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x1, x2, x3, x4])

julia> S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])
```
"""
function graded_polynomial_ring(C::Ring, V::Union{Tuple{Vararg{T}}, AbstractVector{T}},
      W=ones(Int, length(V)); ordering=:lex, cached::Bool=true) where
      T<:VarName
   # HACK: we pass 'cached' only to 'polynomial_ring', but it really should
   # also affect whether the MPolyDecRing gets cached...
   return grade(polynomial_ring(C, V; ordering, cached)[1], W)
end

function graded_polynomial_ring(C::Ring, n::Int, s::VarName,
      W=ones(Int, n); ordering=:lex, cached::Bool=true)
   # HACK: we pass 'cached' only to 'polynomial_ring', but it really should
   # also affect whether the MPolyDecRing gets cached...
   return grade(polynomial_ring(C, n, s; ordering, cached)[1], W)
end

filtrate(R::MPolyRing) = decorate(R)

function show_special_elem_grad(io::IO, a::GrpAbFinGenElem)
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

function filtrate(R::MPolyRing, v::Vector{GrpAbFinGenElem}, lt)
  S = MPolyDecRing(R, v, lt)
  return S, map(S, gens(R))
end

@doc raw"""
    grade(R::MPolyRing, W::Vector{GrpAbFinGenElem})

Given a vector `W` of `ngens(R)` elements of a finitely presented group `G`, say, create a 
`G`-graded ring by assigning the entries of `W` as weights to the variables of `R`. Return
the new ring as an object of type `MPolyDecRing`, together with the vector of variables.

# Examples
```jldoctest
julia> R, (t, x, y) = polynomial_ring(QQ, ["t", "x", "y"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[t, x, y])

julia> typeof(R)
QQMPolyRing

julia>  typeof(x)
QQMPolyRingElem

julia> G = abelian_group([0])
GrpAb: Z

julia> g = gen(G, 1)
Element of G with components [1]

julia> S, (t, x, y) = grade(R, [-g, g, g])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[t, x, y])

julia> typeof(S)
MPolyDecRing{QQFieldElem, QQMPolyRing}

julia> S isa MPolyRing
true

julia> typeof(x)
MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}

julia> R, x, y = polynomial_ring(QQ, "x" => 1:2, "y" => 1:3)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2]], QQMPolyRingElem[y[1], y[2], y[3]])

julia> G = abelian_group([0, 0])
GrpAb: Z^2

julia> g = gens(G)
2-element Vector{GrpAbFinGenElem}:
 Element of G with components [1 0]
 Element of G with components [0 1]

julia> W = [g[1], g[1], g[2], g[2], g[2]];

julia> S, _ = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], y[1], y[2], y[3]])

julia> typeof(x[1])
QQMPolyRingElem

julia> x = map(S, x)
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x[1]
 x[2]

julia> y = map(S, y)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 y[1]
 y[2]
 y[3]

julia> typeof(x[1])
MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}

julia> R, x = polynomial_ring(QQ, "x" => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G)
4-element Vector{GrpAbFinGenElem}:
 Element of G with components [1 0 0 0]
 Element of G with components [0 1 0 0]
 Element of G with components [0 0 1 0]
 Element of G with components [0 0 0 1]

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])
```
"""
function grade(R::MPolyRing, v::Vector{GrpAbFinGenElem})
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

Nemo.symbols(R::MPolyDecRing) = symbols(forget_decoration(R))
Nemo.nvars(R::MPolyDecRing) = nvars(forget_decoration(R))

elem_type(::MPolyDecRing{T, S}) where {T, S} = MPolyDecRingElem{T, elem_type(S)}
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

for T in [:(+), :(-), :(*), :divexact]
  @eval ($T)(a::MPolyDecRingElem,
             b::MPolyDecRingElem) = MPolyDecRingElem($T(forget_decoration(a), forget_decoration(b)), parent(a))
end

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

divexact(a::MPolyDecRingElem, b::RingElem) = MPolyDecRingElem(divexact(forget_decoration(a), b), parent(a))

divexact(a::MPolyDecRingElem, b::Integer) = MPolyDecRingElem(divexact(forget_decoration(a), b), parent(a))

divexact(a::MPolyDecRingElem, b::Rational) = MPolyDecRingElem(divexact(forget_decoration(a), b), parent(a))

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

function factor(x::Oscar.MPolyDecRingElem)
  R = parent(x)
  D = Dict{elem_type(R), Int64}()
  F = factor(forget_decoration(x))
  n=length(F.fac)
  #if n == 1
  #  return Fac(R(F.unit), D)
  #else
    for i in keys(F.fac)
     push!(D, R(i) => Int64(F[i]))
    end
  return Fac(R(F.unit), D)
  #end
end

function gcd(x::Oscar.MPolyDecRingElem, y::Oscar.MPolyDecRingElem)
  R = parent(x)
  return R(gcd(forget_decoration(x), forget_decoration(y)))
end

function div(x::Oscar.MPolyDecRingElem, y::Oscar.MPolyDecRingElem)
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

function mul!(a::MPolyDecRingElem, b::MPolyDecRingElem, c::MPolyDecRingElem)
  return b*c
end

function addeq!(a::MPolyDecRingElem, b::MPolyDecRingElem)
  return a+b
end

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
  w_ord = degrevlex(gens(R)) # dummy, not used
  # This is not meant to be exhaustive, there a probably more gradings which one
  # can meaningfully translate into a monomial ordering
  # However, we want to stick to global orderings.
  if is_z_graded(R)
    w = Int[ R.d[i].coeff[1] for i = 1:ngens(R) ]
    if all(isone, w)
      w_ord = degrevlex(gens(R))
      grading_to_ordering = true
    elseif all(x -> x > 0, w)
      w_ord = wdegrevlex(gens(R), w)
      grading_to_ordering = true
    end
  end
  return grading_to_ordering, w_ord
end

function default_ordering(R::MPolyDecRing)
  return get_attribute!(R, :default_ordering) do
    fl, w_ord = has_weighted_ordering(R)
    fl ? w_ord : degrevlex(gens(R))
  end
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

function jacobi_matrix(f::MPolyDecRingElem)
  R = parent(f)
  n = nvars(R)
  return matrix(R, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_ideal(f::MPolyDecRingElem)
  R = parent(f)
  n = nvars(R)
  return ideal(R, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyDecRingElem})
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
julia> R, x = polynomial_ring(QQ, "x" => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[2]^2+2*x[4]^2
x[2]^2 + 2*x[4]^2

julia> degree(f)
Element of G with components [0 2 0 0]

julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 0]
 [4, 1]

julia> R, x = graded_polynomial_ring(QQ, ["x[1]", "x[2]", "x[3]", "x[4]"], W)
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^4*x[2]+x[4]
x[1]^4*x[2] + x[4]

julia> degree(f)
[4 1]

julia> degree(Vector{Int}, f)
2-element Vector{Int64}:
 4
 1

julia>  R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^6+y^3+z^2
x^6 + y^3 + z^2

julia> degree(f)
[6]

julia> typeof(degree(f))
GrpAbFinGenElem

julia> degree(Int, f)
6

julia> typeof(degree(Int, f))
Int64
```
"""
function degree(a::MPolyDecRingElem)
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

function degree(::Type{Int}, a::MPolyDecRingElem)
  @assert is_z_graded(parent(a))
  return Int(degree(a)[1])
end

function degree(::Type{Vector{Int}}, a::MPolyDecRingElem)
  @assert is_zm_graded(parent(a))
  d = degree(a)
  return Int[d[i] for i=1:ngens(parent(d))]
end

@doc raw"""
    is_homogeneous(f::MPolyDecRingElem)

Given an element `f` of a graded multivariate ring, return `true` if `f` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^2+y*z
x^2 + y*z

julia> is_homogeneous(f)
false

julia> W = [1 2 1 0; 3 4 0 1]
2×4 Matrix{Int64}:
 1  2  1  0
 3  4  0  1

julia> S, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], W)
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
  S = Set{elem_type(D)}()
  for c = MPolyExponentVectors(forget_decoration(F))
    u = parent(F).D[0]
    for i=1:length(c)
      u += c[i]*d[i]
    end
    push!(S, u)
    if length(S) > 1
      return false
    end
  end
  return true
end

@doc raw"""
    homogeneous_components(f::MPolyDecRingElem{T, S}) where {T, S}

Given an element `f` of a graded multivariate ring, return the homogeneous components of `f`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^2+y+z
x^2 + y + z

julia> homogeneous_components(f)
Dict{GrpAbFinGenElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  [2] => x^2 + y
  [3] => z

julia> R, x = polynomial_ring(QQ, "x" => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[1]^2+x[3]^2+x[5]^2
x[1]^2 + x[3]^2 + x[5]^2

julia> homogeneous_components(f)
Dict{GrpAbFinGenElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  [2 2 0 0] => x[5]^2
  [2 0 0 0] => x[1]^2 + x[3]^2
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
  dmat = vcat([d[i].coeff for i in 1:length(d)])
  tmat = zero_matrix(ZZ, 1, nvars(R))
  res_mat = zero_matrix(ZZ, 1, ncols(dmat))
  for (c, e) = Base.Iterators.zip(AbstractAlgebra.coefficients(forget_decoration(a)), AbstractAlgebra.exponent_vectors(forget_decoration(a)))
    # this is non-allocating
    for i in 1:length(e)
      tmat[1, i] = e[i]
    end
    mul!(res_mat, tmat, dmat)
    u = GrpAbFinGenElem(D, res_mat)
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
    homogeneous_component(f::MPolyDecRingElem, g::GrpAbFinGenElem)

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
julia> R, x = polynomial_ring(QQ, "x" => 1:5)
(Multivariate polynomial ring in 5 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[1]^2+x[3]^2+x[5]^2
x[1]^2 + x[3]^2 + x[5]^2

julia> homogeneous_component(f, 2*g[1])
x[1]^2 + x[3]^2

julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 0]
 [4, 1]

julia> R, x = graded_polynomial_ring(QQ, ["x[1]", "x[2]", "x[3]", "x[4]"], W)
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^2*x[2]+x[4]
x[1]^2*x[2] + x[4]

julia> homogeneous_component(f, [2, 1])
x[1]^2*x[2]

julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
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
function homogeneous_component(a::MPolyDecRingElem, g::GrpAbFinGenElem)
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
Nemo.ngens(W::MPolyDecRing) = Nemo.ngens(forget_decoration(W))
Nemo.gens(W::MPolyDecRing) = map(W, gens(forget_decoration(W)))
Nemo.gen(W::MPolyDecRing, i::Int) = W(gen(forget_decoration(W), i))
Base.getindex(W::MPolyDecRing, i::Int) = W(forget_decoration(W)[i])

base_ring(f::MPolyDecRingElem) = base_ring(forget_decoration(f))

function show_homo_comp(io::IO, M)
  (W, d) = get_attribute(M, :data)
  n = get_attribute(W, :name)
  if n != nothing
    print(io, "$(n)_$(d.coeff) of dim $(dim(M))")
  else
    println(io, "homogeneous component of $W of degree $d")
  end
end

@doc raw"""
    monomial_basis(R::MPolyDecRing, g::GrpAbFinGenElem)

Given a polynomial ring `R` over a field which is graded by a free
group of type `GrpAbFinGen`, and given an element `g` of that group,
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
julia> T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> G = grading_group(T)
GrpAb: Z

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
function monomial_basis(W::MPolyDecRing, d::GrpAbFinGenElem)
  #TODO: lazy: ie. no enumeration of points
  #      apparently it is possible to get the number of points faster than the points
  #TODO: in the presence of torsion, this is wrong. The component
  #      would be a module over the deg-0-sub ring.
  @req coefficient_ring(W) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  D = W.D
  is_free(D) || error("Grading group must be free")
  h = hom(free_abelian_group(ngens(W)), W.d)
  fl, p = haspreimage(h, d)
  R = base_ring(W)
  B = elem_type(W)[]
  if fl
     k, im = kernel(h)
     #need the positive elements in there...
     #Ax = b, Cx >= 0
     C = identity_matrix(FlintZZ, ngens(W))
     A = vcat([x.coeff for x = W.d])
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
    homogeneous_component(R::MPolyDecRing, g::GrpAbFinGenElem) 

Given a polynomial ring `R` over a field which is graded by a free
group of type `GrpAbFinGen`, and given an element `g` of that group,
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
    If the component is not finite dimensional, an error message will be thrown.

# Examples
```jldoctest
julia> R, x, y = polynomial_ring(QQ, "x" => 1:2, "y" => 1:3);

julia> W = [1 1 0 0 0; 0 0 1 1 1]
2×5 Matrix{Int64}:
 1  1  0  0  0
 0  0  1  1  1

julia> S, _ = grade(R, W);

julia> G = grading_group(S)
GrpAb: Z^2

julia> L = homogeneous_component(S, [1, 1]);

julia> L[1]
homogeneous component of Graded multivariate polynomial ring in 5 variables over QQ of degree [1 1]

julia> FG = gens(L[1]);

julia> EMB = L[2]
Map defined by a julia-function with inverse
  from s_[1 1] of dim 6
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
function homogeneous_component(W::MPolyDecRing, d::GrpAbFinGenElem)
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

function vector_space(K::AbstractAlgebra.Field, e::Vector{T}; target = nothing) where {T <:MPolyRingElem}
  local R
  if length(e) == 0
    R = target
    @assert R !== nothing
  else
    R = parent(e[1])
  end
  @assert base_ring(R) == K
  mon = Dict{elem_type(R), Int}()
  mon_idx = Vector{elem_type(R)}()
  M = sparse_matrix(K)
  last_pos = 1
  for i = e
    pos = Vector{Int}()
    val = Vector{elem_type(K)}()
    for (c, m) = Base.Iterators.zip(AbstractAlgebra.coefficients(i), AbstractAlgebra.monomials(i))
      if haskey(mon, m)
        push!(pos, mon[m])
        push!(val, c)
      else
        push!(mon_idx, m)
        mon[m] = last_pos
        push!(pos, last_pos)
        push!(val, c)
        last_pos += 1
      end
    end
    push!(M, sparse_row(K, pos, val))
  end
  rref!(M)

  b = Vector{elem_type(R)}()
  for i=1:nrows(M)
    s = zero(e[1])
    for (k,v) = M[i]
      s += v*mon_idx[k]
    end
    push!(b, s)
  end

  F = FreeModule(K, length(b); cached = false)
  function g(x::T)
    @assert parent(x) == R
    v = zero(F)
    for (c, m) = Base.Iterators.zip(AbstractAlgebra.coefficients(x), AbstractAlgebra.monomials(x))
      if !haskey(mon, m)
        error("not in image")
      end
      v += c*gen(F, mon[m])
    end
    return v
  end
  h = MapFromFunc(F, R, x -> sum([x[i] * b[i] for i in 1:length(b) if !is_zero_entry(x.v, 1, i)]; init = zero(R)), g)

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
  D = get_attribute!(() -> Dict{Any, Any}(), R, :relshp)::Dict{Any, Any}
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
  data::Vector{Int32}
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
    
    G = groebner_assure(I)
    h = Singular.hilbert_series(G.S, W)
    return new(h, W, I)
  end
  function HilbertData(B::IdealGens)
    return HilbertData(Oscar.MPolyIdeal(B))
  end
end

function hilbert_series(H::HilbertData)
  Zt, t = ZZ["t"]
  den = prod([1-t^w for w in H.weights])
  h = Zt(map(ZZRingElem, H.data[1:end-1]))
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
  Qt, t = QQ["t"]
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
  Qt, t = QQ["t"]
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
    expand(f::Frac{QQPolyRingElem}, d::Int) -> RelPowerSeries

Given a rational function $f$ over the rationals, expand $f$ as a power series
up to terms of degree $d$.

```jldoctest
julia> Qx, x = QQ["x"];

julia> expand(1//(1 - x^2), 5)
1 + t^2 + t^4 + O(t^6)
```
"""
function expand(f::Generic.Frac{QQPolyRingElem}, d::Int)
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
  print(io, "Hilbert Series for $(h.I), data: $(h.data), weights: $(h.weights)")  ###new
end

############################################################################
### Homogenization and Dehomogenization
############################################################################

## -----------------------------------
## 2023-07-14   Who wrote this code?
## (JAA)        Can it be deleted now?
## -----------------------------------

###old: special and to be applied with care
### needed for: AffinePlaneCurve.jl
### TODO: make adjustments there and omit function below
# function homogenization(f::MPolyRingElem, S::MPolyDecRing)
#     return homogenization(f, S, 1+ngens(parent(f)))
# end
function homogenization(f::MPolyRingElem, S::MPolyDecRing; pos::Union{Int,Nothing} = nothing)
  P = parent(f)
  if pos === nothing
    pos = 1+ngens(P)
  else
  @req  pos in 1:1+ngens(parent(f))  "Homog var index out of range."
  end
  d = total_degree(f)
  B = MPolyBuildCtx(S)
  for (c,e) = zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    insert!(e, pos, d-sum(e))
    push_term!(B, c, e)
  end
  return finish(B)
end

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

function _homogenization(f::MPolyRingElem, W::ZZMatrix, var::VarName, pos::Int)
    # ASSUME pos is in 1:1+ngens(parent(f))
  R = parent(f)
  A = copy(symbols(R))
  l = length(A)
  @req pos in 1:l+1 "Index out of range."
  if size(W, 1) == 1
     insert!(A, pos, Symbol(var))
  else
     for i = 1:size(W, 1)
       insert!(A, pos-1+i, _make_variable(var, i))
     end
  end
  L, _ = polynomial_ring(base_ring(R), A)
  Lgens = gens(L)
  ###ff = evaluate(f, vcat(Lgens[1:pos-1], Lgens[pos + size(W, 1):end]))
  G = abelian_group(zeros(Int, size(W, 1)))
  WH = hcat(Matrix(W[:, 1:(pos-1)]), Matrix(identity_matrix(ZZ, size(W, 1))[:, 1:size(W, 1)]), Matrix(W[:, pos:l]))
  S, _ = grade(L, [G(WH[:, i]) for i = 1:size(WH, 2)])
  ###return _homogenization(ff, S, pos)
  return _homogenization(f, S, pos)
end

function _homogenization(f::MPolyRingElem, W::Matrix{<:IntegerUnion}, var::VarName, pos::Int)
   W = matrix(ZZ, W)
   return _homogenization(f, W, var, pos)
end

@doc raw"""
    homogenization(f::MPolyRingElem, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName; pos::Int)
    homogenization(f::MPolyRingElem, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName)

If $m$ is the number of rows of `W`, extend the parent polynomial ring of `f` by inserting $m$ extra variables, starting at position `pos`
(if `pos` is not specified, it defaults to the position after the last variable).
Correspondingly, extend the integer matrix `W` by inserting the standard unit vectors of size $m$ as new columns, starting at column `pos`.
Grade the extended ring by converting the columns of the extended matrix to elements of the group $\mathbb Z^m$ and assigning these
as weights to the variables. Homogenize `f` with respect to the induced $\mathbb Z^m$-grading on the original ring, using the extra variables as homogenizing
variables. Return the result as an element of the extended ring with its $\mathbb Z^m$-grading. If $m=1$, the extra variable prints as `var`. Otherwise,
the extra variables print as `var[`$i$`]`, for $i = 1 \dots m$.

    homogenization(V::Vector{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName; pos::Int) where {T <: MPolyRingElem}
    homogenization(V::Vector{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName) where {T <: MPolyRingElem}

Given a vector `V` of elements in a common polynomial ring, create an extended ring with $\mathbb Z^m$-grading as above.
Homogenize the elements of `V` correspondingly,  and return the vector of homogenized elements.

    homogenization(I::MPolyIdeal{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName; pos::Int) where {T <: MPolyRingElem}
    homogenization(I::MPolyIdeal{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName) where {T <: MPolyRingElem}

Return the homogenization of `I` in an extended ring with $\mathbb Z^m$-grading as above.

!!! note
    If `W` comprises a single row of positive weights then the method used is essentially the same as
    for the standard-graded case: compute a wdegrevlex Groebner basis then homogenize its elements.
    Otherwise applied to an ideal `I`, the function first homogenizes the generators of `I` in the extended ring.
    It then creates the ideal generated by these homogenizations, and saturates this ideal
    with respect to the ideal which is generated by the product of the homogenizing variables.
    Source: Kreuzer+Robbiano (vol 2) Cor 4.3.8(a), Defn 4.4.1, Cor 4.4.9, Tutorial 37(g)

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> f = x^3+x^2*y+x*y^2+y^3
x^3 + x^2*y + x*y^2 + y^3

julia> W = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> F = homogenization(f, W, "z"; pos=3)
x^3*z[1]^3*z[2]^3 + x^2*y*z[1]^2*z[2]^2 + x*y^2*z[1]*z[2] + y^3

julia> parent(F)
Multivariate polynomial ring in 4 variables over QQ graded by
  x -> [1 3]
  y -> [2 4]
  z[1] -> [1 0]
  z[2] -> [0 1]
```
"""
function homogenization(f::MPolyRingElem, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName; pos::Union{Int,Nothing} = nothing)
  if pos === nothing
    pos = 1+ngens(parent(f))
  else
      @req  pos in 1:1+ngens(parent(f))  "Homog var index out of range."
  end
  return _homogenization(f, W, var, pos)
end

# # pos: default value is 1+ngens(parent(f))
# function homogenization(f::MPolyRingElem, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName)
#    return _homogenization(f, W, var, 1+ngens(parent(f)))
# end

# Without pos: default value determined by 1st elem of V
# function homogenization(V::Vector{T}, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName) where {T <: MPolyRingElem}
#     @assert !isempty(V)
#     # NB: call below checks that all elements of V are in the same ring
#     return homogenization(V, W, var, 1+ngens(parent(V[1])))
# end
function homogenization(V::Vector{T}, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName; pos::Union{Int,Nothing} = nothing) where {T <: MPolyRingElem}
    @assert !isempty(V) # OR ALLOW EMPTY V???   if isempty(V) return V end
  @assert all(x->parent(x) == parent(V[1]), V)
  R = parent(V[1])
  if pos === nothing
    pos = 1+ngens(R)
  end
  @req  pos in 1:ngens(R)+1  "Index out of range."
  A = copy(symbols(R))
  l = length(A)
  if size(W, 1) == 1
     insert!(A, pos, Symbol(var))
  else
     for i = 1:(size(W, 1))
       insert!(A, pos-1+i, _make_variable(var, i))
     end
  end
  G = abelian_group(zeros(Int, size(W, 1)))
  WH = hcat(Matrix(W[:, 1:(pos-1)]), Matrix(identity_matrix(ZZ, size(W, 1))[:, 1:size(W, 1)]), Matrix(W[:, pos:l]))
  S, _ = graded_polynomial_ring(base_ring(R), A, [G(WH[:, i]) for i = 1:size(WH, 2)])
  l = length(V)
  return [_homogenization(V[i], S, pos) for i=1:l]
end


## ------------------------------------------------------------------
## Homogenization of an ideal -- works for all gradings (incl. non-positive),
## but is slow.  If grading is positive over ZZ^1: use Singular (and wdegrevlex);
## if grading is positive over ZZ^m (with m > 1) use new "fast/clever" code.
## Author: John Abbott  2023-07-14
## ------------------------------------------------------------------

# To homogenize an ideal we look for some additional simple elements
# (as redundant generators), and then apply the standard algorithm.
# Seems wasteful, but in many cases leads to faster overall computation.


# Internal function: used only in _homogenization_via_saturation (immediately below)
# This function returns gens(I) and possibly some more "small" elements of I
# (obtained from 1 or more groebner bases of I).  If extra_gens_flag is true,
# we compute one more groebner bases of I in the hope that they contain some small elements.
function _gens_for_homog_via_sat(I::MPolyIdeal{T}, extra_gens_flag::Bool) where {T <: MPolyRingElem}
  ##    @req  !is_zero(I)  "Ideal must be non-zero"
  ##    @req  !is_one(I)   "Ideal must not be whole ring"
  if isempty(gens(I))  throw("Ideal must have at least 1 generator"); end;
  OrigR = parent(gens(I)[1])      # Since I is not zero, it has at least 1 gen.
  # Next few lines: we adjoin some more small, redundant gens.
  # !!HEURISTIC!!  We use 2*AveNumTerms as size limit for redundant gens we shall adjoin.
  AveNumTerms = sum([length(f)  for f in gens(I)])/length(gens(I)) ## floating-point!
  G = gens(I)
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
# Optional kwargs:
#   pos where to insert the homogenizing variables (default is after all other vars)
#   extra_gens_flag: if false, inhibits computing grobner basis speculatively in _gens_for_homog_via_sat
function _homogenization_via_saturation(I::MPolyIdeal{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName; pos::Union{Int,Nothing} = nothing, extra_gens_flag::Bool = false) where {T <: MPolyRingElem}
  P = base_ring(I)
  if pos === nothing
    pos = 1+ngens(P)
  else
    @req  pos in 1:1+ngens(P)  "Homog var index out of range."
  end
  if is_zero(I)  # special handling for ideal(0) incl. when there are no gens (?is there a cleaner way?)
    tmp = homogenization(zero(P), W, var; pos=pos)
    R = parent(tmp)
    return ideal(R) # no gens, so it is the zero ideal
  end # of special handling for zero ideal
  Hgens = homogenization(_gens_for_homog_via_sat(I, extra_gens_flag), W, var; pos=pos)
  R = parent(Hgens[1])
  if length(Hgens) == 1  # short-cut for principal ideal
    return ideal(R, Hgens)
  end
  DoSatByProduct = false  # true means saturate by product of homogenizing vars;
  #false means saturate successively by each homogenizing var (probably faster, in most cases?)
  if DoSatByProduct
    prod_h_vars = prod([gen(R,i)  for i in pos:(pos+size(W, 1)-1)])  # product of the homogenizing variables
    Ih = saturation(ideal(R, Hgens), ideal(R, prod_h_vars))
  else # not DoSatByProduct, so do a cascade of saturations
    Ih = ideal(R, Hgens)
    for i in pos:(pos+size(W, 1)-1)
      Ih = saturation(Ih, ideal(R, gen(R,i)))
    end
  end
  return Ih
  ######  2023-07-14: (this comment is now old, but probably still valid)
  ######  There is a problem/bug with Singular.satstd: so disable this version
  ######  IH = ideal(R, Hgens);  Y = ideal(R, prod_h_vars);
  ######  singular_assure(IH); singular_assure(Y); return Singular.satstd(IH.gens.S, Y.gens.S)
end


# ============================================
# 2023-06-30 START: New impl of homogenization
# By John Abbott, based on K+R Book "Computational Commutative Algebra" (see "justification" below)
# This code delegates to _homogenization_via_saturation for non-positive gradings
# ============================================

# sat_poly and kronecker_delta are internal auxiliary functions
# This homogenization impl seems to be usefully faster with positive ZZ^m-grading for m > 1

# Saturate a polynomial by a variable -- equiv to  f/gcd(f,h^deg(f))
function _sat_poly_by_var(f,h)
  # ASSUMES h is a variable i.e. monic poly of deg 1 with just 1 term
  # Is there a better way to implement this?
  Ih = ideal([h]);
  while ideal_membership(f, Ih)
    f = div(f,h); ##f /= h;  but this shorter form does not work (why?)
  end
  return f;
  # This loop is slower (presumably because it computes the remainder always)
  # while true
  #     q,r = divrem(f,h);
  #     if !is_zero(r)
  #         return q;
  #     end
  #     f = q;
  # end
end

# This function should be a file of utilities (JAA thinks)
function kronecker_delta(i::Int, j::Int)
  return (i == j) ? 1 : 0
end

# Check that W defines a positive grading: namely
#     each col contains a non-zero entry, and
#     first non-zero going down the col is positive.
# Follows defn 4.2.4 from R+K vol 2.
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


# This is the main exported function (for homogenizing ideals).  It has two optional kwargs
# (*) pos where to insert the homogenizing vars (default is right after the normal vars)
# (*) extra_gens_flag (relevant only for non-positive gradings); default is false, but if true
#     it triggers computation of a "needless" groebner basis in _gens_for_homog_via_sat which can speed up homogenization
function homogenization(I::MPolyIdeal{T}, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, h::VarName; pos::Union{Int,Nothing} = nothing, extra_gens_flag::Bool = false) where {T <: MPolyRingElem}
  ## ASSUME NumCols(W) = nvars(P)
  P = base_ring(I) # initial polynomial ring
  if pos === nothing
    pos = 1+ngens(P)
  else
    @req  pos in 1:1+ngens(P)  "Homog var index out of range."
  end
  # Handle zero ideal as special case: (quick and easy)
  if  is_zero(I)
    return ideal([homogenization(zero(P), W, h; pos=pos)])  # ideal(0) in correct ring
  end
  # # The fast method below is valid only for positive gradings; otherwise delegate to _homogenization_via_saturation
  if !is_positive_grading_matrix(W)
    return _homogenization_via_saturation(I, W, h; pos=pos, extra_gens_flag=extra_gens_flag)
  end
  # Henceforth we know that W defines a positive grading
  if size(W,1) == 1
    # Case: grading is over ZZ^1, so delegate to faster special-case function.
    # For correctness see Kreuzer+Robbiano (vol.2) Tutorial 53 part (d).
    if all(a -> (a == 1), W)
      W_ordering = degrevlex(P)
    else
      W_ordering = wdegrevlex(P, W[1,:])
    end
    GB = groebner_basis(I; ordering=W_ordering)
    return ideal(homogenization(elements(GB), W, h; pos=pos))
  end
  # Handle one ideal as special case (NB is_one might be costly!)
  if  is_one(I)
    return ideal([homogenization(one(P), W, h; pos=pos)])  # ideal(1) in correct ring
  end
  # Case: grading is over ZZ^k with k > 1 and is positive
  # >>> JUSTIFICATION <<<  from book Kreuzer+Robbiano "Computational Commutative Algebra"
  # Combination of Cor 4.3.8(a) [main algorithm],  Cor 4.4.9 [saturation by indet], and
  # Tutorial 37(g) [compute saturation by product via "cascade"]
  G = homogenization(gens(I), W, h; pos=1+ngens(P))  # h vars come last -- deliberately ignore pos here!
  Ph = parent(G[1])  # Ph is graded ring, "homog-extn" of P.
  if length(G) == 1  # short-cut for principal ideal
    return ideal(Ph, G);
  end
  N = ngens(Ph)
  num_x = ngens(P)
  num_h = ngens(Ph) - num_x
  # Build ordering matrix: weights matrix followed by identity mat, underneath is a revlex matrix
  Id = reduce(hcat, [[kronecker_delta(i,j)  for i in 1:num_h]  for j in 1:num_h])
  RevLexMat = reduce(hcat, [[-kronecker_delta(i+j, 1+ngens(Ph))  for i in 1:ngens(Ph)]  for j in 1:ngens(Ph)])
  M = hcat(W, Id)

  for j in 1:num_h
    M_complete = vcat(M, RevLexMat)
    JohnsOrdering = matrix_ordering(Ph, M_complete)
    Ih = ideal(G)
    G = groebner_basis(Ih; ordering=JohnsOrdering)  # option complete_reduction seemed to make no measurable difference
    G = elements(G)
    h_j = gen(Ph, ngens(Ph)+1-j)
    G = [_sat_poly_by_var(g, h_j)  for g in G]
    # Loop below: rotate cols in RevLexMat which corr to the "h" variables:
    for i in 1:num_h
      col1 = (i+j > num_h+1) ? (N+2+num_h-i-j) : (N+2-i-j)
      col2 = (col1 == N-num_h+1) ? N : col1-1
      RevLexMat[i, col1] = 0
      RevLexMat[i, col2] = -1
    end
  end
  if pos != 1+num_x
    # Caller wants the homogenizing variables in a non-standard place,
    # so permute the variables appropriately.  Is there a cleaner way to do this?
    junk = homogenization(gen(I,1), W, h; pos=pos)  # I just want the ring
    Ph2 = parent(junk)  # the desired graded ring
    images = [zero(Ph2)  for i in 1:N]
    for i in 1:pos-1
      images[i] = gen(Ph2,i)
    end
    for i in pos:num_x
      images[i] = gen(Ph2,i+num_h)
    end
    for i in 1:num_h
      images[num_x+i] = gen(Ph2, pos+i-1)
    end
    reorder_variables = hom(Ph, Ph2, images)
    G = [reorder_variables(g)  for g in G]
  end
  return ideal(G) # now in correct ring
end

# ============================================
# 2023-06-30 END: New impl of homogenization
# ============================================



@doc raw"""
    homogenization(f::MPolyRingElem, var::VarName)
    homogenization(f::MPolyRingElem, var::VarName; pos::Int)

    homogenization(V::Vector{T}, var::VarName) where {T <: MPolyRingElem}
    homogenization(V::Vector{T}, var::VarName; pos::Int) where {T <: MPolyRingElem}

    homogenization(I::MPolyIdeal{T}, var::VarName; ordering::Symbol = :degrevlex) where {T <: MPolyRingElem}
    homogenization(I::MPolyIdeal{T}, var::VarName; pos::Int, ordering::Symbol = :degrevlex) where {T <: MPolyRingElem}

Homogenize `f`, `V`, or `I` with respect to the standard $\mathbb Z$-grading using a homogenizing variable printing as `var`.
Return the result as an element of a graded polynomial ring with the homogenizing variable at position `pos`;
if `pos` is not specified it defaults to just after the last variable.

!!! note
    Applied to an ideal `I`, the function proceeds by homogenizing the elements of a Gröbner basis of `I` with respect to a degree compatible
    monomial ordering such as `degrevlex` (default).  If a Gröbner basis with respect to the specified ordering has not yet been computed and,
    thus, not yet been cached, executing the `homogenization` function with argument `I` may take some time.
    The degree compatibility of the specified ordering is not checked by the function.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> f = x^3-y^2-z
x^3 - y^2 - z

julia> F = homogenization(f, "w"; pos=4)
x^3 - y^2*w - z*w^2

julia> parent(F)
Multivariate polynomial ring in 4 variables over QQ graded by
  x -> [1]
  y -> [1]
  z -> [1]
  w -> [1]

julia> V = [y-x^2, z-x^3]
2-element Vector{QQMPolyRingElem}:
 -x^2 + y
 -x^3 + z

julia> homogenization(V, "w")
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 -x^2 + y*w
 -x^3 + z*w^2

julia> I = ideal(R, V)
ideal(-x^2 + y, -x^3 + z)

julia> PTC = homogenization(I, "w")
ideal(-x*z + y^2, x*y - z*w, x^2 - y*w)

julia> parent(PTC[1])
Multivariate polynomial ring in 4 variables over QQ graded by
  x -> [1]
  y -> [1]
  z -> [1]
  w -> [1]

julia> homogenization(I, "w"; ordering = deglex(gens(base_ring(I))))
ideal(x*z - y^2, x*y - z*w, x^2 - y*w, y^3 - z^2*w)
```
"""
function homogenization(f::MPolyRingElem, var::VarName; pos::Union{Int,Nothing} = nothing)
  R = parent(f)
  if pos === nothing
    pos = 1+ngens(R)
  else
    @req  pos in 1:ngens(R)+1  "Homog index out of range."
  end
  A = copy(symbols(R))
  l = length(A)
  insert!(A, pos, Symbol(var))
  S, _ = graded_polynomial_ring(base_ring(R), A)
  return _homogenization(f, S, pos)
end
# function homogenization(f::MPolyRingElem, var::VarName)
#     return homogenization(f, var, 1+ngens(parent(f)))
# end

# function homogenization(V::Vector{T}, var::VarName) where {T <: MPolyRingElem}
#     @assert !isempty(V)
#     return homogenization(V, var, 1+ngens(parent(V[1])))
# end
function homogenization(V::Vector{T}, var::VarName; pos::Union{Int,Nothing} = nothing) where {T <: MPolyRingElem}
  @req !isempty(V)  "homogenization of vector: vector must be non-empty"
  @assert all(x->parent(x) == parent(V[1]), V)
  R = parent(V[1])
  if pos === nothing
    pos = 1+ngens(R)
  else
    @req  pos in 1:ngens(R)+1  "Homog index out of range."
  end
  A = copy(symbols(R))
  l = length(A)
  insert!(A, pos, Symbol(var))
  S, _ = graded_polynomial_ring(base_ring(R), A)
  l = length(V)
  return [_homogenization(V[i], S, pos) for i=1:l]
end

function homogenization(I::MPolyIdeal{T}, var::VarName; pos::Union{Int,Nothing} = nothing, ordering::MonomialOrdering = default_ordering(base_ring(I))) where {T <: MPolyRingElem}
  R = base_ring(I)
  if pos === nothing
    pos = 1+ngens(R)
  else
    @req  pos in 1:ngens(R)+1  "Homog index out of range."
  end
  # Handle zero ideal as special case: (quick and easy)
  if  is_zero(I)
    return ideal([homogenization(zero(R), var; pos=pos)])  # zero ideal in correct ring
  end
  # TODO: Adjust as soon as new GB concept is implemented   [[@wdecker: delete this comment?]]
  return ideal(homogenization(gens(groebner_basis(I; ordering)), var; pos=pos))
end
# function homogenization(I::MPolyIdeal{T}, var::VarName,; ordering::MonomialOrdering = default_ordering(base_ring(I))) where {T <: MPolyRingElem}
#   # TODO: Adjust as soon as new GB concept is implemented
#   return ideal(homogenization(gens(groebner_basis(I; ordering)), var, 1+ngens(base_ring(I))))
# end

### needed for: PlaneCurve-test.jl, ProjPlaneCurve.jl
### TODO: make adjustments there and omit function below
function dehomogenization(F::MPolyDecRingElem, R::MPolyRing, pos::Int)
  B = MPolyBuildCtx(R)
  for (c,e) = zip(AbstractAlgebra.coefficients(F), AbstractAlgebra.exponent_vectors(F))
    deleteat!(e, pos)
    push_term!(B, c, e)
  end
  return finish(B)
end



function _dehomogenization(F::MPolyDecRingElem, R::MPolyRing, pos::Int, m::Int)
  B = MPolyBuildCtx(R)
  for (c,e) = zip(AbstractAlgebra.coefficients(F), AbstractAlgebra.exponent_vectors(F))
    deleteat!(e, pos:pos+m-1)
    push_term!(B, c, e)
  end
  return finish(B)
end


@doc raw"""
    dehomogenization(F::MPolyDecRingElem, pos::Int)

Given an element `F` of a $\mathbb Z^m$-graded ring, where the generators of $\mathbb Z^m$
are the assigned weights to the variables at positions `pos`, $\dots$, `pos` $-1+m$, dehomogenize
`F` using the variables at these positions. Return the result as an element of a polynomial ring not
depending on the variables at these positions.

    dehomogenization(V::Vector{T}, pos::Int) where {T <: MPolyDecRingElem}

Given a vector `V` of elements in a common graded polynomial ring, create a polynomial ring not
depending on the variables at positions `pos`, $\dots$, `pos` $-1+m$. Dehomogenize the elements 
of `V` correspondingly,  and return the vector of dehomogenized elements.  

    dehomogenization(I::MPolyIdeal{T}, pos::Int) where {T <: MPolyDecRingElem}

Return the dehomogenization of `I` in a polynomial ring as above.

# Examples
```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> F = x^3-x^2*y-x*z^2
x^3 - x^2*y - x*z^2

julia> f = dehomogenization(F, 1)
-y - z^2 + 1

julia> parent(f)
Multivariate polynomial ring in 2 variables y, z
  over rational field

julia> V = [x*y-z^2, x^2*z-x^3]
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x*y - z^2
 -x^3 + x^2*z

julia> dehomogenization(V, 3)
2-element Vector{QQMPolyRingElem}:
 x*y - 1
 -x^3 + x^2

julia> I = ideal(S, V)
ideal(x*y - z^2, -x^3 + x^2*z)

julia> dehomogenization(I, 3)
ideal(x*y - 1, -x^3 + x^2)

julia> W = [1 2 1 0; 3 4 0 1]
2×4 Matrix{Int64}:
 1  2  1  0
 3  4  0  1

julia> S, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], W)
(Graded multivariate polynomial ring in 4 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w, x, y, z])

julia> F = w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3
w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3

julia> dehomogenization(F, 3)
w^3 + w^2*x + w*x^2 + x^3
```
"""
function dehomogenization(F::MPolyDecRingElem, pos::Int)
  S = parent(F)
  @assert is_zm_graded(S)
  l = ngens(S)
  G = grading_group(S)
  m = ngens(G)
  @req pos in 1:l-m+1 "Index out of range."
  for i = 1:m
      @assert degree(gen(S, pos-1+i)) == gen(G, i)
  end
  A = copy(symbols(S))
  deleteat!(A, pos:pos+m-1)
  R, _ = polynomial_ring(base_ring(S), A)
  return _dehomogenization(F, R, pos, m)
end
function dehomogenization(V::Vector{T}, pos::Int) where {T <: MPolyDecRingElem}
  @assert all(x->parent(x) == parent(V[1]), V)
  S = parent(V[1])
  @assert is_zm_graded(S)
  l = ngens(S)
  G = grading_group(S)
  m = ngens(G)
  @req pos in 1:l-m+1 "Index out of range."
  for i = 1:m
      @assert degree(gen(S, pos-1+i)) == gen(G, i)
  end
  A = copy(symbols(S))
  deleteat!(A, pos:pos+m-1)
  R, _ = polynomial_ring(base_ring(S), A)
  l = length(V)
  return [_dehomogenization(V[i], R, pos, m) for i=1:l]
end
function dehomogenization(I::MPolyIdeal{T}, pos::Int) where {T <: MPolyDecRingElem}
  return ideal(dehomogenization(gens(I), pos))
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

Return the ungraded undecorated ring.
"""
forget_grading(R::MPolyDecRing) = forget_decoration(R)

@doc raw"""
    forget_decoration(f::MPolyDecRingElem)

Return the element in the underlying undecorated ring.
"""
forget_decoration(f::MPolyDecRingElem) = f.f

@doc raw"""
    forget_grading(f::MPolyDecRingElem)

Return the element in the underlying ungraded ring.
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

Return the ideal in the underlying ungraded ring.
"""
forget_grading(I::MPolyIdeal{<:MPolyDecRingElem}) = forget_decoration(I)

### This is a temporary fix that needs to be addressed in AbstractAlgebra, issue #1105.
# TODO: This still seems to be not resolved!!!
Generic.ordering(S::MPolyDecRing) = :degrevlex


#############truncation#############

@doc raw"""
    truncate(I::MPolyIdeal, g::GrpAbFinGenElem)

Given a (homogeneous) ideal `I` in a $\mathbb Z$-graded multivariate polynomial ring
with positive weights, return the truncation of `I` at degree `g`.

    truncate(I::MPolyIdeal, d::Int)

Given an ideal `I` as above, and given an integer `d`, convert `d`
into an element `g` of the grading group of `base_ring(I)` and proceed as above.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = ideal(R, [x, y^4, z^6])
ideal(x, y^4, z^6)

julia> truncate(I, 3)
ideal(x*z^2, x*y*z, x*y^2, x^2*z, x^2*y, x^3, y^4, z^6)
```

```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [3,2,1]);

julia> I = ideal(R, [x, y^4, z^6])
ideal(x, y^4, z^6)

julia> truncate(I, 3)
ideal(x, y^4, z^6)

julia> truncate(I, 4)
ideal(x*z, z^6, y^4)
```
"""
function truncate(I::MPolyIdeal, g::GrpAbFinGenElem)
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
  V = sort(gens(I), lt = (a, b) -> degree(Int, a) <= degree(Int, b))
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
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> V = [x, z^2, x^3+y^3, y^4, y*z^5];

julia> I = ideal(R, V)
ideal(x, z^2, x^3 + y^3, y^4, y*z^5)

julia> minimal_generating_set(I)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x
 z^2
 y^3

julia> I = ideal(R, zero(R))
ideal(0)

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
    G.gens.S.isGB = true
    _, sing_min = Singular.mstd(singular_generators(G, G.ord))
    return filter(!iszero, (R).(gens(sing_min)))
  else
    sing_gb, sing_min = Singular.mstd(singular_generators(I))
    ring = I.gens.Ox
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
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [y^2*z − x^2*w, z^4 − x*w^3]);

julia> cm_regularity(I)
6

julia> minimal_betti_table(I);
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
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> I = ideal(R, [y-x^2, x-z^3])
ideal(-x^2 + y, x - z^3)

julia> degree(I)
6
```
"""
@attr Int function degree(I::MPolyIdeal)
  is_standard_graded(base_ring(I)) || (I = homogenization(I, "_h"))
  return Int(degree(HilbertData(I)))
end
