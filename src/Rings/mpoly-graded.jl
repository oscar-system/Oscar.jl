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

If `R` is, say, `G`-graded, return `G`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

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


function show(io::IO, W::MPolyDecRing)
  Hecke.@show_name(io, W)
  Hecke.@show_special(io, W)
  if is_filtered(W)
    println(io, "$(forget_decoration(W)) filtrated by ")
  else
    println(io, "$(forget_decoration(W)) graded by ")
  end
  R = forget_decoration(W)
  g = gens(R)
  for i = 1:ngens(R)
    if i == ngens(R)
       print(io, "  $(g[i]) -> $(W.d[i].coeff)")
    else  
       println(io, "  $(g[i]) -> $(W.d[i].coeff)")
    end  
  end
#  println(IOContext(io, :compact => true, ), W.d)
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
(Multivariate Polynomial Ring in x, y, z over Rational Field, QQMPolyRingElem[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> T, (x, y, z) = grade(R)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])
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
(Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field, QQMPolyRingElem[x[1], x[2]], QQMPolyRingElem[y[1], y[2], y[3]])

julia> W = [1 1 0 0 0; 0 0 1 1 1]
2×5 Matrix{Int64}:
 1  1  0  0  0
 0  0  1  1  1

julia> grade(R, W)
(Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  y[1] -> [0 1]
  y[2] -> [0 1]
  y[3] -> [0 1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], y[1], y[2], y[3]])
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
(Multivariate Polynomial Ring in x, y, z over Rational Field, QQMPolyRingElem[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

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
(Multivariate Polynomial Ring in x, y, z over Rational Field, QQMPolyRingElem[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

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
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field graded by 
  x[1] -> [1 0 1 1]
  x[2] -> [0 1 0 1]
  x[3] -> [1 0 1 0]
  x[4] -> [0 1 0 0]
  x[5] -> [1 1 0 0], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> is_zm_graded(S)
false

julia> G = abelian_group(ZZMatrix([1 -1]));

julia> g = gen(G, 1)
Element of
(General) abelian group with relation matrix
[1 -1]
with components [0 1]

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
(Multivariate Polynomial Ring in t, x, y over Rational Field, QQMPolyRingElem[t, x, y])

julia> G = abelian_group([0])
GrpAb: Z

julia> S, (t, x, y) = grade(R, [-1, 1, 1])
(Multivariate Polynomial Ring in t, x, y over Rational Field graded by
  t -> [-1]
  x -> [1]
  y -> [1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[t, x, y])

julia> is_positively_graded(S)
false

julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, QQMPolyRingElem[x, y])

julia> G = abelian_group([0, 2])
(General) abelian group with relation matrix
[0 0; 0 2]

julia> W = [gen(G, 1)+gen(G, 2), gen(G, 1)]
2-element Vector{GrpAbFinGenElem}:
 Element of
(General) abelian group with relation matrix
[0 0; 0 2]
with components [1 1]
 Element of
(General) abelian group with relation matrix
[0 0; 0 2]
with components [1 0]

julia> S, (x, y) = grade(R, W)
(Multivariate Polynomial Ring in x, y over Rational Field graded by 
  x -> [1 1]
  y -> [1 0], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> is_positively_graded(S)
false
```
"""
function is_positively_graded(R::MPolyDecRing)
  is_graded(R) || return false
  G = grading_group(R)
  is_free(G) || return false
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
(Multivariate Polynomial Ring in x1, x2, x3, x4 over Rational Field graded by
  x1 -> [1 0]
  x2 -> [0 1]
  x3 -> [1 0]
  x4 -> [4 1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x1, x2, x3, x4])

julia> S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])
```
"""
function graded_polynomial_ring(C::Ring, V::Union{Tuple{Vararg{T}}, AbstractVector{T}},
      W=ones(Int, length(V)); ordering=:lex) where
      T<:VarName
   return grade(polynomial_ring(C, V; ordering)[1], W)
end

function graded_polynomial_ring(C::Ring, n::Int, s::VarName,
      W=ones(Int, n); ordering=:lex)
   return grade(polynomial_ring(C, n, s; ordering)[1], W)
end

filtrate(R::MPolyRing) = decorate(R)

function show_special_elem_grad(io::IO, a::GrpAbFinGenElem)
  if get(io, :compact, false)
    print(io, a.coeff)
  else
    print(io, "graded by $(a.coeff)")
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
(Multivariate Polynomial Ring in t, x, y over Rational Field, QQMPolyRingElem[t, x, y])

julia> typeof(R)
QQMPolyRing

julia>  typeof(x)
QQMPolyRingElem

julia> G = abelian_group([0])
GrpAb: Z

julia> g = gen(G, 1)
Element of
GrpAb: Z
with components [1]

julia> S, (t, x, y) = grade(R, [-g, g, g])
(Multivariate Polynomial Ring in t, x, y over Rational Field graded by
  t -> [-1]
  x -> [1]
  y -> [1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[t, x, y])

julia> typeof(S)
MPolyDecRing{QQFieldElem, QQMPolyRing}

julia> S isa MPolyRing
true

julia> typeof(x)
MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}

julia> R, x, y = polynomial_ring(QQ, "x" => 1:2, "y" => 1:3)
(Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field, QQMPolyRingElem[x[1], x[2]], QQMPolyRingElem[y[1], y[2], y[3]])

julia> G = abelian_group([0, 0])
GrpAb: Z^2

julia> g = gens(G)
2-element Vector{GrpAbFinGenElem}:
 Element of
GrpAb: Z^2
with components [1 0]
 Element of
GrpAb: Z^2
with components [0 1]

julia> W = [g[1], g[1], g[2], g[2], g[2]];

julia> S, _ = grade(R, W)
(Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  y[1] -> [0 1]
  y[2] -> [0 1]
  y[3] -> [0 1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], y[1], y[2], y[3]])

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
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G)
4-element Vector{GrpAbFinGenElem}:
 Element of
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]
with components [1 0 0 0]
 Element of
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]
with components [0 1 0 0]
 Element of
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]
with components [0 0 1 0]
 Element of
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]
with components [0 0 0 1]

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field graded by
  x[1] -> [1 0 1 1]
  x[2] -> [0 1 0 1]
  x[3] -> [1 0 1 0]
  x[4] -> [0 1 0 0]
  x[5] -> [1 1 0 0], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])
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
function (W::MPolyDecRing{S, T})(f::U) where {S, T, U}
  if parent_type(U) === T
    @assert forget_decoration(W) === parent(f)
    return MPolyDecRingElem(f, W)
  else
    return W(forget_decoration(W)(f))
  end
end

function (W::MPolyDecRing{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T
  return W(forget_decoration(W)(c, e))
end

(W::MPolyDecRing)(g::MPolyDecRingElem) = MPolyDecRingElem(forget_decoration(g), W)
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
  fl, w_ord = has_weighted_ordering(R)
  if fl
    return w_ord
  end
  return degrevlex(gens(R))
end

function singular_poly_ring(R::MPolyDecRing; keep_ordering::Bool = false)
  if !keep_ordering
    return singular_poly_ring(forget_decoration(R), default_ordering(R))
  end
  return singular_poly_ring(forget_decoration(R), keep_ordering = keep_ordering)
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
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field graded by
  x[1] -> [1 0 1 1]
  x[2] -> [0 1 0 1]
  x[3] -> [1 0 1 0]
  x[4] -> [0 1 0 0]
  x[5] -> [1 1 0 0], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[2]^2+2*x[4]^2
x[2]^2 + 2*x[4]^2

julia> degree(f)
Element of
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]
with components [0 2 0 0]

julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 0]
 [4, 1]

julia> R, x = graded_polynomial_ring(QQ, ["x[1]", "x[2]", "x[3]", "x[4]"], W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [0 1]
  x[3] -> [1 0]
  x[4] -> [4 1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^4*x[2]+x[4]
x[1]^4*x[2] + x[4]

julia> degree(f)
graded by [4 1]

julia> degree(Vector{Int}, f)
2-element Vector{Int64}:
 4
 1

julia>  R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^6+y^3+z^2
x^6 + y^3 + z^2

julia> degree(f)
graded by [6]

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
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^2+y*z
x^2 + y*z

julia> is_homogeneous(f)
false

julia> W = [1 2 1 0; 3 4 0 1]
2×4 Matrix{Int64}:
 1  2  1  0
 3  4  0  1

julia> S, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"], W)
(Multivariate Polynomial Ring in w, x, y, z over Rational Field graded by
  w -> [1 3]
  x -> [2 4]
  y -> [1 0]
  z -> [0 1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w, x, y, z])

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
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> f = x^2+y+z
x^2 + y + z

julia> homogeneous_components(f)
Dict{GrpAbFinGenElem, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} with 2 entries:
  [2] => x^2 + y
  [3] => z

julia> R, x = polynomial_ring(QQ, "x" => 1:5)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field graded by 
  x[1] -> [1 0 1 1]
  x[2] -> [0 1 0 1]
  x[3] -> [1 0 1 0]
  x[4] -> [0 1 0 0]
  x[5] -> [1 1 0 0], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

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
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, QQMPolyRingElem[x[1], x[2], x[3], x[4], x[5]])

julia> G = abelian_group([0, 0, 2, 2])
(General) abelian group with relation matrix
[0 0 0 0; 0 0 0 0; 0 0 2 0; 0 0 0 2]

julia> g = gens(G);

julia> W = [g[1]+g[3]+g[4], g[2]+g[4], g[1]+g[3], g[2], g[1]+g[2]];

julia> S, x = grade(R, W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field graded by 
  x[1] -> [1 0 1 1]
  x[2] -> [0 1 0 1]
  x[3] -> [1 0 1 0]
  x[4] -> [0 1 0 0]
  x[5] -> [1 1 0 0], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4], x[5]])

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
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4] over Rational Field graded by 
  x[1] -> [1 0]
  x[2] -> [0 1]
  x[3] -> [1 0]
  x[4] -> [4 1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^2*x[2]+x[4]
x[1]^2*x[2] + x[4]

julia> homogeneous_component(f, [2, 1])
x[1]^2*x[2]

julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

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
    homogeneous_component(R::MPolyDecRing, g::GrpAbFinGenElem) 

Given a polynomial ring `R` over a field which is graded by a free
group of type `GrpAbFinGen`, and given an element `g` of that group,
return the homogeneous component of `R` of degree `g`. Additionally, return
the embedding of the component into `R`.

    homogeneous_component(R::MPolyDecRing, g::Vector{<:IntegerUnion})

Given a $\mathbb  Z^m$-graded polynomial ring `R` over a field, and given
a vector `g` of $m$ integers, convert `g` into an element of the grading
group of `R`, and return the homogeneous component of `R` whose degree 
is that element. Additionally, return the embedding of the component into `R`.

    homogeneous_component(R::MPolyDecRing, g::IntegerUnion)

Given a $\mathbb  Z$-graded polynomial ring `R` over a field, and given
an integer `g`, convert `g` into an element of the grading group of `R`, 
and return the homogeneous component of `R` whose degree is that element.
Additionally, return the embedding of the component into `R`.

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
homogeneous component of Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  y[1] -> [0 1]
  y[2] -> [0 1]
  y[3] -> [0 1] of degree graded by [1 1]

julia> FG = gens(L[1]);

julia> EMB = L[2]
Map from
homogeneous component of Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  y[1] -> [0 1]
  y[2] -> [0 1]
  y[3] -> [0 1] of degree graded by [1 1]
 to Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [1 0]
  y[1] -> [0 1]
  y[2] -> [0 1]
  y[3] -> [0 1] defined by a julia-function with inverse

julia> for i in 1:length(FG) println(EMB(FG[i])) end
x[2]*y[3]
x[2]*y[2]
x[2]*y[1]
x[1]*y[3]
x[1]*y[2]
x[1]*y[1]

julia> T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> G = grading_group(T)
GrpAb: Z

julia> L = homogeneous_component(T, 2)
(homogeneous component of Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] of degree graded by [2]
, Map from
homogeneous component of Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] of degree graded by [2]
 to Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1] defined by a julia-function with inverse)

julia> FG = gens(L[1]);

julia> EMB = L[2];

julia> for i in 1:length(FG) println(EMB(FG[i])) end
z^2
y*z
y^2
x*z
x*y
x^2
```
"""
function homogeneous_component(W::MPolyDecRing, d::GrpAbFinGenElem)
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

  F = FreeModule(K, length(b), cached = false)
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
  h = MapFromFunc(x -> sum(x[i] * b[i] for i in 1:length(b); init = zero(R)), g, F, R)

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
    @req is_graded(R) "The base ring must be graded"

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

function hilbert_series(H::HilbertData, i::Int= 1)
  Zt, t = ZZ["t"]
  den = prod([1-t^w for w in H.weights])
  if i==1   ### the Hilbert series with denominator prod (1-t^H.weights[i])
    return Zt(map(ZZRingElem, H.data[1:end-1])), den
  elseif i==2   ### the reduced Hilbert series
    h = hilbert_series(H, 1)[1]
    c = gcd(h, den)
    h = divexact(h, c)
    den = divexact(den, c)
    return den(0)*h, den(0)*den
  end
  error("2nd parameter must be 1 or 2")
end

#Decker-Lossen, p23/24
function hilbert_polynomial(H::HilbertData)

  @req all(isone, H.weights) "All weights must be 1"
  
  q, dn = hilbert_series(H, 2)
  a = QQFieldElem[]
  nf = QQFieldElem(1)
  d = degree(dn)-1
  for i=1:d+1
    push!(a, q(1)//nf)    
    nf *= i
    q = derivative(q)
  end
  Qt, t = QQ["t"]
  t = gen(Qt)
  bin = one(parent(t))
  b = QQPolyRingElem[]
  if d==-1 return zero(parent(t)) end
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
     q, _ = hilbert_series(H, 2)
     return q(1)
  end
  deg = leading_coefficient(P)*factorial(ZZ(degree(P)))
  @assert isone(denominator(deg))
  return numerator(deg)
end

function (P::QQRelPowerSeriesRing)(H::HilbertData)
  n, d = hilbert_series(H, 2)
  Qt, t = QQ["t"]
  nn = map_coefficients(QQ, n, parent = Qt)
  dd = map_coefficients(QQ, d, parent = Qt)
  gg, ee, _ = gcdx(dd, gen(Qt)^max_precision(P))
  @assert isone(gg)
  nn = Hecke.mullow(nn, ee, max_precision(P))
  c = collect(coefficients(nn))
  return P(map(QQFieldElem, c), length(c), max_precision(P), 0)
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

########################################################################
# Compute the Hilbert series from the monomial leading ideal using 
# different bisection strategies along the lines of 
# [Kreuzer, Robbiano: Computational Commutative Algebra 2, Springer]
# Section 5.3.
########################################################################

_divides(a::Vector{Int}, b::Vector{Int}) = all(k->(a[k]>=b[k]), 1:length(a))

function _gcd(a::Vector{Int}, b::Vector{Int})
  return [a[k] > b[k] ? b[k] : a[k] for k in 1:length(a)]
end

function _are_pairwise_coprime(a::Vector{Vector{Int}})
  length(a) <= 1 && return true
  n = length(a)
  m = length(a[1])
  for i in 1:n-1
    for j in i+1:n
      all(k -> a[i][k] == 0 || a[j][k] == 0, 1:m) || return false
    end
  end
  return true
end

### Assume that a is minimal. We divide by `pivot` and return 
# a minimal set of generators for the resulting monomial ideal.
function _divide_by(a::Vector{Vector{Int}}, pivot::Vector{Int})
  # The good ones will contribute to a minimal generating 
  # set of the lhs ideal.
  #
  # The bad monomials come from those which hop over the boundaries of 
  # the monomial diagram by the shift. Their span has a new 
  # generator which is collected in `bad` a priori. It is checked 
  # whether they become superfluous and if not, they are added to 
  # the good ones.
  good = sizehint!(Vector{Vector{Int}}(), length(a))
  bad = sizehint!(Vector{Vector{Int}}(), length(a))
  for e in a
    if _divides(e, pivot)
      push!(good, e - pivot)
    else
      push!(bad, e)
    end
  end

  # pre-allocate m so that don't need to allocate it again each loop iteration
  m = similar(pivot)
  for e in bad
    # the next line computers   m = [k < 0 ? 0 : k for k in e]
    # but without allocations
    for i in 1:length(e)
      m[i] = max(e[i] - pivot[i], 0)
    end

    # check whether the new monomial m is already in the span 
    # of the good ones. If yes, discard it. If not, discard those 
    # elements of the good ones that are in the span of m and put 
    # m in the list of good ones instead. 
    if all(x->!_divides(m, x), good)
      # Remove those 'good' elements which are multiples of m
      filter!(x -> !_divides(x, m), good)
      push!(good, copy(m))
    end
  end

  return good
end

### This implements the speedup from the CoCoA strategy, 
# see Exercise 4 in Section 5.3
function _divide_by_monomial_power(a::Vector{Vector{Int}}, j::Int, k::Int)
  divides(a::Vector{Int}, b::Vector{Int}) = all(k->(k>=0), b[1:j-1] - a[1:j-1]) && all(k->(k>=0), b[j+1:end] - a[j+1:end])
  kept = [e for e in a if e[j] == 0]
  for i in 1:k
    next_slice = [e for e in a if e[j] == i]
    still_kept = next_slice
    for e in kept 
      all(v->!divides(v, e), next_slice) && push!(still_kept, e)
    end
    kept = still_kept
  end
  still_kept = vcat(still_kept, [e for e in a if e[j] > k])

  result = Vector{Vector{Int}}()
  for e in still_kept
    v = copy(e)
    v[j] = (e[j] > k ? e[j] - k : 0)
    push!(result, v)
  end
  return result
end

function _find_maximum(a::Vector{Int})
  m = a[1]
  j = 1
  for i in 1:length(a)
    if a[i] > m 
      m = a[i]
      j = i
    end
  end
  return j
end

#######################################################################
# 06.10.2022
# The following internal routine can probably still be tuned.
#
# Compared to the implementation in Singular, it is still about 10 
# times slower. We spotted the following possible deficits:
#
#   1) The returned polynomials are not yet built with MPolyBuildCtx.
#      We tried that, but as for now, the code for the build context 
#      is even slower than the direct implementation that is in place 
#      now. Once MPolyBuildCtx is tuned, we should try that again; 
#      the code snippets are still there. In total, building the return 
#      values takes up more than one third of the computation time 
#      at this moment.
#
#   2) Numerous allocations for integer vectors. The code below 
#      performs lots of iterative allocations for lists of integer 
#      vectors. If we constructed a container data structure to 
#      maintain this list internally and do the allocations at once  
#      (for instance in a big matrix), this could significantly 
#      speed up the code. However, that is too much work for 
#      the time being, as long as a high-performing Hilbert series 
#      computation is not of technical importance.
#
#   3) Singular uses bitmasking for exponent vectors to decide on 
#      divisibility more quickly. This is particularly important in 
#      the method _divide_by. For singular, this leads to a speedup 
#      of factor 5. Here, only less than one third of the time is 
#      spent in `_divide_by`, so it does not seem overly important 
#      at this point, but might become relevant in the future. 
#      A particular modification of the singular version of bitmasking 
#      is to compute means for the exponents occurring for each variable 
#      and set the bits depending on whether a given exponent is greater 
#      or less than that mean value. 

function _hilbert_numerator_from_leading_exponents(
    a::Vector{Vector{Int}},
    weight_matrix::Matrix{Int},
    return_ring::Ring,
    #algorithm=:generator
    #algorithm=:custom
    #algorithm=:gcd
    #algorithm=:indeterminate
    #algorithm=:cocoa
    algorithm::Symbol # =:BayerStillmanA, # This is by far the fastest strategy. Should be used.
    # short exponent vectors where the k-th bit indicates that the k-th 
    # exponent is non-zero.
  )
  t = gens(return_ring)
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, return_ring, t)
  ret !== nothing && return ret

  if algorithm == :BayerStillmanA
    return _hilbert_numerator_bayer_stillman(a, weight_matrix, return_ring, t)
  elseif algorithm == :custom
    return _hilbert_numerator_custom(a, weight_matrix, return_ring, t)
  elseif algorithm == :gcd # see Remark 5.3.11
    return _hilbert_numerator_gcd(a, weight_matrix, return_ring, t)
  elseif algorithm == :generator # just choosing on random generator, cf. Remark 5.3.8
    return _hilbert_numerator_generator(a, weight_matrix, return_ring, t)
  elseif algorithm == :indeterminate # see Remark 5.3.8
    return _hilbert_numerator_indeterminate(a, weight_matrix, return_ring, t)
  elseif algorithm == :cocoa # see Remark 5.3.14
    return _hilbert_numerator_cocoa(a, weight_matrix, return_ring, t)
  end
  error("invalid algorithm")
end

# compute t ^ (weight_matrix * expvec), where t == gens(S)
function _expvec_to_poly(S::Ring, t::Vector, weight_matrix::Matrix{Int}, expvec::Vector{Int})
  o = one(coefficient_ring(S))  # TODO: cache this ?!
  return S([o], [weight_matrix * expvec])
end

# special case for univariate polynomial ring
function _expvec_to_poly(S::PolyRing, t::Vector, weight_matrix::Matrix{Int}, expvec::Vector{Int})
  @assert length(t) == 1
  @assert size(weight_matrix) == (1, length(expvec))

  # compute the dot-product of weight_matrix[1,:] and expvec, but faster than dot
  # TODO: what about overflows?
  s = weight_matrix[1,1] * expvec[1]
  for i in 2:length(expvec)
    @inbounds s += weight_matrix[1,i] * expvec[i]
  end
  return t[1]^s
end

function _hilbert_numerator_trivial_cases(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector, oneS = one(S)
  )
  length(a) == 0 && return oneS

  # See Proposition 5.3.6
  if _are_pairwise_coprime(a)
    return prod(oneS - _expvec_to_poly(S, t, weight_matrix, e) for e in a)
  end

  return nothing
end


function _hilbert_numerator_cocoa(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  n = length(a)
  m = length(a[1])
  counters = [0 for i in 1:m]
  for e in a
    counters += [iszero(k) ? 0 : 1 for k in e]
  end
  j = _find_maximum(counters)

  p = a[rand(1:n)]
  while p[j] == 0
    p = a[rand(1:n)]
  end

  q = a[rand(1:n)]
  while q[j] == 0 || p == q
    q = a[rand(1:n)]
  end

  pivot = [0 for i in 1:m]
  pivot[j] = minimum([p[j], q[j]])


  ### Assembly of the quotient ideal with less generators
  rhs = [e for e in a if !_divides(e, pivot)]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by_monomial_power(a, j, pivot[j])

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_cocoa(rhs, weight_matrix, S, t) + f*_hilbert_numerator_cocoa(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_indeterminate(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  e = first(a)
  found_at = findfirst(!iszero, e)::Int64
  pivot = zero(e)
  pivot[found_at] = 1

  ### Assembly of the quotient ideal with less generators
  rhs = [e for e in a if e[found_at] == 0]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by(a, pivot)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_indeterminate(rhs, weight_matrix, S, t) + f*_hilbert_numerator_indeterminate(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_generator(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  b = copy(a)
  pivot = pop!(b)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  c = _divide_by(b, pivot)
  p1 = _hilbert_numerator_generator(b, weight_matrix, S, t)
  p2 = _hilbert_numerator_generator(c, weight_matrix, S, t)

  return p1 - f * p2
end

function _hilbert_numerator_gcd(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  n = length(a)
  counters = [0 for i in 1:length(a[1])]
  for e in a
    counters += [iszero(k) ? 0 : 1 for k in e]
  end
  j = _find_maximum(counters)

  p = a[rand(1:n)]
  while p[j] == 0
    p = a[rand(1:n)]
  end

  q = a[rand(1:n)]
  while q[j] == 0 || p == q
    q = a[rand(1:n)]
  end

  pivot = _gcd(p, q)

  ### Assembly of the quotient ideal with less generators
  rhs = [e for e in a if !_divides(e, pivot)]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by(a, pivot)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_gcd(rhs, weight_matrix, S, t) + f*_hilbert_numerator_gcd(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_custom(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector
  )
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t)
  ret !== nothing && return ret

  p = Vector{Int}()
  q = Vector{Int}()
  max_deg = 0
  for i in 1:length(a)
    b = a[i]
    for j in i+1:length(a)
      c = a[j]
      r = _gcd(b, c)
      if sum(r) > max_deg
        max_deg = sum(r)
        p = b
        q = c
      end
    end
  end

  ### Assembly of the quotient ideal with less generators
  pivot = _gcd(p, q)
  rhs = [e for e in a if !_divides(e, pivot)]
  push!(rhs, pivot)

  ### Assembly of the division ideal with less total degree
  lhs = _divide_by(a, pivot)

  f = one(S)
  for i in 1:nvars(S)
    z = t[i]
    f *= z^(sum([pivot[j]*weight_matrix[i, j] for j in 1:length(pivot)]))
  end

  return _hilbert_numerator_custom(rhs, weight_matrix, S, t) + f*_hilbert_numerator_custom(lhs, weight_matrix, S, t)
end

function _hilbert_numerator_bayer_stillman(
    a::Vector{Vector{Int}}, weight_matrix::Matrix{Int},
    S::Ring, t::Vector, oneS = one(S)
  )
  ###########################################################################
  # For this strategy see
  #
  # Bayer, Stillman: Computation of Hilber Series
  # J. Symbolic Computation (1992) No. 14, pp. 31--50
  #
  # Algorithm 2.6, page 35
  ###########################################################################
  ret = _hilbert_numerator_trivial_cases(a, weight_matrix, S, t, oneS)
  ret !== nothing && return ret

  # make sure we have lexicographically ordered monomials
  sort!(a, alg=QuickSort)

  # initialize the result 
  h = oneS - _expvec_to_poly(S, t, weight_matrix, a[1])
  linear_mons = Vector{Int}()
  for i in 2:length(a)
    J = _divide_by(a[1:i-1], a[i])
    empty!(linear_mons)
    J1 = filter!(J) do m
      k = findfirst(!iszero, m)
      k === nothing && return false # filter out zero vector
      if m[k] == 1 && findnext(!iszero, m, k + 1) === nothing
        push!(linear_mons, k)
        return false
      end
      return true
    end
    q = _hilbert_numerator_bayer_stillman(J1, weight_matrix, S, t, oneS)
    for k in linear_mons
      q *= (oneS - prod(t[i]^weight_matrix[i, k] for i in 1:length(t)))
    end

    h -= q * _expvec_to_poly(S, t, weight_matrix, a[i])
  end
  return h
end

############################################################################
### Homogenization and Dehomogenization
############################################################################


###old: special and to be applied with care
### needed for: AffinePlaneCurve.jl
### TODO: make adjustments there and omit function below
function homogenization(f::MPolyRingElem, S::MPolyDecRing, pos::Int = 1)
  d = total_degree(f)
  B = MPolyBuildCtx(S)
  for (c,e) = zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    insert!(e, pos, d-sum(e))
    push_term!(B, c, e)
  end
  return finish(B)
end

function _homogenization(f::MPolyRingElem, S::MPolyDecRing, start_pos::Int = 1)
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
   tdeg  = vcat(maximum(hcat(degs...), dims=2)...)
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

function _homogenization(f::MPolyRingElem, W::ZZMatrix, var::VarName, pos::Int = 1)
  R = parent(f)
  A = symbols(R)
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

function _homogenization(f::MPolyRingElem, W::Matrix{<:IntegerUnion}, var::VarName, pos::Int = 1)
   W = matrix(ZZ, W)
   return _homogenization(f, W, var, pos)
end

@doc raw"""
    homogenization(f::MPolyRingElem, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName, pos::Int = 1)

If $m$ is the number of rows of `W`, extend the parent polynomial ring of `f` by inserting $m$ extra variables, starting at position `pos`.
Correspondingly, extend the integer matrix `W` by inserting the standard unit vectors of size $m$ as new columns, starting at column `pos`.
Grade the extended ring by converting the columns of the extended matrix to elements of the group $\mathbb Z^m$ and assigning these
as weights to the variables. Homogenize `f` with respect to the induced $\mathbb Z^m$-grading on the original ring, using the extra variables as homogenizing
variables. Return the result as an element of the extended ring with its $\mathbb Z^m$-grading. If $m=1$, the extra variable prints as `var`. Otherwise,
the extra variables print as `var[`$i$`]`, for $i = 1 \dots m$.

    homogenization(V::Vector{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName, pos::Int = 1) where {T <: MPolyRingElem}

Given a vector `V` of elements in a common polynomial ring, create an extended ring with $\mathbb Z^m$-grading as above.
Homogenize the elements of `V` correspondingly,  and return the vector of homogenized elements.

    homogenization(I::MPolyIdeal{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName, pos::Int = 1) where {T <: MPolyRingElem}

Return the homogenization of `I` in an extended ring with $\mathbb Z^m$-grading as above.

!!! note
    Applied to an ideal `I`, the function first homogenizes the generators of `I` in the extended ring.
    It then creates the ideal generated by these homogenizations, and saturates this ideal
    with respect to the ideal which is generated by the homogenizing variables.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, QQMPolyRingElem[x, y])

julia> f = x^3+x^2*y+x*y^2+y^3
x^3 + x^2*y + x*y^2 + y^3

julia> W = [1 2; 3 4]
2×2 Matrix{Int64}:
 1  2
 3  4

julia> F = homogenization(f, W, "z", 3)
x^3*z[1]^3*z[2]^3 + x^2*y*z[1]^2*z[2]^2 + x*y^2*z[1]*z[2] + y^3

julia> parent(F)
Multivariate Polynomial Ring in x, y, z[1], z[2] over Rational Field graded by
  x -> [1 3]
  y -> [2 4]
  z[1] -> [1 0]
  z[2] -> [0 1]
```
"""
function homogenization(f::MPolyRingElem, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName, pos::Int = 1)
   return _homogenization(f, W, var, pos)
end

function homogenization(V::Vector{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName, pos::Int = 1) where {T <: MPolyRingElem}
  @assert all(x->parent(x) == parent(V[1]), V)
  R = parent(V[1])
  A = copy(symbols(R))
  l = length(A)
  @req pos in 1:l+1 "Index out of range."
  if size(W, 1) == 1
     insert!(A, pos, Symbol(var))
  else
     for i = 1:(size(W, 1))
       insert!(A, pos-1+i, _make_variable(var, i))
     end
  end
  L, _ = polynomial_ring(base_ring(R), A)
  G = abelian_group(zeros(Int, size(W, 1)))
  WH = hcat(Matrix(W[:, 1:(pos-1)]), Matrix(identity_matrix(ZZ, size(W, 1))[:, 1:size(W, 1)]), Matrix(W[:, pos:l]))
  S, _ = grade(L, [G(WH[:, i]) for i = 1:size(WH, 2)])
  l = length(V)
  return [_homogenization(V[i], S, pos) for i=1:l]
end

function homogenization(I::MPolyIdeal{T},  W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, var::VarName, pos::Int = 1) where {T <: MPolyRingElem}
  # TODO: Adjust for ZZ-gradings as soon as weighted orderings are available
  V = homogenization(gens(I), W, var, pos)
  R = parent(V[1])
  IH = ideal(R, V)
  Y = ideal(R, [gen(R, i) for i = pos:(pos+size(W, 1)-1)])
  return saturation(IH, Y)
end

@doc raw"""
    homogenization(f::MPolyRingElem, var::VarName, pos::Int = 1)

    homogenization(V::Vector{T}, var::VarName, pos::Int = 1) where {T <: MPolyRingElem}

    homogenization(I::MPolyIdeal{T}, var::VarName, pos::Int = 1; ordering::Symbol = :degrevlex) where {T <: MPolyRingElem}

Homogenize `f`, `V`, or `I` with respect to the standard $\mathbb Z$-grading using a homogenizing variable printing as `var`.
Return the result as an element of a graded polynomial ring with the homogenizing variable at position `pos`.

!!! note
    Applied to an ideal `I`, the function proceeds by homogenizing the elements of a Gröbner basis of `I` with respect to a degree compatible
    monomial ordering such as `degrevlex` (default).  If a Gröbner basis with respect to the specified ordering has not yet been computed and,
    thus, not yet been cached, executing the `homogenization` function with argument `I` may take some time.
    The degree compatibility of the specified ordering is not checked by the function.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, QQMPolyRingElem[x, y, z])

julia> f = x^3-y^2-z
x^3 - y^2 - z

julia> F = homogenization(f, "w", 4)
x^3 - y^2*w - z*w^2

julia> parent(F)
Multivariate Polynomial Ring in x, y, z, w over Rational Field graded by
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
 w*y - x^2
 w^2*z - x^3

julia> I = ideal(R, V)
ideal(-x^2 + y, -x^3 + z)

julia> PTC = homogenization(I, "w")
ideal(-x*z + y^2, -w*z + x*y, -w*y + x^2)

julia> parent(PTC[1])
Multivariate Polynomial Ring in w, x, y, z over Rational Field graded by
  w -> [1]
  x -> [1]
  y -> [1]
  z -> [1]

julia> homogenization(I, "w", ordering = deglex(gens(base_ring(I))))
ideal(x*z - y^2, -w*z + x*y, -w*y + x^2, -w*z^2 + y^3)
```
"""
function homogenization(f::MPolyRingElem, var::VarName, pos::Int = 1)
  R = parent(f)
  A = copy(symbols(R))
  l = length(A)
  @req pos in 1:l+1 "Index out of range."
  insert!(A, pos, Symbol(var))
  L, _ = polynomial_ring(base_ring(R), A)
  S, = grade(L)
  return _homogenization(f, S, pos)
end

function homogenization(V::Vector{T}, var::VarName, pos::Int = 1) where {T <: MPolyRingElem}
  @assert all(x->parent(x) == parent(V[1]), V)
  R = parent(V[1])
  A = copy(symbols(R))
  l = length(A)
  @req pos in 1:l+1 "Index out of range."
  insert!(A, pos, Symbol(var))
  L, _ = polynomial_ring(base_ring(R), A)
  S, = grade(L)
  l = length(V)
  return [_homogenization(V[i], S, pos) for i=1:l]
end

function homogenization(I::MPolyIdeal{T}, var::VarName, pos::Int = 1; ordering::MonomialOrdering = default_ordering(base_ring(I))) where {T <: MPolyRingElem}
  # TODO: Adjust as soon as new GB concept is implemented
  return ideal(homogenization(gens(groebner_basis(I, ordering=ordering)), var, pos))
end

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
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> F = x^3-x^2*y-x*z^2
x^3 - x^2*y - x*z^2

julia> f = dehomogenization(F, 1)
-y - z^2 + 1

julia> parent(f)
Multivariate Polynomial Ring in y, z over Rational Field

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
(Multivariate Polynomial Ring in w, x, y, z over Rational Field graded by
  w -> [1 3]
  x -> [2 4]
  y -> [1 0]
  z -> [0 1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w, x, y, z])

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
