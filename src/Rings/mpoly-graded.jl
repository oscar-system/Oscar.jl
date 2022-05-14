export weight, decorate, is_homogeneous, homogeneous_components, filtrate,
grade, GradedPolynomialRing, homogeneous_component, jacobi_matrix, jacobi_ideal,
HilbertData, hilbert_series, hilbert_series_reduced, hilbert_series_expanded, hilbert_function, hilbert_polynomial, grading,
homogenization, dehomogenization, grading_group, is_z_graded, is_zm_graded, is_positively_graded, is_standard_graded
export MPolyRing_dec, MPolyElem_dec, is_homogeneous, isgraded

@attributes mutable struct MPolyRing_dec{T, S} <: AbstractAlgebra.MPolyRing{T}
  R::S
  D::GrpAbFinGen
  d::Vector{GrpAbFinGenElem}
  lt
  function MPolyRing_dec(R::S, d::Vector{GrpAbFinGenElem}) where {S}
    @assert length(d) == ngens(R)
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    return r
  end
  function MPolyRing_dec(R::S, d::Vector{GrpAbFinGenElem}, lt) where {S}
    @assert length(d) == ngens(R)
    r = new{elem_type(base_ring(R)), S}()
    r.R = R
    r.D = parent(d[1])
    r.d = d
    r.lt = lt
    return r
  end
end

@doc Markdown.doc"""
    grading_group(R::MPolyRing_dec)

If `R` is, say, `G`-graded, return `G`.

# Examples
```jldoctest
julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> grading_group(R)
GrpAb: Z
```
"""
grading_group(R::MPolyRing_dec) = R.D

@doc Markdown.doc"""

    isgraded(R::MPolyRing_dec)

Return `true` if `R` is graded, `false` otherwise.
"""
function isgraded(R::MPolyRing_dec)
   return !isdefined(R, :lt)
end
isfiltered(R::MPolyRing_dec) = isdefined(R, :lt)

function show(io::IO, W::MPolyRing_dec)
  Hecke.@show_name(io, W)
  Hecke.@show_special(io, W)
  if isfiltered(W)
    println(io, "$(W.R) filtrated by ")
  else
    println(io, "$(W.R) graded by ")
  end
  R = W.R
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
  S = MPolyRing_dec(R, [1*A[1] for i = 1: ngens(R)], (x,y) -> x[1] < y[1])
  return S, map(R, gens(R))
end

@doc Markdown.doc"""
    grade(R::MPolyRing, W::Vector{<:IntegerUnion})

Given a vector `W` of `ngens(R)` integers, define a $\mathbb Z$-grading on `R` by creating 
a free abelian group of type `GrpAbFinGen` given by one free generator, converting the entries 
of `W` to elements of that group, and assigning these elements as weights 
to the variables. Return the graded ring as an object of type `MPolyRing_dec`, 
together with the vector of variables.

    grade(R::MPolyRing)

As above, where the grading is the standard $\mathbb Z$-grading on `R`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> T, (x, y, z) = grade(R)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])
```
"""
function grade(R::MPolyRing, W::Vector{<:IntegerUnion})
  @assert length(W) == ngens(R)
  A = abelian_group([0])
  set_attribute!(A, :show_elem => show_special_elem_grad) 
  S = MPolyRing_dec(R, [i*A[1] for i = W])
  return S, map(S, gens(R))
end

function grade(R::MPolyRing)
  A = abelian_group([0])
  S = MPolyRing_dec(R, [1*A[1] for i = 1: ngens(R)])
  return S, map(S, gens(R))
end

@doc Markdown.doc"""
    grade(R::MPolyRing, W::Vector{<:Vector{<:IntegerUnion}})
 
Given a vector `W` of `ngens(R)` integer vectors of the same size `m`, say, define a $\mathbb Z^m$-grading on `R` 
by creating a free abelian group of type `GrpAbFinGen` given by `m` free generators, converting the vectors in
`W` to elements of that group, and assigning these elements as weights to the variables. Return the graded
ring as an object of type `MPolyRing_dec`, together with the vector of variables.

    grade(R::MPolyRing, W::Union{fmpz_mat, Matrix{<:IntegerUnion}})

As above, converting the columns of `W`.

# Examples
```jldoctest
julia> R, x, y = PolynomialRing(QQ, "x" => 1:2, "y" => 1:3)
(Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field, fmpq_mpoly[x[1], x[2]], fmpq_mpoly[y[1], y[2], y[3]])

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
  y[3] -> [0 1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], y[1], y[2], y[3]])
```
"""
function grade(R::MPolyRing, W::Vector{<:Vector{<:IntegerUnion}})
  @assert length(W) == ngens(R)
  n = length(W[1])
  @assert all(x->length(x) == n, W)

  A = abelian_group(zeros(Int, n))
  set_attribute!(A, :show_elem => show_special_elem_grad) 
  S = MPolyRing_dec(R, [A(w) for w = W])
  return S, map(S, gens(R))
end

function grade(R::MPolyRing, W::Union{fmpz_mat, Matrix{<:IntegerUnion}})
  @assert size(W, 2) == ngens(R)
  A = abelian_group(zeros(Int, size(W, 1)))
  set_attribute!(A, :show_elem => show_special_elem_grad) 
  ###S = MPolyRing_dec(R, [A(view(W, :, i)) for i = 1:size(W, 2)])
  S = MPolyRing_dec(R, [A(W[:, i]) for i = 1:size(W, 2)])
  return S, map(S, gens(R))
end



@doc Markdown.doc"""
    is_standard_graded(R::MPolyRing_dec)

Return `true` if `R` is standard $\mathbb Z$-graded, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> is_standard_graded(S)
false
```
"""
function is_standard_graded(R::MPolyRing_dec)
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

@doc Markdown.doc"""
    is_z_graded(R::MPolyRing_dec)

Return `true` if `R` is $\mathbb Z$-graded, `false` otherwise.

!!! note
    Writing `G = grading_group(R)`, we say that `R` is $\mathbb Z$-graded if `is_free(G) && ngens(G) == rank(G) == 1` evaluates to `true`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> W = [1, 2, 3];

julia> S, (x, y, z) = grade(R, W)
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> is_z_graded(S)
true
```
"""
function is_z_graded(R::MPolyRing_dec)
  isgraded(R) || return false
  A = grading_group(R)
  return ngens(A) == 1 && rank(A) == 1 && is_free(A)
end

@doc Markdown.doc"""
    is_zm_graded(R::MPolyRing_dec)

Return `true` if `R` is $\mathbb Z^m$-graded for some $m$, `false` otherwise.

!!! note
    Writing `G = grading_group(R)`, we say that `R` is $\mathbb Z^m$-graded if `is_free(G) && ngens(G) == rank(G) == m` evaluates to `true`.

# Examples
```jldoctest
julia> R, x = PolynomialRing(QQ, "x" => 1:5)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, fmpq_mpoly[x[1], x[2], x[3], x[4], x[5]])

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
  x[5] -> [1 1 0 0], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4], x[5]])

julia> is_zm_graded(S)
false

julia> G = abelian_group(fmpz_mat([1 -1]));

julia> g = gen(G, 1)
Element of
(General) abelian group with relation matrix
[1 -1]
with components [0 1]

julia> W = [g, g, g, g];

julia> R, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"], W);

julia> is_free(G)
true

julia> is_zm_graded(R)
false
```
"""
function is_zm_graded(R::MPolyRing_dec)
  isgraded(R) || return false
  A = grading_group(R)
  return is_free(A) && ngens(A) == rank(A)
end

@doc Markdown.doc"""
    is_positively_graded(R::MPolyRing_dec)

Return `true` if `R` is positively graded, `false` otherwise.

!!! note
    We say that `R` is *positively graded* by a finitely generated abelian group $G$ if the coefficient ring of `R` is a field, 
    $G$ is free, and each graded part $R_g$, $g\in G$, has finite dimension.

# Examples
```jldoctest
julia> R, (t, x, y) = PolynomialRing(QQ, ["t", "x", "y"])
(Multivariate Polynomial Ring in t, x, y over Rational Field, fmpq_mpoly[t, x, y])

julia> G = abelian_group([0])
GrpAb: Z

julia> S, (t, x, y) = grade(R, [-1, 1, 1])
(Multivariate Polynomial Ring in t, x, y over Rational Field graded by
  t -> [-1]
  x -> [1]
  y -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[t, x, y])

julia> is_positively_graded(S)
false

julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

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
  y -> [1 0], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y])

julia> is_positively_graded(S)
false
```
"""
function is_positively_graded(R::MPolyRing_dec)
  isgraded(R) || return false
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

@doc Markdown.doc"""
    GradedPolynomialRing(C::Ring, V::Vector{String}, W; ordering=:lex)

Create a multivariate polynomial ring with coefficient ring `C` and variables which
print according to the strings in `V`, and grade this ring according to the
data provided by `W` (see the documentation of the `grade`-function for what is
possible). Return the graded ring as an object of type `MPolyRing_dec`, together 
with the vector of variables.

    GradedPolynomialRing(C::Ring, V::Vector{String}; ordering=:lex)

As above, where the grading is the standard $\mathbb Z$-grading. 

# Examples
```jldoctest
julia> W = [[1, 0], [0, 1], [1, 0], [4, 1]]
4-element Vector{Vector{Int64}}:
 [1, 0]
 [0, 1]
 [1, 0]
 [4, 1]

julia> R, x = GradedPolynomialRing(QQ, ["x[1]", "x[2]", "x[3]", "x[4]"], W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4] over Rational Field graded by 
  x[1] -> [1 0]
  x[2] -> [0 1]
  x[3] -> [1 0]
  x[4] -> [4 1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4]])

julia> S, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> T, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])
```
"""
function GradedPolynomialRing(C::Ring, V::Vector{String}, W; ordering=:lex)
   return grade(PolynomialRing(C, V; ordering = ordering)[1], W)
end
function GradedPolynomialRing(C::Ring, V::Vector{String}; ordering=:lex)
   W = ones(Int, length(V))
   return GradedPolynomialRing(C, V, W; ordering = ordering)
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
  S = MPolyRing_dec(R, [i*A[1] for i = v], (x,y) -> x[1] < y[1])
  return S, map(S, gens(R))
end

function filtrate(R::MPolyRing, v::Vector{GrpAbFinGenElem}, lt)
  S = MPolyRing_dec(R, v, lt)
  return S, map(S, gens(R))
end

@doc Markdown.doc"""
    grade(R::MPolyRing, W::Vector{GrpAbFinGenElem})

Given a vector `W` of `ngens(R)` elements of a group `G` of type `GrpAbFinGen`,  
define a  `G`-grading on `R` by assigning weights to the variables according to the entries of `W`.
Return the graded ring as an object of type `MPolyRing_dec`, together with the
vector of variables.

# Examples
```jldoctest
julia> R, (t, x, y) = PolynomialRing(QQ, ["t", "x", "y"])
(Multivariate Polynomial Ring in t, x, y over Rational Field, fmpq_mpoly[t, x, y])

julia> typeof(R)
FmpqMPolyRing

julia>  typeof(x)
fmpq_mpoly

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
  y -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[t, x, y])

julia> typeof(S)
MPolyRing_dec{fmpq, FmpqMPolyRing}

julia> S isa MPolyRing
true

julia> typeof(x)
MPolyElem_dec{fmpq, fmpq_mpoly}

julia> R, x, y = PolynomialRing(QQ, "x" => 1:2, "y" => 1:3)
(Multivariate Polynomial Ring in x[1], x[2], y[1], y[2], y[3] over Rational Field, fmpq_mpoly[x[1], x[2]], fmpq_mpoly[y[1], y[2], y[3]])

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
  y[3] -> [0 1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], y[1], y[2], y[3]])

julia> typeof(x[1])
fmpq_mpoly

julia> x = map(S, x)
2-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x[1]
 x[2]

julia> y = map(S, y)
3-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 y[1]
 y[2]
 y[3]

julia> typeof(x[1])
MPolyElem_dec{fmpq, fmpq_mpoly}

julia> R, x = PolynomialRing(QQ, "x" => 1:5)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, fmpq_mpoly[x[1], x[2], x[3], x[4], x[5]])

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
  x[5] -> [1 1 0 0], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4], x[5]])
```
"""
function grade(R::MPolyRing, v::Vector{GrpAbFinGenElem})
  S = MPolyRing_dec(R, v)
  return S, map(S, gens(R))
end

mutable struct MPolyElem_dec{T, S} <: MPolyElem{T}
  f::S
  parent
  function MPolyElem_dec(f::S, p) where {S}
    r = new{elem_type(base_ring(f)), S}(f, p)
#    if isgraded(p) && length(r) > 1
#      if !is_homogeneous(r)
#        error("element not homogeneous")
#      end
#both wrong and undesired.
#    end
    return r
  end
end

function Base.deepcopy_internal(f::MPolyElem_dec{T, S}, dict::IdDict) where {T, S}
  return MPolyElem_dec(Base.deepcopy_internal(f.f, dict), f.parent)
end

function show(io::IO, w::MPolyElem_dec)
  show(io, w.f)
end

parent(a::MPolyElem_dec{T, S}) where {T, S} = a.parent::MPolyRing_dec{T, parent_type(S)}

Nemo.symbols(R::MPolyRing_dec) = symbols(R.R)
Nemo.nvars(R::MPolyRing_dec) = nvars(R.R)

elem_type(::MPolyRing_dec{T, S}) where {T, S} = MPolyElem_dec{T, elem_type(S)}
elem_type(::Type{MPolyRing_dec{T, S}}) where {T, S} = MPolyElem_dec{T, elem_type(S)}
parent_type(::Type{MPolyElem_dec{T, S}}) where {T, S} = MPolyRing_dec{T, parent_type(S)}

(W::MPolyRing_dec)() = MPolyElem_dec(W.R(), W)
(W::MPolyRing_dec)(i::Int) = MPolyElem_dec(W.R(i), W)
(W::MPolyRing_dec)(f::Singular.spoly) = MPolyElem_dec(W.R(f), W)

function (W::MPolyRing_dec{S, T})(f::U) where {S, T, U}
  if parent_type(U) === T
    @assert W.R === parent(f)
    return MPolyElem_dec(f, W)
  else
    return W(W.R(f))
  end
end

function (W::MPolyRing_dec{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T
  return W(W.R(c, e))
end

(W::MPolyRing_dec)(g::MPolyElem_dec) = MPolyElem_dec(g.f, W)
one(W::MPolyRing_dec) = MPolyElem_dec(one(W.R), W)
zero(W::MPolyRing_dec) = MPolyElem_dec(zero(W.R), W)

################################################################################
#
#  Binary operations
#
################################################################################

for T in [:(+), :(-), :(*), :divexact]
  @eval ($T)(a::MPolyElem_dec,
             b::MPolyElem_dec) = MPolyElem_dec($T(a.f, b.f), parent(a))
end

################################################################################
#
#  Unitary operations
#
################################################################################

-(a::MPolyElem_dec)   = MPolyElem_dec(-a.f, parent(a))

################################################################################
#
#  Binary ad hoc operations
#
################################################################################

divexact(a::MPolyElem_dec, b::RingElem) = MPolyElem_dec(divexact(a.f, b), parent(a))

divexact(a::MPolyElem_dec, b::Integer) = MPolyElem_dec(divexact(a.f, b), parent(a))

divexact(a::MPolyElem_dec, b::Rational) = MPolyElem_dec(divexact(a.f, b), parent(a))

for T in [:(-), :(+)]
  @eval ($T)(a::MPolyElem_dec,
             b::RingElem) = MPolyElem_dec($(T)(a.f, b), parent(a))

  @eval ($T)(a::MPolyElem_dec,
             b::Integer) = MPolyElem_dec($(T)(a.f, b), parent(a))

  @eval ($T)(a::MPolyElem_dec,
             b::Rational) = MPolyElem_dec($(T)(a.f, b), parent(a))

  @eval ($T)(a::RingElem,
             b::MPolyElem_dec) = MPolyElem_dec($(T)(a, b.f), b.parent)

  @eval ($T)(a::Integer,
             b::MPolyElem_dec) = MPolyElem_dec($(T)(a, b.f), b.parent)

  @eval ($T)(a::Rational,
             b::MPolyElem_dec) = MPolyElem_dec($(T)(a, b.f), b.parent)
end

function *(a::MPolyElem_dec{T, S}, b::T) where {T <: RingElem, S}
  return MPolyElem_dec(a.f * b, parent(a))
end

function *(a::T, b::MPolyElem_dec{T, S}) where {T <: RingElem, S}
  return b * a
end

################################################################################
#
#  Factoring, division, ...
#
################################################################################

function factor(x::Oscar.MPolyElem_dec)
  R = parent(x)
  D = Dict{elem_type(R), Int64}()
  F = factor(x.f)
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

function gcd(x::Oscar.MPolyElem_dec, y::Oscar.MPolyElem_dec)
  R = parent(x)
  return R(gcd(x.f, y.f))
end

function div(x::Oscar.MPolyElem_dec, y::Oscar.MPolyElem_dec)
  R = parent(x)
  return R(div(x.f, y.f))
end

function divrem(x::MPolyElem_dec, y::MPolyElem_dec)
  R = parent(x)
  q, r = divrem(x.f, y.f)
  return R(q), R(r)
end

################################################################################
#
#  Equality
#
################################################################################

==(a::MPolyElem_dec, b::MPolyElem_dec) = a.f == b.f

^(a::MPolyElem_dec, i::Int) = MPolyElem_dec(a.f^i, parent(a))

function mul!(a::MPolyElem_dec, b::MPolyElem_dec, c::MPolyElem_dec)
  return b*c
end

function addeq!(a::MPolyElem_dec, b::MPolyElem_dec)
  return a+b
end

length(a::MPolyElem_dec) = length(a.f)
total_degree(a::MPolyElem_dec) = total_degree(a.f)
monomial(a::MPolyElem_dec, i::Int) = parent(a)(monomial(a.f, i))
coeff(a::MPolyElem_dec, i::Int) = coeff(a.f, i)
term(a::MPolyElem_dec, i::Int) = parent(a)(term(a.f, i))
exponent_vector(a::MPolyElem_dec, i::Int) = exponent_vector(a.f, i)
exponent_vector(a::MPolyElem_dec, i::Int, ::Type{T}) where T = exponent_vector(a.f, i, T)
exponent(a::MPolyElem_dec, i::Int, j::Int) = exponent(a.f, i, j)
exponent(a::MPolyElem_dec, i::Int, j::Int, ::Type{T}) where T = exponent(a.f, i, j, T)

function has_weighted_ordering(R::MPolyRing_dec)
  grading_to_ordering = false
  w_ord = degrevlex(gens(R.R)) # dummy, not used
  # This is not meant to be exhaustive, there a probably more gradings which one
  # can meaningfully translate into a monomial ordering
  # However, we want to stick to global orderings.
  if is_z_graded(R)
    w = Int[ R.d[i].coeff[1] for i = 1:ngens(R) ]
    if all(isone, w)
      w_ord = degrevlex(gens(R.R))
      grading_to_ordering = true
    elseif all(x -> x > 0, w)
      w_ord = wdegrevlex(gens(R.R), w)
      grading_to_ordering = true
    end
  end
  return grading_to_ordering, w_ord
end

function default_ordering(R::MPolyRing_dec)
  fl, w_ord = has_weighted_ordering(R)
  if fl
    return w_ord
  end
  return degrevlex(gens(R))
end

function singular_poly_ring(R::MPolyRing_dec; keep_ordering::Bool = false)
  if !keep_ordering
    return singular_poly_ring(R.R, default_ordering(R).o)
  end
  return singular_poly_ring(R.R, keep_ordering = keep_ordering)
end

MPolyCoeffs(f::MPolyElem_dec) = MPolyCoeffs(f.f)
MPolyExponentVectors(f::MPolyElem_dec) = MPolyExponentVectors(f.f)

function push_term!(M::MPolyBuildCtx{<:MPolyElem_dec{T, S}}, c::T, expv::Vector{Int}) where {T <: RingElement, S}
  if iszero(c)
    return M
  end
  len = length(M.poly.f) + 1
  set_exponent_vector!(M.poly.f, len, expv)
  setcoeff!(M.poly.f, len, c)
  return M
end

function set_exponent_vector!(f::MPolyElem_dec, i::Int, exps::Vector{Int})
  f.f = set_exponent_vector!(f.f, i, exps)
  return f
end

function finish(M::MPolyBuildCtx{<:MPolyElem_dec})
  f = sort_terms!(M.poly.f)
  f = combine_like_terms!(M.poly.f)
  return parent(M.poly)(f)
end

# constructor for ideals#######################################################

function ideal(g::Vector{T}) where {T <: MPolyElem_dec}
  @assert length(g) > 0
  @assert all(x->parent(x) == parent(g[1]), g)
  if isgraded(parent(g[1]))
     if !(all(is_homogeneous, g))
       throw(ArgumentError("The generators of the ideal must be homogeneous."))
     end
  end
  return MPolyIdeal(g)
end

function jacobi_matrix(f::MPolyElem_dec)
  R = parent(f)
  n = nvars(R)
  return matrix(R, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_ideal(f::MPolyElem_dec)
  R = parent(f)
  n = nvars(R)
  return ideal(R, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyElem_dec})
  R = parent(g[1])
  n = nvars(R)
  @assert all(x->parent(x) === R, g)
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

@doc Markdown.doc"""
    degree(f::MPolyElem_dec)

Given a homogeneous element `f` of a graded multivariate ring, return the degree of `f`.

    degree(::Type{Vector{Int}}, f::MPolyElem_dec)

Given a homogeneous element `f` of a $\mathbb Z^m$-graded multivariate ring, return the degree of `f`, converted to a vector of integer numbers.

    degree(::Type{Int}, f::MPolyElem_dec)

Given a homogeneous element `f` of a $\mathbb Z$-graded multivariate ring, return the degree of `f`, converted to an integer number.

# Examples
```jldoctest
julia> R, x = PolynomialRing(QQ, "x" => 1:5)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, fmpq_mpoly[x[1], x[2], x[3], x[4], x[5]])

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
  x[5] -> [1 1 0 0], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4], x[5]])

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

julia> R, x = GradedPolynomialRing(QQ, ["x[1]", "x[2]", "x[3]", "x[4]"], W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4] over Rational Field graded by
  x[1] -> [1 0]
  x[2] -> [0 1]
  x[3] -> [1 0]
  x[4] -> [4 1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^4*x[2]+x[4]
x[1]^4*x[2] + x[4]

julia> degree(f)
graded by [4 1]

julia> degree(Vector{Int}, f)
2-element Vector{Int64}:
 4
 1

julia>  R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

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
function degree(a::MPolyElem_dec)
  @req !iszero(a) "Element must be non-zero"
  W = parent(a)
  w = W.D[0]
  first = true
  d = W.d
  for c = MPolyExponentVectors(a.f)
    u = W.D[0]
    for i=1:length(c)
      u += c[i]*d[i]
    end
    if first
      first = false
      w = u
    elseif isfiltered(W)
      w = W.lt(w, u) ? u : w
    else
      w == u || error("element not homogeneous")
    end
  end
  return w
end

function degree(::Type{Int}, a::MPolyElem_dec)
  @assert is_z_graded(parent(a))
  return Int(degree(a)[1])
end

function degree(::Type{Vector{Int}}, a::MPolyElem_dec)
  @assert is_zm_graded(parent(a))
  d = degree(a)
  return Int[d[i] for i=1:ngens(parent(d))]
end

@doc Markdown.doc"""
    is_homogeneous(f::MPolyElem_dec)

Given an element `f` of a graded multivariate ring, return `true` if `f` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> f = x^2+y*z
x^2 + y*z

julia> is_homogeneous(f)
false

julia> W = [1 2 1 0; 3 4 0 1]
2×4 Matrix{Int64}:
 1  2  1  0
 3  4  0  1

julia> S, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"], W)
(Multivariate Polynomial Ring in w, x, y, z over Rational Field graded by
  w -> [1 3]
  x -> [2 4]
  y -> [1 0]
  z -> [0 1], MPolyElem_dec{fmpq, fmpq_mpoly}[w, x, y, z])

julia> F = w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3
w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3

julia> is_homogeneous(F)
true
```
"""
function is_homogeneous(F::MPolyElem_dec)
  D = parent(F).D
  d = parent(F).d
  S = Set{elem_type(D)}()
  for c = MPolyExponentVectors(F.f)
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

@doc Markdown.doc"""
    homogeneous_components(f::MPolyElem_dec{T, S}) where {T, S}

Given an element `f` of a graded multivariate ring, return the homogeneous components of `f`.

# Examples
```jldoctest
julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> f = x^2+y+z
x^2 + y + z

julia> homogeneous_components(f)
Dict{GrpAbFinGenElem, MPolyElem_dec{fmpq, fmpq_mpoly}} with 2 entries:
  [2] => x^2 + y
  [3] => z

julia> R, x = PolynomialRing(QQ, "x" => 1:5)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, fmpq_mpoly[x[1], x[2], x[3], x[4], x[5]])

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
  x[5] -> [1 1 0 0], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4], x[5]])

julia> f = x[1]^2+x[3]^2+x[5]^2
x[1]^2 + x[3]^2 + x[5]^2

julia> homogeneous_components(f)
Dict{GrpAbFinGenElem, MPolyElem_dec{fmpq, fmpq_mpoly}} with 2 entries:
  [2 2 0 0] => x[5]^2
  [2 0 0 0] => x[1]^2 + x[3]^2
```
"""
function homogeneous_components(a::MPolyElem_dec{T, S}) where {T, S}
  D = parent(a).D
  d = parent(a).d
  h = Dict{elem_type(D), typeof(a)}()
  W = parent(a)
  R = W.R
  # First assemble the homogeneous components into the build contexts.
  # Afterwards compute the polynomials.
  hh = Dict{elem_type(D), MPolyBuildCtx{S, DataType}}()
  dmat = vcat([d[i].coeff for i in 1:length(d)])
  tmat = zero_matrix(ZZ, 1, nvars(R))
  res_mat = zero_matrix(ZZ, 1, ncols(dmat))
  for (c, e) = Base.Iterators.zip(coefficients(a.f), exponent_vectors(a.f))
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

@doc Markdown.doc"""
    homogeneous_component(f::MPolyElem_dec, g::GrpAbFinGenElem)

Given an element `f` of a graded multivariate ring, and given an element 
`g` of the grading group of that ring, return the homogeneous component of `f` of degree `g`.

    homogeneous_component(f::MPolyElem_dec, g::Vector{<:IntegerUnion})

Given an element `f` of a $\mathbb  Z^m$-graded multivariate ring `R`, say, and given
a vector `g` of $m$ integers, convert `g` into an element of the grading group of `R`,
and return the homogeneous component of `f` whose degree is that element.

    homogeneous_component(f::MPolyElem_dec, g::IntegerUnion)

Given an element `f` of a $\mathbb  Z$-graded multivariate ring `R`, say, and given
an integer `g`, convert `g` into an element of the grading group of `R`, and return the 
homogeneous component of `f` whose degree is that element.

# Examples
```jldoctest
julia> R, x = PolynomialRing(QQ, "x" => 1:5)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4], x[5] over Rational Field, fmpq_mpoly[x[1], x[2], x[3], x[4], x[5]])

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
  x[5] -> [1 1 0 0], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4], x[5]])

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

julia> R, x = GradedPolynomialRing(QQ, ["x[1]", "x[2]", "x[3]", "x[4]"], W)
(Multivariate Polynomial Ring in x[1], x[2], x[3], x[4] over Rational Field graded by 
  x[1] -> [1 0]
  x[2] -> [0 1]
  x[3] -> [1 0]
  x[4] -> [4 1], MPolyElem_dec{fmpq, fmpq_mpoly}[x[1], x[2], x[3], x[4]])

julia> f = x[1]^2*x[2]+x[4]
x[1]^2*x[2] + x[4]

julia> homogeneous_component(f, [2, 1])
x[1]^2*x[2]

julia> R, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"], [1, 2, 3])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by 
  x -> [1]
  y -> [2]
  z -> [3], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

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
function homogeneous_component(a::MPolyElem_dec, g::GrpAbFinGenElem)
  R = parent(a).R
  r = R(0)
  d = parent(a).d
  for (c, m) = Base.Iterators.zip(MPolyCoeffs(a.f), Generic.MPolyMonomials(a.f))
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

function homogeneous_component(a::MPolyElem_dec, g::IntegerUnion)
  @assert is_z_graded(parent(a))
  return homogeneous_component(a, grading_group(parent(a))([g]))
end

function homogeneous_component(a::MPolyElem_dec, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(parent(a))
  return homogeneous_component(a, grading_group(parent(a))(g))
end

base_ring(W::MPolyRing_dec) = base_ring(W.R)
Nemo.ngens(W::MPolyRing_dec) = Nemo.ngens(W.R)
Nemo.gens(W::MPolyRing_dec) = map(W, gens(W.R))
Nemo.gen(W::MPolyRing_dec, i::Int) = W(gen(W.R, i))
Base.getindex(W::MPolyRing_dec, i::Int) = W(W.R[i])

base_ring(f::MPolyElem_dec) = base_ring(f.f)

function show_homo_comp(io::IO, M)
  (W, d) = get_attribute(M, :data)
  n = get_attribute(W, :name)
  if n != nothing
    print(io, "$(n)_$(d.coeff) of dim $(dim(M))")
  else
    println(io, "homogeneous component of $W of degree $d")
  end
end

@doc Markdown.doc"""
    homogeneous_component(R::MPolyRing_dec, g::GrpAbFinGenElem) 

Given a polynomial ring `R` over a field which is graded by a free
group of type `GrpAbFinGen`, and given an element `g` of that group,
return the homogeneous component of `R` of degree `g`. Additionally, return
the embedding of the component into `R`.

    homogeneous_component(R::MPolyRing_dec, g::Vector{<:IntegerUnion})

Given a $\mathbb  Z^m$-graded polynomial ring `R` over a field, and given
a vector `g` of $m$ integers, convert `g` into an element of the grading
group of `R`, and return the homogeneous component of `R` whose degree 
is that element. Additionally, return the embedding of the component into `R`.

    homogeneous_component(R::MPolyRing_dec, g::IntegerUnion)

Given a $\mathbb  Z$-graded polynomial ring `R` over a field, and given
an integer `g`, convert `g` into an element of the grading group of `R`, 
and return the homogeneous component of `R` whose degree is that element.
Additionally, return the embedding of the component into `R`.

!!! note
    If the component is not finite dimensional, an error message will be thrown.

# Examples
```jldoctest
julia> R, x, y = PolynomialRing(QQ, "x" => 1:2, "y" => 1:3);

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

julia> T, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

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
function homogeneous_component(W::MPolyRing_dec, d::GrpAbFinGenElem)
  #TODO: lazy: ie. no enumeration of points
  #      aparently it is possible to get the number of points faster than the points
  #TODO: in the presence of torsion, this is wrong. The component
  #      would be a module over the deg-0-sub ring.
  if !(coefficient_ring(W) isa AbstractAlgebra.Field)
       throw(ArgumentError("The coefficient ring must be a field."))
  end
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
       a = MPolyBuildCtx(W.R)
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

function homogeneous_component(R::MPolyRing_dec, g::Vector{<:IntegerUnion})
  @assert is_zm_graded(R)
  return homogeneous_component(R, grading_group(R)(g))
end

function homogeneous_component(R::MPolyRing_dec, g::IntegerUnion)
  @assert is_z_graded(R)
  return homogeneous_component(R, grading_group(R)([g]))
end

function vector_space(K::AbstractAlgebra.Field, e::Vector{T}; target = nothing) where {T <:MPolyElem}
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
    for (c, m) = Base.Iterators.zip(coefficients(i), monomials(i))
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
    for (c, m) = Base.Iterators.zip(coefficients(x), monomials(x))
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
function (W::MPolyRing_dec)(m::Generic.FreeModuleElem) 
  h = hasrelshp(parent(m), W)
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

function hasrelshp(R, S)
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
    W = R.d
    W = [Int(W[i][1]) for i = 1:ngens(R)]
    
    if !(minimum(W)>0)
       throw(ArgumentError("The weights must be positive."))
    end
    
    if !(coefficient_ring(R) isa AbstractAlgebra.Field)
       throw(ArgumentError("The coefficient ring must be a field."))
    end

    if !((R isa Oscar.MPolyRing_dec) && (isgraded(R)))
       throw(ArgumentError("The base ring must be graded."))
    end
    
    if !(all(is_homogeneous, gens(I)))
       throw(ArgumentError("The generators of the ideal must be homogeneous."))
    end
    
    G = groebner_assure(I)
    h = Singular.hilbert_series(G.S, W)
    return new(h, W, I)
  end
  function HilbertData(B::BiPolyArray)
    return HilbertData(Oscar.MPolyIdeal(B))
  end
end

function hilbert_series(H::HilbertData, i::Int= 1)
  Zt, t = ZZ["t"]
  den = prod([1-t^H.weights[i] for i = 1:ngens(base_ring(H.I))])
  if i==1   ### the Hilbert series with denominator prod (1-t^H.weights[i])
    return Zt(map(fmpz, H.data[1:end-1])), den
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

  for i = 1:ngens(base_ring(H.I))
     if H.weights[i] != 1
       throw(ArgumentError("All weights must be 1."))
     end
  end
  
  q, dn = hilbert_series(H, 2)
  a = fmpq[]
  nf = fmpq(1)
  d = degree(dn)-1
  for i=1:d+1
    push!(a, q(1)//nf)    
    nf *= i
    q = derivative(q)
  end
  Qt, t = QQ["t"]
  t = gen(Qt)
  bin = one(parent(t))
  b = fmpq_poly[]
  if d==-1 return zero(parent(t)) end
  for i=0:d
    push!(b, (-1)^(d-i)*a[d-i+1]*bin)
    bin *= (t+i+1)*fmpq(1, i+1)
  end
  return sum(b)
end

function Oscar.degree(H::HilbertData)

  for i = 1:ngens(base_ring(H.I))
     if H.weights[i] != 1
       throw(ArgumentError("All weights must be 1."))
     end
  end
  
  P = hilbert_polynomial(H)
  if P==zero(parent(P))
     q, _ = hilbert_series(H, 2)
     return q(1)
  end
  return leading_coefficient(P)*factorial(degree(P))
end

function (P::FmpqRelSeriesRing)(H::HilbertData)
  n, d = hilbert_series(H, 2)
  Qt, t = QQ["t"]
  nn = map_coefficients(QQ, n, parent = Qt)
  dd = map_coefficients(QQ, d, parent = Qt)
  gg, ee, _ = gcdx(dd, gen(Qt)^max_precision(P))
  @assert isone(gg)
  nn = Hecke.mullow(nn, ee, max_precision(P))
  c = collect(coefficients(nn))
  return P(map(fmpq, c), length(c), max_precision(P), 0)
end

function hilbert_series_expanded(H::HilbertData, d::Int)
   T, t = PowerSeriesRing(QQ, d+1, "t")   
   return T(H)
end

function hilbert_function(H::HilbertData, d::Int)
   if d<0 return 0 end
   HS = hilbert_series_expanded(H,d)
   return coeff(hilbert_series_expanded(H, d), d)
end

function Base.show(io::IO, h::HilbertData)
  print(io, "Hilbert Series for $(h.I), data: $(h.data), weights: $(h.weights)")  ###new
end


############################################################################
### Homogenization and Dehomogenization
############################################################################


###old: special and to be applied with care
### needed for: AffinePlaneCurve.jl
### TODO: make adjustments there and omitt fuction below
function homogenization(f::MPolyElem, S::MPolyRing_dec, pos::Int = 1)
  d = total_degree(f)
  B = MPolyBuildCtx(S)
  for (c,e) = zip(coefficients(f), exponent_vectors(f))
    insert!(e, pos, d-sum(e))
    push_term!(B, c, e)
  end
  return finish(B)
end

function _homogenization(f::MPolyElem, S::MPolyRing_dec, start_pos::Int = 1)
   len_gg  = ngens(S.D)
      #= get exponent vectors and enlarge them for S =#
   exps = Vector{Int}[]
   for e in exponent_vectors(f)
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

   for (i, (c, e)) = enumerate(zip(coefficients(f), exps))
       for j in 1:len_gg
           e[start_pos+j-1] = tdeg[j]-degs[i][j]
       end
       push_term!(F, c, e)
   end

   return finish(F)
end

function _homogenization(f::MPolyElem, W::fmpz_mat, var::String, pos::Int = 1)
  R = parent(f)
  A = String.(symbols(R))
  l = length(A)
  pos in 1:l+1 || throw(ArgumentError("Index out of range."))
  if size(W, 1) == 1
     insert!(A, pos, var)
  else
     for i = 1:(size(W, 1))
       insert!(A, pos-1+i, _make_variable(var, i))
     end
  end
  L, _ = PolynomialRing(R.base_ring, A)
  Lgens = gens(L)
  ###ff = evaluate(f, vcat(Lgens[1:pos-1], Lgens[pos + size(W, 1):end]))
  G = abelian_group(zeros(Int, size(W, 1)))
  WH = hcat(Matrix(W[:, 1:(pos-1)]), Matrix(identity_matrix(ZZ, size(W, 1))[:, 1:size(W, 1)]), Matrix(W[:, pos:l]))
  S, _ = grade(L, [G(WH[:, i]) for i = 1:size(WH, 2)])
  ###return _homogenization(ff, S, pos)
  return _homogenization(f, S, pos)
end

function _homogenization(f::MPolyElem, W::Matrix{<:IntegerUnion}, var::String, pos::Int = 1)
   W = matrix(ZZ, W)
   return _homogenization(f, W, var, pos)
end

@doc Markdown.doc"""
    homogenization(f::MPolyElem, W::Union{fmpz_mat, Matrix{<:IntegerUnion}}, var::String, pos::Int = 1)

If $m$ is the number of rows of `W`, extend the parent polynomial ring of `f` by inserting $m$ extra variables, starting at position `pos`.
Correspondingly, extend the integer matrix `W` by inserting the standard unit vectors of size $m$ as new columns, starting at column `pos`.
Grade the extended ring by converting the columns of the extended matrix to elements of the group $\mathbb Z^m$ and assigning these
as weights to the variables. Homogenize `f` with respect to the induced $\mathbb Z^m$-grading on the original ring, using the extra variables as homogenizing
variables. Return the result as an element of the extended ring with its $\mathbb Z^m$-grading. If $m=1$, the extra variable prints as `var`. Otherwise,
the extra variables print as `var[`$i$`]`, for $i = 1 \dots m$.

    homogenization(V::Vector{T},  W::Union{fmpz_mat, Matrix{<:IntegerUnion}}, var::String, pos::Int = 1) where {T <: MPolyElem}

Given a vector `V` of elements in a common polynomial ring, create an extended ring with $\mathbb Z^m$-grading as above.
Homogenize the elements of `V` correspondingly,  and return the vector of homogenized elements.

    homogenization(I::MPolyIdeal{T},  W::Union{fmpz_mat, Matrix{<:IntegerUnion}}, var::String, pos::Int = 1; ordering::Symbol = :degrevlex) where {T <: MPolyElem}

Return the homogenization of `I` in an extended ring with $\mathbb Z^m$-grading as above.

!!! note
    Applied to an ideal `I`, the function first homogenizes the generators of `I` in the extended ring.
    It then creates the ideal generated by these homogenizations, and saturates this ideal
    with respect to the ideal which is generated by the homogenizing variables.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

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
function homogenization(f::MPolyElem, W::Union{fmpz_mat, Matrix{<:IntegerUnion}}, var::String, pos::Int = 1)
   return _homogenization(f, W, var, pos)
end

function homogenization(V::Vector{T},  W::Union{fmpz_mat, Matrix{<:IntegerUnion}}, var::String, pos::Int = 1) where {T <: MPolyElem}
  @assert all(x->parent(x) == parent(V[1]), V)
  R = parent(V[1])
  A = String.(symbols(R))
  l = length(A)
  pos in 1:l+1 || throw(ArgumentError("Index out of range."))
  if size(W, 1) == 1
     insert!(A, pos, var)
  else
     for i = 1:(size(W, 1))
       insert!(A, pos-1+i, _make_variable(var, i))
     end
  end
  L, _ = PolynomialRing(R.base_ring, A)
  G = abelian_group(zeros(Int, size(W, 1)))
  WH = hcat(Matrix(W[:, 1:(pos-1)]), Matrix(identity_matrix(ZZ, size(W, 1))[:, 1:size(W, 1)]), Matrix(W[:, pos:l]))
  S, _ = grade(L, [G(WH[:, i]) for i = 1:size(WH, 2)])
  l = length(V)
  return [_homogenization(V[i], S, pos) for i=1:l]
end

function homogenization(I::MPolyIdeal{T},  W::Union{fmpz_mat, Matrix{<:IntegerUnion}}, var::String, pos::Int = 1) where {T <: MPolyElem}
  # TODO: Adjust for ZZ-gradings as soon as weighted orderings are available
  V = homogenization(gens(I), W, var, pos)
  R = parent(V[1])
  IH = ideal(R, V)
  Y = ideal(R, [gens(R)[i] for i = pos:(pos+size(W, 1)-1)])
  return saturation(IH, Y)
end

@doc Markdown.doc"""
    homogenization(f::MPolyElem, var::String, pos::Int = 1)

    homogenization(V::Vector{T}, var::String, pos::Int = 1) where {T <: MPolyElem}

    homogenization(I::MPolyIdeal{T}, var::String, pos::Int = 1; ordering::Symbol = :degrevlex) where {T <: MPolyElem}

Homogenize `f`, `V`, or `I` with respect to the standard $\mathbb Z$-grading using a homogenizing variable printing as `var`.
Return the result as an element of a graded polynomial ring with the homogenizing variable at position `pos`.

!!! note
    Applied to an ideal `I`, the function proceeds by homogenizing the elements of a Gröbner basis of `I` with respect to a degree compatible
    monomial ordering such as `degrevlex` (default).  If a Gröbner basis with respect to the specified ordering has not yet been computed and,
    thus, not yet been cached, executing the `homogenization` function with argument `I` may take some time.
    The degree compatibility of the specified ordering is not checked by the function.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

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
2-element Vector{fmpq_mpoly}:
 -x^2 + y
 -x^3 + z

julia> homogenization(V, "w")
2-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
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
function homogenization(f::MPolyElem, var::String, pos::Int = 1)
  R = parent(f)
  A = String.(symbols(R))
  l = length(A)
  pos in 1:l+1 || throw(ArgumentError("Index out of range."))
  insert!(A, pos, var)
  L, _ = PolynomialRing(R.base_ring, A)
  S, = grade(L)
  return _homogenization(f, S, pos)
end
function homogenization(V::Vector{T}, var::String, pos::Int = 1) where {T <: MPolyElem}
  @assert all(x->parent(x) == parent(V[1]), V)
  R = parent(V[1])
  A = String.(symbols(R))
  l = length(A)
  pos in 1:l+1 || throw(ArgumentError("Index out of range."))
  insert!(A, pos, var)
  L, _ = PolynomialRing(R.base_ring, A)
  S, = grade(L)
  l = length(V)
  return [_homogenization(V[i], S, pos) for i=1:l]
end
function homogenization(I::MPolyIdeal{T}, var::String, pos::Int = 1; ordering::MonomialOrdering = default_ordering(base_ring(I))) where {T <: MPolyElem}
  # TODO: Adjust as soon as new GB concept is implemented
  return ideal(homogenization(groebner_basis(I, ordering=ordering), var, pos))
end

### needed for: PlaneCurve-test.jl, ProjPlaneCurve.jl
### TODO: make adjustments there and omitt fuction below
function dehomogenization(F::MPolyElem_dec, R::MPolyRing, pos::Int)
  B = MPolyBuildCtx(R)
  for (c,e) = zip(coefficients(F), exponent_vectors(F))
    deleteat!(e, pos)
    push_term!(B, c, e)
  end
  return finish(B)
end



function _dehomogenization(F::MPolyElem_dec, R::MPolyRing, pos::Int, m::Int)
  B = MPolyBuildCtx(R)
  for (c,e) = zip(coefficients(F), exponent_vectors(F))
    deleteat!(e, pos:pos+m-1)
    push_term!(B, c, e)
  end
  return finish(B)
end


@doc Markdown.doc"""
    dehomogenization(F::MPolyElem_dec, pos::Int)

Given an element `F` of a $\mathbb Z^m$-graded ring, where the generators of $\mathbb Z^m$
are the assigned weights to the variables at positions `pos`, $\dots$, `pos` $-1+m$, dehomogenize
`F` using the variables at these positions. Return the result as an element of a polynomial ring not
depending on the variables at these positions.

    dehomogenization(V::Vector{T}, pos::Int) where {T <: MPolyElem_dec}

Given a vector `V` of elements in a common graded polynomial ring, create a polynomial ring not
depending on the variables at positions `pos`, $\dots$, `pos` $-1+m$. Dehomogenize the elements 
of `V` correspondingly,  and return the vector of dehomogenized elements.  

    dehomogenization(I::MPolyIdeal{T}, pos::Int) where {T <: MPolyElem_dec}

Return the dehomogenization of `I` in a polynomial ring as above.

# Examples
```jldoctest
julia> S, (x, y, z) = GradedPolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field graded by
  x -> [1]
  y -> [1]
  z -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y, z])

julia> F = x^3-x^2*y-x*z^2
x^3 - x^2*y - x*z^2

julia> f = dehomogenization(F, 1)
-y - z^2 + 1

julia> parent(f)
Multivariate Polynomial Ring in y, z over Rational Field

julia> V = [x*y-z^2, x^2*z-x^3]
2-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x*y - z^2
 -x^3 + x^2*z

julia> dehomogenization(V, 3)
2-element Vector{fmpq_mpoly}:
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

julia> S, (w, x, y, z) = GradedPolynomialRing(QQ, ["w", "x", "y", "z"], W)
(Multivariate Polynomial Ring in w, x, y, z over Rational Field graded by
  w -> [1 3]
  x -> [2 4]
  y -> [1 0]
  z -> [0 1], MPolyElem_dec{fmpq, fmpq_mpoly}[w, x, y, z])

julia> F = w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3
w^3*y^3*z^3 + w^2*x*y^2*z^2 + w*x^2*y*z + x^3

julia> dehomogenization(F, 3)
w^3 + w^2*x + w*x^2 + x^3
```
"""
function dehomogenization(F::MPolyElem_dec, pos::Int)
  S = parent(F)
  @assert is_zm_graded(S)
  l = ngens(S)
  G = grading_group(S)
  m = ngens(G)
  pos in 1:l-m+1 || throw(ArgumentError("Index out of range."))
  for i = 1:m
      @assert degree(gens(S)[pos-1+i]) == gen(G, i)
  end
  A = String.(symbols(S))
  deleteat!(A, pos:pos+m-1)
  R, _ = PolynomialRing(base_ring(S), A)
  return _dehomogenization(F, R, pos, m)
end
function dehomogenization(V::Vector{T}, pos::Int) where {T <: MPolyElem_dec}
  @assert all(x->parent(x) == parent(V[1]), V)
  S = parent(V[1])
  @assert is_zm_graded(S)
  l = ngens(S)
  G = grading_group(S)
  m = ngens(G)
  pos in 1:l-m+1 || throw(ArgumentError("Index out of range."))
  for i = 1:m
      @assert degree(gens(S)[pos-1+i]) == gen(G, i)
  end
  A = String.(symbols(S))
  deleteat!(A, pos:pos+m-1)
  R, _ = PolynomialRing(base_ring(S), A)
  l = length(V)
  return [_dehomogenization(V[i], R, pos, m) for i=1:l]
end
function dehomogenization(I::MPolyIdeal{T}, pos::Int) where {T <: MPolyElem_dec}
  return ideal(dehomogenization(gens(I), pos))
end

################################################################################
#
#  Evaluation
#
################################################################################

(f::MPolyElem_dec)(x...) = evaluate(f, collect(x))

################################################################################
#
#  Promote rule
#
################################################################################

function AbstractAlgebra.promote_rule(::Type{MPolyElem_dec{S, T}}, ::Type{U}) where {S, T, U}
  if AbstractAlgebra.promote_rule(T, U) === T
    return MPolyElem_dec{S, T}
  else
    return Union{}
  end
end
