

###############################################################################
# Graded Modules constructors
###############################################################################

@doc raw"""
    graded_free_module(R::Ring, p::Int, W::Vector{GrpAbFinGenElem}=[grading_group(R)[0] for i in 1:p], name::String="e")

Given a graded ring `R` with grading group `G`, say,
and given a vector `W` with `p` elements of `G`, create the free module $R^p$ 
equipped with its basis of standard unit vectors, and assign weights to these 
vectors according to the entries of `W`. Return the resulting graded free module.

    graded_free_module(R::Ring, W::Vector{GrpAbFinGenElem}, name::String="e")

As above, with `p = length(W)`.

!!! note
    The function applies to graded multivariate polynomial rings and their quotients.

The string `name` specifies how the basis vectors are printed. 

# Examples
```jldoctest
julia> R, (x,y) = graded_polynomial_ring(QQ, ["x", "y"])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> graded_free_module(R,3)
Graded free module R^3([0]) of rank 3 over R

julia> G = grading_group(R)
GrpAb: Z

julia> graded_free_module(R, [G[1], 2*G[1]])
Graded free module R^1([-1]) + R^1([-2]) of rank 2 over R
```
"""
function graded_free_module(R::Ring, p::Int, W::Vector{GrpAbFinGenElem}=[grading_group(R)[0] for i in 1:p], name::String="e")
  @assert length(W) == p
  @assert is_graded(R)
  all(x -> parent(x) == grading_group(R), W) || error("entries of W must be elements of the grading group of the base ring")
  M = FreeMod(R, p, name)
  M.d = W
  return M
end

function graded_free_module(R::Ring, p::Int, W::Vector{Any}, name::String="e")
  @assert length(W) == p
  @assert is_graded(R)
  p == 0 || error("W should be either an empty array or a Vector{GrpAbFinGenElem}")
  W = GrpAbFinGenElem[]
  return graded_free_module(R, p, W, name)
end

function graded_free_module(R::Ring, W::Vector{GrpAbFinGenElem}, name::String="e")
  p = length(W)
  return graded_free_module(R, p, W, name)
end

function graded_free_module(R::Ring, W::Vector{Any}, name::String="e")
  p = length(W)
  @assert is_graded(R)
  p == 0 || error("W should be either an empty array or a Vector{GrpAbFinGenElem}")
  W = GrpAbFinGenElem[]
  return graded_free_module(R, p, W, name)
end

@doc raw"""
    graded_free_module(R::Ring, W::Vector{<:Vector{<:IntegerUnion}}, name::String="e")

Given a graded ring `R` with grading group $G = \mathbb Z^m$, 
and given a vector `W` of integer vectors of the same size `p`, say, create the free 
module $R^p$ equipped with its basis of standard unit vectors, and assign weights to these 
vectors according to the entries of `W`, converted to elements of `G`. Return the 
resulting graded free module.

    graded_free_module(R::Ring, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, name::String="e")

As above, converting the columns of `W`.

    graded_free_module(R::Ring, W::Vector{<:IntegerUnion}, name::String="e")

Given a graded ring `R` with grading group $G = \mathbb Z$, 
and given a vector `W` of integers, set `p = length(W)`, create the free module $R^p$ 
equipped with its basis of standard unit vectors, and assign weights to these 
vectors according to the entries of `W`, converted to elements of `G`. Return 
the resulting graded free module.

The string `name` specifies how the basis vectors are printed. 

!!! note
    The function applies to graded multivariate polynomial rings and their quotients.

# Examples
```jldoctest
julia> R, (x,y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> F = graded_free_module(R, [1, 2])
Graded free module R^1([-1]) + R^1([-2]) of rank 2 over R
```

```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1 0 1; 0 1 1]);

julia> FF = graded_free_module(S, [[1, 2], [-1, 3]])
Graded free module S^1([-1 -2]) + S^1([1 -3]) of rank 2 over S

julia> FFF = graded_free_module(S, [1 -1; 2 3])
Graded free module S^1([-1 -2]) + S^1([1 -3]) of rank 2 over S

julia> FF == FFF
true
```
"""
function graded_free_module(R::Ring, W::Vector{<:Vector{<:IntegerUnion}}, name::String="e")
  @assert is_zm_graded(R)
  n = length(W[1])
  @assert all(x->length(x) == n, W)
  A = grading_group(R)
  d = [A(w) for w = W]
  return graded_free_module(R, length(W), d, name)
end

function graded_free_module(R::Ring, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, name::String="e")
  @assert is_zm_graded(R)
  A = grading_group(R)
  d = [A(W[:, i]) for i = 1:size(W, 2)]
  return graded_free_module(R, size(W, 2), d, name)
end

function graded_free_module(R::Ring, W::Vector{<:IntegerUnion}, name::String="e")
  @assert is_graded(R)
  A = grading_group(R)
  d = [W[i] * A[1] for i in 1:length(W)]
  return graded_free_module(R, length(W), d, name)
end

@doc raw"""
    grade(F::FreeMod, W::Vector{GrpAbFinGenElem})

Given a free module `F` over a graded ring with grading group `G`, say, and given
a vector `W` of `ngens(F)` elements of `G`, create a `G`-graded free module
by assigning the entries of `W` as weights to the generators of `F`. Return
the new module. 

    grade(F::FreeMod)

As above, with all weights set to `zero(G)`.

!!! note
    The function applies to free modules over both graded multivariate polynomial rings and their quotients.

# Examples
```jldoctest
julia> R, x, y = polynomial_ring(QQ, "x" => 1:2, "y" => 1:3);

julia> G = abelian_group([0, 0])
GrpAb: Z^2

julia> g = gens(G)
2-element Vector{GrpAbFinGenElem}:
 Element of G with components [1 0]
 Element of G with components [0 1]

julia> W = [g[1], g[1], g[2], g[2], g[2]];

julia> S, _ = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], y[1], y[2], y[3]])

julia> F = free_module(S, 3)
Free module of rank 3 over S

julia> FF = grade(F)
Graded free module S^3([0 0]) of rank 3 over S

julia> F
Free module of rank 3 over S
```
"""
function grade(F::FreeMod, W::Vector{GrpAbFinGenElem})
  @assert length(W) == ngens(F)
  @assert is_graded(base_ring(F))
  R = base_ring(F)
  all(x -> parent(x) == grading_group(R), W) || error("entries of W must be elements of the grading group of the base ring")
  N = free_module(R, length(W))
  N.d = W
  N.S = F.S
  return N
end

function grade(F::FreeMod)
  @assert is_graded(base_ring(F))
  R = base_ring(F)
  G = grading_group(R)
  W = [zero(G) for i = 1: ngens(F)]
  return grade(F, W)
end

@doc raw"""
    grade(F::FreeMod, W::Vector{<:Vector{<:IntegerUnion}})

Given a free module `F` over a graded ring with grading group $G = \mathbb Z^m$, and given
a vector `W` of `ngens(F)` integer vectors of the same size `m`, say, define a $G$-grading on `F` 
by converting the vectors in `W` to elements of $G$, and assigning these elements as weights to 
the variables. Return the new module.

    grade(F::FreeMod, W::Union{ZZMatrix, Matrix{<:IntegerUnion}})

As above, converting the columns of `W`.

    grade(F::FreeMod, W::Vector{<:IntegerUnion})

Given a free module `F` over a graded ring with grading group $G = \mathbb Z$, and given
a vector `W` of `ngens(F)` integers, define a $G$-grading on `F` converting the entries 
of `W` to elements of `G`, and assigning these elements as weights to the variables. 
Return the new module.

!!! note
    The function applies to free modules over both graded multivariate polynomial 
    rings and their quotients.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"],  [1 0 1; 0 1 1])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> F = free_module(R, 2)
Free module of rank 2 over R

julia> FF = grade(F,  [[1, 0], [0, 1]])
Graded free module R^1([-1 0]) + R^1([0 -1]) of rank 2 over R

julia> FFF = grade(F,  [1 0; 0 1])
Graded free module R^1([-1 0]) + R^1([0 -1]) of rank 2 over R
```

```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> S, _ = quo(R, [x*y])
(Quotient of multivariate polynomial ring by ideal(x*y), Map: graded multivariate polynomial ring -> quotient of multivariate polynomial ring)

julia> F = free_module(S, 2)
Free module of rank 2 over S

julia> FF = grade(F, [1, 2])
Graded free module S^1([-1]) + S^1([-2]) of rank 2 over S
```
"""
function grade(F::FreeMod, W::Vector{<:Vector{<:IntegerUnion}})
  @assert length(W) == ngens(F)
  R = base_ring(F)
  @assert is_zm_graded(R)
  n = length(W[1])
  @assert all(x->length(x) == n, W)
  A = grading_group(R)
  return grade(F, [A(w) for w = W])
end

function grade(F::FreeMod, W::Union{ZZMatrix, Matrix{<:IntegerUnion}})
  @assert size(W, 2) == ngens(F)
  R = base_ring(F)
  @assert is_zm_graded(R)
  A = grading_group(R)
  return grade(F, [A(W[:, i]) for i = 1:size(W, 2)])
end

function grade(F::FreeMod, W::Vector{<:IntegerUnion})
  @assert length(W) == ngens(F)
  R = base_ring(F)
  @assert is_z_graded(R)
  A = grading_group(R)
  N = free_module(R, length(W))
  N.d = [W[i] * A[1] for i in 1:length(W)]
  N.S = F.S
  return N
end

@doc raw"""
    grading_group(F::FreeMod)

Return the grading group of `base_ring(F)`.

# Examples
```jldoctest
julia> R, (x,y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> F = graded_free_module(R, 3)
Graded free module R^3([0]) of rank 3 over R

julia> grading_group(F)
GrpAb: Z
```
"""
function grading_group(M::FreeMod)
  return grading_group(base_ring(M))
end


# Dangerous: Only for internal use with care!!!
@doc raw"""
    set_grading!(F::FreeMod, W::Vector{GrpAbFinGenElem})

    set_grading!(F::FreeMod, W::Vector{<:Vector{<:IntegerUnion}})

    set_grading!(F::FreeMod, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}) 

    set_grading!(F::FreeMod, W::Vector{<:IntegerUnion})

Assign weights to the generators of `F` according to the entries of `W`.

See the `grade` and `graded_free_module` functions.
```
"""
function set_grading!(M::FreeMod, W::Vector{GrpAbFinGenElem})
  @assert length(W) == ngens(M)
  @assert is_graded(base_ring(M))
  R = base_ring(F)
  all(x -> parent(x) == grading_group(R), W) || error("entries of W must be elements of the grading group of the base ring")
  M.d = W
end

function set_grading!(M::FreeMod, W::Vector{<:Vector{<:IntegerUnion}})
  @assert length(W) == ngens(M)
  R = base_ring(M)
  @assert is_zm_graded(R)
  n = length(W[1])
  @assert all(x->length(x) == n, W)
  A = grading_group(R)
  M.d = [A(w) for w = W]
end

function set_grading!(M::FreeMod, W::Union{ZZMatrix, Matrix{<:IntegerUnion}})
  @assert size(W, 2) == ngens(M)
  R = base_ring(M)
  @assert is_zm_graded(R)
  A = grading_group(R)
  M.d = [A(W[:, i]) for i = 1:size(W, 2)]
end

function set_grading!(M::FreeMod, W::Vector{<:IntegerUnion})
  @assert length(W) == ngens(M)
  R = base_ring(M)
  @assert is_z_graded(R)
  A = grading_group(R)
  M.d = [W[i] * A[1] for i in 1:length(W)]
end

function degrees(M::FreeMod)
  @assert is_graded(M)
  return M.d
end

@doc raw"""
    degrees_of_generators(F::FreeMod)

Return the degrees of the generators of `F`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(R, 2)
Graded free module R^2([0]) of rank 2 over R

julia> degrees_of_generators(F)
2-element Vector{GrpAbFinGenElem}:
 [0]
 [0]
```
"""
function degrees_of_generators(F::FreeMod)
  return degrees(F)
end

###############################################################################
# Graded Free Modules functions
###############################################################################

function swap!(A::Vector{T}, i::Int, j::Int) where T
  A[i], A[j] = A[j], A[i]
end

function generate(k::Int, A::Vector{T}) where T
  if k == 1
    return [copy(A)]
  else
    perms = generate(k - 1, A)
    for i in 0:(k - 2)
      if k % 2 == 0
        swap!(A, i + 1, k)
      else
        swap!(A, 1, k)
      end
      perms = vcat(perms, generate(k - 1, A))
    end
    return perms
  end
end

function permute(v::Vector{T}) where T
  return generate(length(v), v)
end

function find_bijections(v_dict::Dict{T,Vector{Int}}, w_dict::Dict{T,Vector{Int}}, v_key::Int, bijections::Vector{Dict{Int,Int}}, current_bijection::Dict{Int,Int}) where T
  if v_key > length(keys(v_dict))
    push!(bijections, deepcopy(current_bijection))
    return nothing
  end
  element = collect(keys(v_dict))[v_key]
  v_indices = v_dict[element]
  w_indices = w_dict[element]
  if length(v_indices) == length(w_indices)
    for w_perm in permute(w_indices)
      next_bijection = deepcopy(current_bijection)
      for (i, j) in zip(v_indices, w_perm)
        next_bijection[i] = j
      end
      find_bijections(v_dict, w_dict, v_key + 1, bijections, next_bijection)
    end
  end
end

function get_multiset_bijection(
    v::Vector{T},
    w::Vector{T},
    all_bijections::Bool=false
) where {T<:Any}
  v_dict = Dict{T,Vector{Int}}()
  w_dict = Dict{T,Vector{Int}}()
  for (i, x) in enumerate(v)
    push!(get!(v_dict, x, []), i)
  end
  for (i, x) in enumerate(w)
    push!(get!(w_dict, x, []), i)
  end
  bijections = Vector{Dict{Int,Int}}()
  find_bijections(v_dict, w_dict, 1, bijections, Dict{Int,Int}())
  return all_bijections ? bijections : (isempty(bijections) ? nothing : bijections[1])
end

###############################################################################
# Graded Free Module elements functions
###############################################################################

@doc raw"""
    is_homogeneous(f::FreeModElem)

Given an element `f` of a graded free module, return `true` if `f` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"], [1, 2, 3]);

julia> F = free_module(R, 2)
Free module of rank 2 over R

julia> FF = grade(F, [1,4])
Graded free module R^1([-1]) + R^1([-4]) of rank 2 over R

julia> f = y^2*2*FF[1]-x*FF[2]
2*y^2*e[1] - x*e[2]

julia> is_homogeneous(f)
true
```
"""
function is_homogeneous(el::FreeModElem)
  !is_graded(parent(el)) && error("The parent module is not graded.")
  iszero(el) && return true
  el.d = isa(el.d, GrpAbFinGenElem) ? el.d : determine_degree_from_SR(coordinates(el), degrees(parent(el)))
  return isa(el.d, GrpAbFinGenElem)
end

@doc raw"""
    degree(f::FreeModElem)

Given a homogeneous element `f` of a graded free module, return the degree of `f`.

    degree(::Type{Vector{Int}}, f::FreeModElem)

Given a homogeneous element `f` of a $\mathbb Z^m$-graded free module, return the degree of `f`, converted to a vector of integer numbers.

    degree(::Type{Int}, f::FreeModElem)

Given a homogeneous element `f` of a $\mathbb Z$-graded free module, return the degree of `f`, converted to an integer number.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> f = y^2*z − x^2*w
-w*x^2 + y^2*z

julia> degree(f)
[3]

julia> typeof(degree(f))
GrpAbFinGenElem

julia> degree(Int, f)
3

julia> typeof(degree(Int, f))
Int64
```
"""
function degree(f::FreeModElem)
  !is_graded(parent(f)) && error("The parent module is not graded.")
  A = grading_group(base_ring(parent(f)))
  iszero(f) && return A[0]
  f.d = isa(f.d, GrpAbFinGenElem) ? f.d : determine_degree_from_SR(coordinates(f), degrees(parent(f)))
  isa(f.d, GrpAbFinGenElem) || error("The specified element is not homogeneous.")
  return f.d
end

function degree(::Type{Vector{Int}}, f::FreeModElem)
  @assert is_zm_graded(parent(f))
  d = degree(f)
  return Int[d[i] for i=1:ngens(parent(d))]
end

function degree(::Type{Int}, f::FreeModElem)
  @assert is_z_graded(parent(f))
  return Int(degree(f)[1])
end

function determine_degree_from_SR(coords::SRow, unit_vector_degrees::Vector{GrpAbFinGenElem})
  element_degree = nothing
  for (position, coordval) in coords
      if !is_homogeneous(coordval)
          return nothing
      end
      current_degree = degree(coordval) + unit_vector_degrees[position]
      if element_degree === nothing
          element_degree = current_degree
      elseif element_degree != current_degree
          return nothing
      end
  end
  return element_degree
end

###############################################################################
# Graded Free Module homomorphisms constructors
###############################################################################

function graded_map(A::MatElem)
  R = base_ring(A)
  G = grading_group(R)
  Fcdm = graded_free_module(R, [G[0] for _ in 1:ncols(A)])
  return graded_map(Fcdm, A)
end

function graded_map(F::FreeMod{T}, A::MatrixElem{T}) where {T <: RingElement}
  R = base_ring(F)
  G = grading_group(R)
  source_degrees = Vector{eltype(G)}()
  for i in 1:nrows(A)
      for j in 1:ncols(A)
          if A[i, j] != R[0]
              push!(source_degrees, degree(A[i, j]) + degree(F[j]))
              break
          end
      end
  end
  Fcdm = graded_free_module(R, source_degrees)
  phi = hom(Fcdm, F, A)
  return phi
end

function graded_map(F::FreeMod{T}, V::Vector{<:AbstractFreeModElem{T}}) where {T <: RingElement}
  R = base_ring(F)
  G = grading_group(R)
  nrows = length(V)
  ncols = rank(F)
  
  source_degrees = Vector{eltype(G)}()
  for i in 1:nrows
    for j in 1:ncols
      if coordinates(V[i])[j] != R[0]
        push!(source_degrees, degree(coordinates(V[i])[j]) + degree(F[j]))
        break
      end
    end
  end
  Fcdm = graded_free_module(R, source_degrees)
  phi = hom(Fcdm, F, V)
  return phi
end


function graded_map(F::SubquoModule{T}, V::Vector{<:ModuleFPElem{T}}) where {T <: RingElement}
  R = base_ring(F)
  G = grading_group(R)
  nrows = length(V)
  source_degrees = Vector{eltype(G)}()
  for i in 1:nrows
    for (j, coord_val) in coordinates(V[i])
      if coord_val != R[0]
        push!(source_degrees, degree(coord_val) + degree(F[j]))
        break
      end
    end
  end
  Fcdm = graded_free_module(R, source_degrees)
  phi = hom(Fcdm, F, V)
  return phi
end

###############################################################################
# Graded Free Module homomorphisms functions
###############################################################################

function set_grading(f::FreeModuleHom{T1, T2}) where {T1 <: FreeMod, T2 <: Union{FreeMod, SubquoModule, Oscar.SubModuleOfFreeModule}}
  if !is_graded(domain(f)) || !is_graded(codomain(f))
      return f
  end
  f.d = degree(f)
  return f
end

function set_grading(f::FreeModuleHom{T1, T2}) where {T1 <: FreeMod_dec, T2 <: FreeMod_dec}
  return f
end
# for decorations: add SubquoModule_dec for codomain once it exists

@doc raw"""
    degree(a::FreeModuleHom)

If `a` is graded, return the degree of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(R, 3)
Graded free module R^3([0]) of rank 3 over R

julia> G = graded_free_module(R, 2)
Graded free module R^2([0]) of rank 2 over R

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
3-element Vector{FreeModElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y*e[1]
 x*e[1] + y*e[2]
 z*e[2]

julia> a = hom(F, G, V)
F -> G
e[1] -> y*e[1]
e[2] -> x*e[1] + y*e[2]
e[3] -> z*e[2]
Graded module homomorphism of degree [1]

julia> degree(a)
[1]
```
"""
function degree(f::FreeModuleHom)
  if isdefined(f, :d)
    return f.d
  end
  T1 = domain(f)
  T2 = codomain(f)
  if !is_graded(T1) || !is_graded(T2)
    error("Both domain and codomain must be graded.")
  end
  domain_degrees = degrees(T1)
  df = nothing
  for i in 1:length(domain_degrees)
    image_vector = f(T1[i])
    if isempty(coordinates(image_vector))
      continue
    end
    current_df = degree(image_vector) - domain_degrees[i]
    if df === nothing
      df = current_df
    elseif df != current_df
      error("The homomorphism is not graded")
    end
  end
  if df === nothing
    R = base_ring(T1)
    G = grading_group(R)
    return G[0]
  end
  return df
end

@doc raw"""
    is_graded(a::FreeModuleHom)

Return `true` if `a` is graded, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(R, 3)
Graded free module R^3([0]) of rank 3 over R

julia> G = graded_free_module(R, 2)
Graded free module R^2([0]) of rank 2 over R

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
3-element Vector{FreeModElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y*e[1]
 x*e[1] + y*e[2]
 z*e[2]

julia> a = hom(F, G, V)
F -> G
e[1] -> y*e[1]
e[2] -> x*e[1] + y*e[2]
e[3] -> z*e[2]
Graded module homomorphism of degree [1]

julia> is_graded(a)
true
```
"""
function is_graded(f::FreeModuleHom)
  return isdefined(f, :d)
end

@doc raw"""
    grading_group(a::FreeModuleHom)

If `a` is graded, return the grading group of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(R, 3)
Graded free module R^3([0]) of rank 3 over R

julia> G = graded_free_module(R, 2)
Graded free module R^2([0]) of rank 2 over R

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
3-element Vector{FreeModElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y*e[1]
 x*e[1] + y*e[2]
 z*e[2]

julia> a = hom(F, G, V)
F -> G
e[1] -> y*e[1]
e[2] -> x*e[1] + y*e[2]
e[3] -> z*e[2]
Graded module homomorphism of degree [1]

julia> is_graded(a)
true

julia> grading_group(a)
GrpAb: Z
```
"""
function grading_group(f::FreeModuleHom)
  return grading_group(base_ring(domain(f)))
end

@doc raw"""
    is_homogeneous(a::FreeModuleHom)

Return `true` if `a` is homogeneous, `false` otherwise

Here, if `G` is the grading group of `a`, `a` is homogeneous if `a`
is graded of degree `zero(G)`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(R, 3)
Graded free module R^3([0]) of rank 3 over R

julia> G = graded_free_module(R, 2)
Graded free module R^2([0]) of rank 2 over R

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
3-element Vector{FreeModElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y*e[1]
 x*e[1] + y*e[2]
 z*e[2]

julia> a = hom(F, G, V)
F -> G
e[1] -> y*e[1]
e[2] -> x*e[1] + y*e[2]
e[3] -> z*e[2]
Graded module homomorphism of degree [1]

julia> is_homogeneous(a)
false
```
"""
function is_homogeneous(f::FreeModuleHom)
  A = grading_group(f)
  return isdefined(f, :d) && degree(f)==A[0]
end

###############################################################################
# Graded submodules
###############################################################################

function is_graded(M::SubModuleOfFreeModule)
  is_graded(M.F) && all(is_homogeneous, M.gens)
end


function degrees_of_generators(M::SubModuleOfFreeModule{T}) where T
  return map(gen -> degree(gen), gens(M))
end

###############################################################################
# Graded subquotient constructors
###############################################################################

# mostly automatic, just needed for matrices

function graded_cokernel(A::MatElem)
  return cokernel(graded_map(A))
end

function graded_cokernel(F::FreeMod{R}, A::MatElem{R}) where R
  @assert is_graded(F)
  cokernel(graded_map(F,A))
end

function graded_image(F::FreeMod{R}, A::MatElem{R}) where R
  @assert is_graded(F)
  image(graded_map(F,A))[1]
end

function graded_image(A::MatElem)
  return image(graded_map(A))[1]
end

###############################################################################
# Graded subquotients
###############################################################################

@doc raw"""
    grading_group(M::SubquoModule)

Return the grading group of `base_ring(M)`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F1 = graded_free_module(Rg, [2,2,2]);

julia> F2 = graded_free_module(Rg, [2]);

julia> G = graded_free_module(Rg, [1,1]);

julia> V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]];

julia> V2 = [z*G[2]+y*G[1]];

julia> a1 = hom(F1, G, V1);

julia> a2 = hom(F2, G, V2);

julia> M = subquotient(a1,a2);

julia> grading_group(M)
GrpAb: Z
```
"""
function grading_group(M::SubquoModule)
  return grading_group(base_ring(M))
end

@doc raw"""
    degrees_of_generators(M::SubquoModule)

Return the degrees of the generators of `M`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F1 = graded_free_module(Rg, [2,2,2]);

julia> F2 = graded_free_module(Rg, [2]);

julia> G = graded_free_module(Rg, [1,1]);

julia> V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]];

julia> V2 = [z*G[2]+y*G[1]];

julia> a1 = hom(F1, G, V1);

julia> a2 = hom(F2, G, V2);

julia> M = subquotient(a1,a2);

julia> degrees_of_generators(M)
3-element Vector{GrpAbFinGenElem}:
 [2]
 [2]
 [2]

julia> gens(M)
3-element Vector{SubquoModuleElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y*e[1]
 (x + y)*e[1] + y*e[2]
 z*e[2]
```
"""
function degrees_of_generators(M::SubquoModule{T}) where T
  isempty(gens(M)) ? GrpAbFinGenElem[] : map(gen -> degree(repres(gen)), gens(M))
end

###############################################################################
# Graded subquotient elements
###############################################################################

 @doc raw"""
    is_homogeneous(m::SubquoModuleElem)

Return  `true` if `m` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F1 = graded_free_module(Rg, [2,2,2]);

julia> F2 = graded_free_module(Rg, [2]);

julia> G = graded_free_module(Rg, [1,1]);

julia> V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]];

julia> V2 = [z*G[2]+y*G[1]];

julia> a1 = hom(F1, G, V1);

julia> a2 = hom(F2, G, V2);

julia> M = subquotient(a1,a2);

julia> m1 = x*M[1]+y*M[2]+z*M[3]
(2*x*y + y^2)*e[1] + (y^2 + z^2)*e[2]

julia> is_homogeneous(m1)
true

julia> is_homogeneous(zero(M))
true

julia> m2 = M[1]+x*M[2]
(x^2 + x*y + y)*e[1] + x*y*e[2]

julia> is_homogeneous(m2)
false

julia> m3 = x*M[1]+M[2]+x*M[3]
(x*y + x + y)*e[1] + (x*z + y)*e[2]

julia> is_homogeneous(m3)
true

julia> simplify(m3)
x*e[1] + (y - z)*e[2]
```
"""
function is_homogeneous(el::SubquoModuleElem)
  if iszero(el.coeffs)
      return is_homogeneous(repres(el))
  else
      degree = determine_degree_from_SR(el.coeffs, degrees_of_generators(parent(el)))
      if degree === nothing
          reduced_el = simplify(el)
          degree_reduced = determine_degree_from_SR(reduced_el.coeffs, degrees_of_generators(parent(reduced_el)))
          return degree_reduced !== nothing
      else
          return true
      end
  end
end

@doc raw"""
    degree(m::SubquoModuleElem)

Given a homogeneous element `m` of a graded subquotient, return the degree of `m`.

    degree(::Type{Vector{Int}}, m::SubquoModuleElem)

Given a homogeneous element `m` of a $\mathbb Z^m$-graded subquotient, return the degree of `m`, converted to a vector of integer numbers.

    degree(::Type{Int}, m::SubquoModuleElem)

Given a homogeneous element `m` of a $\mathbb Z$-graded subquotient, return the degree of `m`, converted to an integer number.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F1 = graded_free_module(Rg, [2,2,2]);

julia> F2 = graded_free_module(Rg, [2]);

julia> G = graded_free_module(Rg, [1,1]);

julia> V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]];

julia> V2 = [z*G[2]+y*G[1]];

julia> a1 = hom(F1, G, V1);

julia> a2 = hom(F2, G, V2);

julia> M = subquotient(a1,a2);

julia> m = x*y*z*M[1]
x*y^2*z*e[1]

julia> degree(m)
[5]

julia> degree(Int, m)
5

julia> m3 = x*M[1]+M[2]+x*M[3]
(x*y + x + y)*e[1] + (x*z + y)*e[2]

julia> degree(m3)
[2]
```
"""
function degree(el::SubquoModuleElem)
  if !iszero(el.coeffs)
      result = determine_degree_from_SR(el.coeffs, degrees_of_generators(parent(el)))
      if result === nothing
          reduced_el = simplify(el)
          result_reduced = determine_degree_from_SR(reduced_el.coeffs, degrees_of_generators(parent(reduced_el)))
          @assert result_reduced !== nothing "The specified element is not homogeneous."
          return result_reduced
      else
          return result
      end
  else
      return degree(repres(el))
  end
end

function degree(::Type{Vector{Int}}, el::SubquoModuleElem)
  @assert is_zm_graded(parent(el))
  d = degree(el)
  return Int[d[i] for i=1:ngens(parent(d))]
end

function degree(::Type{Int}, el::SubquoModuleElem)
  @assert is_z_graded(parent(el))
  return Int(degree(el)[1])
end


###############################################################################
# Graded subquotient homomorphisms functions
###############################################################################

function set_grading(f::SubQuoHom)
  if !is_graded(domain(f)) || !is_graded(codomain(f))
    return(f)
  end
  f.d = degree(f)
  return f
end

@doc raw"""
    degree(a::SubQuoHom)

If `a` is graded, return the degree of `a`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
M -> M
x*e[1] -> x*y^2*e[1]
y*e[1] -> x^2*y*e[1]
Graded module homomorphism of degree [2]

julia> degree(a)
[2]
```
"""
function degree(f::SubQuoHom)
  if isdefined(f, :d)
    return f.d
  end
  T1 = domain(f)
  T2 = codomain(f)
  if !is_graded(T1) || !is_graded(T2)
    error("Both domain and codomain must be graded.")
  end
  if iszero(T1)
    R = base_ring(T1)
    G = grading_group(R)
    return G[0]
  end
  domain_degrees = degrees_of_generators(T1)
  df = nothing
  for i in 1:length(domain_degrees)
    image_vector = f(T1[i])
    if isempty(coordinates(image_vector))
      continue
    end
    current_df = degree(image_vector) - domain_degrees[i]
    if df === nothing
      df = current_df
    elseif df != current_df
      error("The homomorphism is not graded.")
    end
  end
  if df === nothing
    R = base_ring(T1)
    G = grading_group(R)
    return G[0]
  end
  return df
end

@doc raw"""
    is_graded(a::SubQuoHom)

Return `true` if `a` is graded, `false` otherwise.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
M -> M
x*e[1] -> x*y^2*e[1]
y*e[1] -> x^2*y*e[1]
Graded module homomorphism of degree [2]

julia> is_graded(a)
true
```
"""
function is_graded(f::SubQuoHom)
  return isdefined(f, :d)
end

@doc raw"""
    grading_group(a::SubQuoHom)

If `a` is graded, return the grading group of `a`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
M -> M
x*e[1] -> x*y^2*e[1]
y*e[1] -> x^2*y*e[1]
Graded module homomorphism of degree [2]

julia> grading_group(a)
GrpAb: Z
```
"""
function grading_group(f::SubQuoHom)
  return grading_group(base_ring(domain(f)))
end

@doc raw"""
    is_homogeneous(a::SubQuoHom)

Return `true` if `a` is homogeneous, `false` otherwise

Here, if `G` is the grading group of `a`, `a` is homogeneous if `a`
is graded of degree `zero(G)`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
M -> M
x*e[1] -> x*y^2*e[1]
y*e[1] -> x^2*y*e[1]
Graded module homomorphism of degree [2]

julia> is_homogeneous(a)
false
```
"""
function is_homogeneous(f::SubQuoHom)
  A = grading_group(f)
  return isdefined(f, :d) && degree(f)==A[0]
end


###############################################################################
# Graded free resolutions
###############################################################################

function is_graded(resolution::FreeResolution{T}) where T
  C = resolution.C
 return all(is_graded(C[i]) for i in reverse(Hecke.range(C))) && all(is_graded(map(C, i)) for i in reverse(Hecke.map_range(C)))
end

###############################################################################
# Betti table
###############################################################################

@doc raw"""
    betti_table(F::FreeResolution)

Given a $\mathbb Z$-graded free resolution `F`, return the graded Betti numbers 
of `F` in form of a Betti table.

Alternatively, use `betti`.

# Examples
```julia
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [x*z, y*z, x*w^2, y*w^2])
ideal(x*z, y*z, w^2*x, w^2*y)

julia> A, _= quo(R, I)
(Quotient of multivariate polynomial ring by ideal with 4 generators, Map from
R to A defined by a julia-function with inverse)

julia> FA  = free_resolution(A)
Free resolution of A
R^1 <---- R^4 <---- R^4 <---- R^1 <---- 0
0         1         2         3         4

julia> betti_table(FA)
       0  1  2  3
------------------
0    : 1  -  -  -
1    : -  2  1  -
2    : -  2  3  1
------------------
total: 1  4  4  1

julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x, y, x+y]);

julia> M = quotient_ring_as_module(I);

julia> FM = free_resolution(M, algorithm = :nres);

julia> betti_table(FM)
       0  1  2
---------------
-1   : -  -  1
0    : 1  3  1
---------------
total: 1  3  2
```
"""
function betti_table(F::FreeResolution; project::Union{GrpAbFinGenElem, Nothing} = nothing, reverse_direction::Bool=false)
  generator_count = Dict{Tuple{Int, Any}, Int}()
  C = F.C
  rng = Hecke.map_range(C)
  n = first(rng)
  for i in 0:n
    module_degrees = F[i].d
    module_degrees === nothing && error("One of the modules in the graded free resolution is not graded.")
    for degree in module_degrees
      idx = (i, degree)
      generator_count[idx] = get(generator_count, idx, 0) + 1
    end
  end
  return BettiTable(generator_count, project = project, reverse_direction = reverse_direction)
end

function betti(b::FreeResolution; project::Union{GrpAbFinGenElem, Nothing} = nothing, reverse_direction::Bool = false)
  return betti_table(b; project, reverse_direction)
end

function as_dictionary(b::BettiTable)
  return b.B
end

function reverse_direction!(b::BettiTable)
  b.reverse_direction = !b.reverse_direction
  return
end

function induce_shift(B::Dict{Tuple{Int, Any}, Int})
  A = parent(first(keys(B))[2])
  new_B = Dict{Tuple{Int, Any}, Int}()
  for ((i, key), value) in B
    new_key = (i, key-i*A[1])
    new_B[new_key] = value
  end
  return new_B
end

function Base.show(io::IO, b::BettiTable)
  T = induce_shift(b.B)
  x = collect(keys(T))
  if isempty(x)
    println(io, "Empty table")
    return
  end
  step, min, maxv = b.reverse_direction ? (-1, maximum(first, x), minimum(first, x)) : (1, minimum(first, x), maximum(first, x))
  column_widths = Dict()
  for j in min:step:maxv
    sum_col = sum(getindex(T, x[m]) for m in 1:length(x) if x[m][1] == j)
    col_width_from_sum = ndigits(abs(sum_col))
    col_width_from_header = ndigits(abs(j)) + (j < 0 ? 1 : 0)
    column_widths[j] = max(col_width_from_sum, col_width_from_header) + 2
  end

  if b.project === nothing
    for i in 1:ngens(parent(x[1][2]))
      ngens(parent(x[1][2])) > 1 && println(io, "Betti Table for component ", i)
      L = sort(unique(collect(x[k][2][i] for k in 1:length(x))))
      mi = minimum(L)
      mx = maximum(L)
      initial_padding = max(ndigits(mi) + mi < 0 ? 0 : 1, 7)
      print(io, " "^initial_padding)
      total_space_count = initial_padding
      for j in min:step:maxv
        adjustment = j < 0 ? 1 : 0
        space_count = max(0, column_widths[j] - ndigits(j) - adjustment)
        print(io, j)
        if j != maxv
          print(io, " "^space_count)
        end
        total_space_count = total_space_count + space_count + ndigits(j) + adjustment
      end
      total_space_count = total_space_count - 1
      print(io, "\n")
      print(io, "-"^total_space_count)
      print(io, "\n")
      for j in mi:mx
        adjustment = j < 0 ? 1 : 0
        print(io, j, " "^(5 - ndigits(j) - adjustment))
        print(io, ":")
        for h in min:step:maxv
          sum_current = sum([getindex(T, x[k]) for k in 1:length(x) if x[k][1] == h && x[k][2][i] == j])
          print(io, " ", sum_current == 0 ? "-" : sum_current)
          if h != maxv
            print(io, " "^(column_widths[h] - ndigits(sum_current) - 1))
          end
        end
        print(io,"\n")
      end
      print(io, "-" ^ total_space_count)
      print(io, "\n", "total:")
      for i_total in min:step:maxv
        sum_row = sum(getindex(T, x[j]) for j in 1:length(x) if x[j][1] == i_total)
        print(io, " ", sum_row)
        if i_total != maxv
          print(io, " "^(column_widths[i_total] - ndigits(sum_row) - 1))
        end
      end
      print(io, "\n")
    end
  else
    parent(b.project) == parent(x[1][2]) || error("projection vector has wrong type")
    print(io, "Betti Table for scalar product of grading with ", b.project.coeff, "\n")
    print(io, "  ")
    L = Vector{fmpz}(undef,0)
    for i in 1:length(x)
        temp_sum = (b.project.coeff * transpose(x[i][2].coeff))[1]
        Base.push!(L, temp_sum)
    end
    L1 = sort(unique(L))
    s2 = ndigits(maximum(L1))
    spaces = maximum([s1, s2, s3])
    print(io, " " ^ (s2 + (5 - s2)))
    for j in min:step:max
        print(io, j, " " ^ (spaces - ndigits(j) + 1))
    end
    print(io, "\n")
    for k in 1:length(L1)
        print(io, L1[k], " " ^ (s2 - ndigits(L1[k]) + (5 - s2)), ": ")
        for h in min:step:max
            partial_sum = 0
            for i in 1:length(x)
                current_sum = (b.project.coeff * transpose(x[i][2].coeff))[1]
                if current_sum == L1[k] && x[i][1] == h
                    partial_sum += getindex(T, x[i])
                end
            end
            if partial_sum == 0
                print(io, "-", " " ^ spaces)
            else
                print(io, partial_sum, " " ^ (spaces - ndigits(partial_sum) + 1))
            end
        end
        print(io, "\n")
    end
    divider_width = 6 + (b.reverse_direction ? (min - max) + 1 : (max - min) + 1) * (spaces + 1)
    print(io, "-" ^ divider_width)
    print(io, "\n", "total: ")
    for i in min:step:max
        total_sum = 0
        for j in 1:length(x)
            if x[j][1] == i
                total_sum += getindex(T, x[j])
            end
        end
        print(io, total_sum, " " ^ (spaces - ndigits(total_sum) + 1))
    end
  end
end




###############################################################################
# data structure for table as in
# https://www.singular.uni-kl.de/Manual/4-3-1/sing_1827.htm#SEC1908
###############################################################################

mutable struct sheafCohTable
  twist_range::UnitRange{Int}
  values::Matrix{Int}
end

function Base.getindex(st::sheafCohTable, ind...)
  row_ind = ind[1] + 1
  col_ind = ind[2] - first(st.twist_range) + 1
  return st.values[row_ind, col_ind]
end

function Base.show(io::IO, table::sheafCohTable)
  chi = [any(v -> v == -1, col) ? -1 : sum(col) for col in eachcol(table.values)]
  # pad every value in the table to this length
  val_space_length = max(maximum(_ndigits, table.values), maximum(_ndigits, chi))
  nrows = size(table.values, 1)

  # rows to print
  print_rows = [[_shcoh_string_rep(v, val_space_length) for v in row]
                for row in eachrow(table.values)]
  chi_print =  [_shcoh_string_rep(v, val_space_length) for v in chi]

  # row labels
  row_label_length = max(_ndigits(nrows - 1), 3) + 3
  for i in 1:nrows
    pushfirst!(print_rows[i], rpad("$(i-1): ", row_label_length, " "))
  end
  pushfirst!(chi_print, rpad("chi: ", row_label_length, " "))  

  # header
  header = [lpad(v, val_space_length, " ") for v in table.twist_range]
  pushfirst!(header, rpad("twist:", row_label_length, " "))

  println(io, header...)
  size_row = sum(length, first(print_rows))
  println(io, repeat("-", size_row))
  for rw in print_rows
    println(io, rw...)
  end
  println(io, repeat("-", size_row))
  print(io, chi_print...)
end

@doc raw"""
    sheaf_cohomology_bgg(M::ModuleFP{T}, l::Int, h::Int) where {T <: MPolyDecRingElem}

Compute the cohomology of twists of of the coherent sheaf on projective
space associated to `M`. The range of twists is between `l` and `h`.
In the displayed result, '-' refers to a zero enty and '*' refers to a
negative entry (= dimension not yet determined). To determine all values
in the desired range between `l` and `h` use `sheafCoh_BGG_regul(M, l-ngens(base_ring(M)), h+ngens(base_ring(M)))`.
The values of the returned table can be accessed by indexing it
with a cohomological index and a value between `l` and `h` as shown
in the example below.

```jldoctest
julia> R, x = polynomial_ring(QQ, "x" => 1:4);

julia> S, _= grade(R);

julia> I = ideal(S, gens(S))
ideal(x[1], x[2], x[3], x[4])

julia> FI = free_resolution(I)
Free resolution of I
S^4 <---- S^6 <---- S^4 <---- S^1 <---- 0
0         1         2         3         4

julia> M = cokernel(map(FI, 2));

julia> tbl = sheaf_cohomology_bgg(M, -6, 2)
twist:  -6  -5  -4  -3  -2  -1   0   1   2
------------------------------------------
0:      70  36  15   4   -   -   -   -   *
1:       *   -   -   -   -   -   -   -   -
2:       *   *   -   -   -   -   1   -   -
3:       *   *   *   -   -   -   -   -   6
------------------------------------------
chi:     *   *   *   4   -   -   1   -   *

julia> tbl[0, -6]
70

julia> tbl[2, 0]
1

julia> R, x = polynomial_ring(QQ, "x" => 1:5);

julia> R, x = grade(R);

julia> F = graded_free_module(R, 1);

julia> sheaf_cohomology_bgg(F, -7, 2)
twist:  -7  -6  -5  -4  -3  -2  -1   0   1   2
----------------------------------------------
0:      15   5   1   -   -   -   *   *   *   *
1:       *   -   -   -   -   -   -   *   *   *
2:       *   *   -   -   -   -   -   -   *   *
3:       *   *   *   -   -   -   -   -   -   *
4:       *   *   *   *   -   -   -   1   5  15
----------------------------------------------
chi:     *   *   *   *   -   -   *   *   *   *
```
"""
function sheaf_cohomology_bgg(M::ModuleFP{T},
                              l::Int,
                              h::Int) where {T <: MPolyDecRingElem}


  free_mod = ambient_free_module(M)
  @assert is_graded(M) "Module must be graded."
  @assert is_standard_graded(free_mod) "Only supported for the standard grading of polynomials."

  reg=Int(cm_regularity(M))
  # get a cokernel presentation of M
  p = presentation(M)
  cokern_repr = image(map(p, 1))[1]
  cokern_gens = ambient_representatives_generators(cokern_repr)
  if isempty(cokern_gens)
    cokern_gens = [zero(ambient_free_module(cokern_repr))]
  end
  sing_mod = singular_generators(ModuleGens(cokern_gens))

  weights = [Int(d[1]) for d in degrees_of_generators(free_mod)]

  values = Singular.LibSheafcoh.sheafCohBGGregul_w(sing_mod,
                                                   l, h, reg,
                                                   weights)
  return sheafCohTable(l:h, values)
end

# helper functions
function _shcoh_string_rep(val::Int, padlength::Int)
  iszero(val) && return lpad("-", padlength, " ")
  val == -1 && return lpad("*", padlength, " ")
  return lpad(val, padlength, " ")
end

function _ndigits(val::Int)
  iszero(val) && return 3
  val == -1 && return 3
  return Int(floor(log10(val))) + 3
end

##################################
### Tests on graded modules
##################################

function is_graded(M::FreeMod)
  return isa(M.d, Vector{GrpAbFinGenElem})
end

function is_graded(M::SubquoModule)
  if isdefined(M, :quo)
    return is_graded(M.sub) && is_graded(M.quo) && is_graded(M.sum)
  else
    return is_graded(M.sub)
  end
end

function is_standard_graded(M::FreeMod)
  return is_graded(M) && is_standard_graded(base_ring(M))
end

function is_standard_graded(M::SubquoModule)
  return  is_graded(M) && is_standard_graded(base_ring(M))
end

function is_z_graded(M::FreeMod)
  return  is_graded(M) && is_z_graded(base_ring(M))
end

function is_z_graded(M::SubquoModule)
  return  is_graded(M) && is_z_graded(base_ring(M))
end

function is_zm_graded(M::FreeMod)
  return  is_graded(M) && is_zm_graded(base_ring(M))
end

function is_zm_graded(M::SubquoModule)
  return  is_graded(M) && is_zm_graded(base_ring(M))
end


###############################################################################
# FreeMod_dec constructors
###############################################################################

@doc raw"""
    FreeMod_dec(R::CRing_dec, n::Int, name::VarName = :e; cached::Bool = false) 

Construct a decorated (graded or filtered) free module over the ring `R` with rank `n`
with the standard degrees, that is the standard unit vectors have degree 0.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for 
the standard unit vectors.
"""
function FreeMod_dec(R::CRing_dec, n::Int, name::VarName = :e; cached::Bool = false) 
  return FreeMod_dec{elem_type(R)}(R, [Symbol("$name[$i]") for i=1:n], [decoration(R)[0] for i=1:n])
end

@doc raw"""
    free_module_dec(R::CRing_dec, n::Int, name::VarName = :e; cached::Bool = false)

Create the decorated free module $R^n$ equipped with its basis of standard unit vectors
and standard degrees, that is the standard unit vectors have degree 0.

The string `name` specifies how the basis vectors are printed. 

# Examples
```jldoctest
julia> R, (x,y) = graded_polynomial_ring(QQ, ["x", "y"])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> free_module_dec(R,3)
Decorated free module of rank 3 over RR^3([0])

```
"""
free_module_dec(R::CRing_dec, n::Int, name::VarName = :e; cached::Bool = false) = FreeMod_dec(R, n, name, cached = cached)


@doc raw"""
    FreeMod_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::VarName = :e; cached::Bool = false) 

Construct a decorated (graded or filtered) free module over the ring `R` 
with rank `n` where `n` is the length of `d`. `d` is the vector of degrees for the 
components, i.e. `d[i]` is the degree of `e[i]` where `e[i]` is the `i`th standard unit
vector of the free module.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for 
the standard unit vectors.
"""
function FreeMod_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::VarName = :e; cached::Bool = false) 
  return FreeMod_dec{elem_type(R)}(R, [Symbol("$name[$i]") for i=1:length(d)],d)
end

@doc raw"""
    free_module_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::VarName = :e; cached::Bool = false)

Create the decorated free module $R^n$ (`n` is the length of `d`)
equipped with its basis of standard unit vectors where the 
i-th standard unit vector has degree `d[i]`.

The string `name` specifies how the basis vectors are printed. 
"""
free_module_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::VarName = :e; cached::Bool = false) = FreeMod_dec(R, d, name, cached = cached)


function FreeMod_dec(F::FreeMod, d::Vector{GrpAbFinGenElem})
  return FreeMod_dec{elem_type(base_ring(F))}(F, d)
end


function AbstractAlgebra.extra_name(F::FreeMod_dec)
  t = get_attribute(F, :twist)
  if t !== nothing
    n = get_attribute(t[1], :name)
    if n !== nothing
      return "$n($(t[2]))"
    end
  end
  if length(Set(F.d)) == 1
    n = get_attribute(forget_decoration(F).R, :name)
    if n !== nothing
      return "$n^$(ngens(F))($(-F.d[1]))"
    end
  end
  return nothing
end

function show(io::IO, F::FreeMod_dec)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Decorated free module of rank $(rank(F)) over ")
  print(IOContext(io, :compact =>true), base_ring(F))

  i = 1
  while i < dim(F)
    d = F.d[i]
    j = 1
    while i+j <= dim(F) && d == F.d[i+j]
      j += 1
    end
    print(IOContext(io, :compact => true), base_ring(F), "^$j")
    print(IOContext(io, :compact => true), "(", -d, ")")
    if i+j < dim(F)
      print(io, " + ")
    end
    i += j
  end
end

function forget_decoration(F::FreeMod_dec)
  return F.F
end

@doc raw"""
    base_ring(F::FreeMod_dec)

Return the underlying ring of `F`.
"""
base_ring(F::FreeMod_dec) = forget_decoration(F).R

@doc raw"""
    rank(F::FreeMod_dec)

Return the rank of `F`.
"""
rank(F::FreeMod_dec) = rank(forget_decoration(F))

@doc raw"""
    decoration(F::FreeMod_dec)

Return the vector of degrees of the standard unit vectors.
"""
decoration(F::FreeMod_dec) = F.d
decoration(R::MPolyDecRing) = R.D

@doc raw"""
    is_graded(F::FreeMod_dec)

Check if `F` is graded.
"""
is_graded(F::FreeMod_dec) = is_graded(base_ring(F))

@doc raw"""
    is_filtered(F::FreeMod_dec)

Check if `F` is filtered.
"""
is_filtered(F::FreeMod_dec) = is_filtered(base_ring(F))

is_decorated(F::FreeMod_dec) = true

@doc raw"""
    ==(F::FreeMod_dec, G::FreeMod_dec)

Return  `true` if `F` and `G` are equal, `false` otherwise.

Here, `F` and `G` are equal iff their base rings, ranks, decorations 
and names for printing the basis elements are equal.
"""
function Base.:(==)(F::FreeMod_dec, G::FreeMod_dec)
  return forget_decoration(F) == forget_decoration(G) && F.d == G.d
end

function Base.hash(F::FreeMod_dec, h::UInt)
  b = 0x13d6e1b453cb661a % UInt
  return xor(hash(forget_decoration(F), hash(F.d, h)), b)
end

###############################################################################
# FreeModElem_dec constructors
###############################################################################

@doc raw"""
    FreeModElem_dec(c::SRow{T}, parent::FreeMod_dec{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
FreeModElem_dec(c::SRow{T}, parent::FreeMod_dec{T}) where T = FreeModElem_dec{T}(c, parent)

@doc raw"""
    FreeModElem_dec(c::Vector{T}, parent::FreeMod_dec{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
function FreeModElem_dec(c::Vector{T}, parent::FreeMod_dec{T}) where T
  @assert length(c) == rank(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:rank(parent)), c)
  return FreeModElem_dec{T}(sparse_coords,parent)
end

#@doc raw"""
#    (F::FreeMod_dec{T})(c::SRow{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#"""
function (F::FreeMod_dec{T})(c::SRow{T}) where T
  return FreeModElem_dec(c, F)
end

#@doc raw"""
#    (F::FreeMod_dec{T})(c::Vector{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#"""
function (F::FreeMod_dec{T})(c::Vector{T}) where T 
  return FreeModElem_dec(c, F)
end

@doc raw"""
    (F::FreeMod_dec)()

Return the zero element of `F`.
"""
function (F::FreeMod_dec)()
  return FreeModElem_dec(sparse_row(base_ring(F)), F)
end

@doc raw"""
    FreeModElem_dec(v::FreeModElem{T}, parent::FreeMod_dec{T}) where T <: CRingElem_dec

Lift `v` to the decorated module `parent`.
"""
function FreeModElem_dec(v::FreeModElem{T}, p::FreeMod_dec{T}) where T <: CRingElem_dec
  @assert forget_decoration(p) === parent(v)
  return FreeModElem_dec(coordinates(v), p)
end


elem_type(::Type{FreeMod_dec{T}}) where {T} = FreeModElem_dec{T}
parent_type(::Type{FreeModElem_dec{T}}) where {T} = FreeMod_dec{T}
elem_type(::FreeMod_dec{T}) where {T} = FreeModElem_dec{T}
parent_type(::FreeModElem_dec{T}) where {T} = FreeMod_dec{T}

@doc raw"""
"""
function forget_decoration(v::FreeModElem_dec)
  return FreeModElem(coordinates(v),forget_decoration(parent(v)))
end


@doc raw"""
    generator_symbols(F::FreeMod_dec)

Return the list of symbols of the standard unit vectors.
"""
function generator_symbols(F::FreeMod_dec)
  return generator_symbols(forget_decoration(F))
end
@enable_all_show_via_expressify FreeModElem_dec


@doc raw"""
    degree_homogeneous_helper(u::FreeModElem_dec)

Compute the degree and homogeneity of `u` (simultaneously).
A tuple is returned: The first entry is the degree (or nothing if
the element has no degree); the second entry is true if `u` is 
homogeneous and false otherwise.
"""
function degree_homogeneous_helper(u::FreeModElem_dec)
  if iszero(u)
    return nothing, true
  end
  first = true
  homogeneous = true #only needed in filtered case
  F = parent(u)
  W = base_ring(F)
  ww = W.D[0]
  local w
  for (p,v) in coordinates(u)
    if !is_homogeneous(v)
      if is_graded(W)
        return nothing,false
      else
        homogeneous = false
      end
    end
    w = degree(v)+F.d[p]
    if first
      ww = w
      first = false
    elseif is_graded(W)
      if ww != w
        return nothing, false
      end
    else
      if ww != w
        homogeneous = false
      end
      if W.lt(ww, w) 
        ww = w
      end
    end
  end
  return ww, homogeneous
end

@doc raw"""
    degree(a::FreeModElem_dec)

Return the degree of `a`. If `a` has no degree an error is thrown.
"""
function degree(a::FreeModElem_dec)
  d,_ = degree_homogeneous_helper(a)
  d === nothing ? error("elem has no degree") : return d
end

@doc raw"""
    homogeneous_components(a::FreeModElem_dec)

Return the homogeneous components of `a` in a dictionary.
The keys are those group elements for which `a` has a component
having this element as its degree.
"""
function homogeneous_components(a::FreeModElem_dec)
  res = Dict{GrpAbFinGenElem, FreeModElem_dec}()
  F = parent(a)
  for (p,v) in coordinates(a)
    c = homogeneous_components(v)
    for (pp, vv) in c
      w = pp + F.d[p]
      if haskey(res, w)
        res[w] += vv*gen(F, p)
      else
        res[w] = vv*gen(F, p)
      end
    end
  end
  return res
end

@doc raw"""
    homogeneous_component(a::FreeModElem_dec, g::GrpAbFinGenElem)

Return the homogeneous component of `a` which has degree `g`.
"""
function homogeneous_component(a::FreeModElem_dec, g::GrpAbFinGenElem)
  F = parent(a)
  x = zero(F)
  for (p,v) in coordinates(a)
    x += homogeneous_component(v, g-F.d[p])*gen(F, p)
  end
  return x
end

@doc raw"""
    is_homogeneous(a::FreeModElem_dec)

Check if `a` is homogeneous.
"""
function is_homogeneous(a::FreeModElem_dec)
  return degree_homogeneous_helper(a)[2]
end

# Weight vector or function?
# Should we already grade ModuleGens?
# Should it be possible to construct ungraded SubquoModule with graded elements? (I.e. should the constructors
# accept AbstractFreeMod and AbstractFreeModElem instead of FreeMod and FreeModElem?)
# proceed with FreeModHom_dec?



@doc raw"""
    tensor_product(G::FreeMod_dec...; task::Symbol = :none)

Given decorated free modules $G_i$ compute the decorated tensor product 
$G_1\otimes \cdots \otimes G_n$.
If `task` is set to ":map", a map $\phi$ is returned that
maps tuples in $G_1 \times \cdots \times G_n$ to pure tensors
$g_1 \otimes \cdots \otimes g_n$. The map admits a preimage as well.
"""
function tensor_product(G::FreeMod_dec...; task::Symbol = :none)
  undecorated_tensor_product, tuple_to_pure = tensor_product(map(forget_decoration, G)...; task=:map)
  pure_to_tuple = inv(tuple_to_pure)
  d = [sum(map(degree, [FreeModElem_dec(elem,parent) for (elem,parent) in zip(pure_to_tuple(v),G)])) 
                                                    for v in gens(undecorated_tensor_product)]
  F = FreeMod_dec(undecorated_tensor_product, d)

  function pure(T::Tuple)
    return FreeModElem_dec(tuple_to_pure(map(forget_decoration, T)), F)
  end

  function inv_pure(e::FreeModElem_dec)
    a = pure_to_tuple(forget_decoration(e))
    return Tuple(FreeModElem_dec(elem,parent) for (elem,parent) in zip(a,G))
  end

  set_attribute!(F, :tensor_pure_function => pure, :tensor_generator_decompose_function => inv_pure)

  if task == :none
    return F
  end
  return F, MapFromFunc(Hecke.TupleParent(Tuple([g[0] for g = G])), F, pure, inv_pure)
end


###############################################################################
# FreeModuleHom_dec constructors
###############################################################################

FreeModuleHom_dec(F::FreeMod_dec{T}, G::ModuleFP_dec, a::Vector) where {T} = FreeModuleHom_dec{T}(F, G, a)

FreeModuleHom_dec(F::FreeMod_dec{T}, G::ModuleFP_dec, mat::MatElem{T}) where {T} = FreeModuleHom{T}(F, G, mat)

function forget_decoration_on_morphism(f::FreeModuleHom_dec)
  return f.f
end

function forget_decoration(f::FreeModuleHom_dec)
  F = forget_decoration(domain(f))
  G = forget_decoration(codomain(f))
  return hom(F, G, [forget_decoration(f(v)) for v in gens(domain(f))])
end

function matrix(a::FreeModuleHom_dec)
  return matrix(forget_decoration_on_morphism(a))
end

(h::FreeModuleHom_dec)(a::FreeModElem_dec) = image(h, a)

hom(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, V::Vector{<:FreeModElem_dec}) where T = FreeModuleHom_dec(F, G, V) 
hom(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, A::MatElem{T}) where T = FreeModuleHom_dec(F, G, A)


function hom(F::FreeMod_dec, G::FreeMod_dec)
  undecorated_hom, elem_to_hom = hom(forget_decoration(F), forget_decoration(G))
  d = [y-x for x in decoration(F) for y in decoration(G)]
  GH = FreeMod_dec(undecorated_hom, d)
  X = Hecke.MapParent(F, G, "homomorphisms")

  function im(v::FreeModElem_dec)
    return hom(F, G, [FreeModElem_dec(elem_to_hom(forget_decoration(v))(forget_decoration(u)),G) for u in gens(F)])
  end

  function pre(f::FreeModuleHom_dec)
    undecorated_v = inv(elem_to_hom)(forget_decoration(f))
    return FreeModElem_dec(undecorated_v, GH)
  end

  to_hom_map = MapFromFunc(GH, X, im, pre)
  set_attribute!(GH, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map)
  return GH, to_hom_map
end

                                                                                                      
########################################################################
# Minimal Betti tables
########################################################################

# TODO: Are the signatures sufficient to assure that the modules are graded?

@doc raw"""
    minimal_betti_table(M::ModuleFP{T}) where {T<:MPolyDecRingElem}
    minimal_betti_table(A::MPolyQuoRing{T}) where {T<:MPolyDecRingElem}
    minimal_betti_table(I::MPolyIdeal{T}) where {T<:MPolyDecRingElem} 

Given a finitely presented graded module `M` over a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the Betti Table
of the minimal free resolution of `M`. Similarly for `A` and `I`.

# Examples
```julia
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [w^2-x*z, w*x-y*z, x^2-w*y, x*y-z^2, y^2-w*z]);

julia> A, _ = quo(R, I)
(Quotient of multivariate polynomial ring by ideal with 5 generators, Map from
R to A defined by a julia-function with inverse)

julia> minimal_betti_table(A)
       0  1  2  3
------------------
0    : 1  -  -  -
1    : -  5  5  -
2    : -  -  -  1
------------------
total: 1  5  5  1
```
"""
function minimal_betti_table(M::ModuleFP{T}) where {T<:MPolyDecRingElem}
 error("Not implemented for the given type")
end

function minimal_betti_table(M::SubquoModule{T}) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(M))
end

function minimal_betti_table(F::FreeMod{T}) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(M))
end

function minimal_betti_table(A::MPolyQuoRing{T}) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(A))
end

function minimal_betti_table(I::MPolyIdeal{T}) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(I))
end



@doc raw"""
    minimal_betti_table(F::FreeResolution{T}) where {T<:ModuleFP}

Given a graded free resolution `F` over a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the
Betti table of the minimal free resolution arising from `F`.

!!! note
    The algorithm proceeds without actually minimizing the resolution.

# Examples
```julia
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [w^2-x*z, w*x-y*z, x^2-w*y, x*y-z^2, y^2-w*z]);

julia> A, _ = quo(R, I)
(Quotient of multivariate polynomial ring by ideal with 5 generators, Map from
R to A defined by a julia-function with inverse)

julia> FA = free_resolution(A)
Free resolution of A
R^1 <---- R^5 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> betti_table(FA)
       0  1  2  3  
------------------
0    : 1  -  -  -  
1    : -  5  5  1  
2    : -  -  1  1  
------------------
total: 1  5  6  2  


julia> minimal_betti_table(FA)
       0  1  2  3  
------------------
0    : 1  -  -  -  
1    : -  5  5  -  
2    : -  -  -  1  
------------------
total: 1  5  5  1  
```
"""
function minimal_betti_table(res::FreeResolution{T}) where {T<:ModuleFP}
  @assert is_standard_graded(base_ring(res)) "resolution must be defined over a standard graded ring"
  @assert is_graded(res) "resolution must be graded"
  C = complex(res)
  @assert is_complete(res) "resolution must be complete"
  rng = range(C)
  # The following needs the resolution to be complete to be true
  res_length = first(rng)-1
  offsets = Dict{GrpAbFinGenElem, Int}()
  betti_hash_table = Dict{Tuple{Int, Any}, Int}()
  for i in 1:res_length+1
    phi = map(C, i)
    F = domain(phi)
    G = codomain(phi)
    dom_degs = unique!([degree(g) for g in gens(F)])
    cod_degs = unique!([degree(g) for g in gens(G)])
    for d in cod_degs
      d::GrpAbFinGenElem
      if d in dom_degs
        _, _, sub_mat = _constant_sub_matrix(phi, d)
        r = rank(sub_mat)
        c = ncols(sub_mat) - r - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
        offsets[d] = r
      else
        c = length(_indices_of_generators_of_degree(G, d)) - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
      end
    end
  end
  return BettiTable(betti_hash_table)
end

function hash_table(B::BettiTable) 
  return B.B
end

function generators_of_degree(
    C::FreeResolution{T},
    i::Int,
    d::GrpAbFinGenElem
  ) where {T<:ModuleFP}
  F = C[i]
  return [g for g in gens(F) if degree(g) == d]
end

function _indices_of_generators_of_degree(F::FreeMod{T}, d::GrpAbFinGenElem) where {T<:MPolyDecRingElem}
  result = Vector{Int}()
  for (i, g) in enumerate(gens(F))
    if degree(g) == d
      push!(result, i)
    end
  end
  return result
end

function _constant_sub_matrix(
    phi::FreeModuleHom{T, T},
    d::GrpAbFinGenElem
  ) where {RET<:MPolyDecRingElem{<:FieldElem}, T<:FreeMod{RET}}
  S = base_ring(domain(phi))::MPolyDecRing
  kk = coefficient_ring(S)::Field
  F = domain(phi)
  G = codomain(phi)
  ind_dom = _indices_of_generators_of_degree(F, d)
  ind_cod = _indices_of_generators_of_degree(G, d)
  m = length(ind_dom)
  n = length(ind_cod)
  result = zero(MatrixSpace(kk, m, n))
  for i in 1:m
    for j in 1:n
      c = phi(F[ind_dom[i]])[ind_cod[j]]
      result[i, j] = iszero(c) ? zero(kk) : first(coefficients(c))
    end
  end
  return ind_dom, ind_cod, result
end

# TODO: This will be provided soon from different sources.
function complex(F::FreeResolution) 
  return F.C
end

function base_ring(res::FreeResolution{T}) where {T<:ModuleFP}
  return base_ring(res[-1])
end
                                                                                                                                    
#############truncation#############

@doc raw"""
    truncate(M::ModuleFP, g::GrpAbFinGenElem, task::Symbol = :with_morphism)

Given a finitely presented graded module `M` over a $\mathbb Z$-graded multivariate 
polynomial ring with positive weights, return the truncation of `M` at degree `g`.

Put more precisely, return the truncation as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the inclusion map `N` $\to$ `M` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `M` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.

    truncate(M::ModuleFP, d::Int, task::Symbol = :with_morphism)

Given a module `M` as above, and given an integer `d`, convert `d` into an element `g`
of the grading group of `base_ring(I)` and proceed as above.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(R, 1)
Graded free module R^1([0]) of rank 1 over R

julia> V = [x*F[1]; y^4*F[1]; z^5*F[1]];

julia> M, _ = quo(F, V);

julia> M[1]
e[1]

julia> MT = truncate(M, 3);

julia> MT[1]
Graded subquotient of submodule of F generated by
1 -> z^3*e[1]
2 -> y*z^2*e[1]
3 -> y^2*z*e[1]
4 -> y^3*e[1]
5 -> x*z^2*e[1]
6 -> x*y*z*e[1]
7 -> x*y^2*e[1]
8 -> x^2*z*e[1]
9 -> x^2*y*e[1]
10 -> x^3*e[1]
by submodule of F generated by
1 -> x*e[1]
2 -> y^4*e[1]
3 -> z^5*e[1]
```
"""
function truncate(I::ModuleFP, g::GrpAbFinGenElem, task::Symbol = :with_morphism)
  return truncate(I, Int(g[1]), task)
end

function  truncate(I::ModuleFP, d::Int, task::Symbol = :with_morphism)
  @req I isa FreeMod || I isa SubquoModule "Not implemented for the given type"
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
  RES = elem_type(I)[]
  s = dmin
  B = monomial_basis(R, d-s)
  for i = 1:length(V)
    if degree(Int, V[i]) < d
      if degree(Int, V[i]) > s
        s = degree(Int, V[i])
        B = monomial_basis(R, d-s)
      end
      append!(RES, [x*V[i] for x in B])
    else
      push!(RES, V[i])
    end
  end
  return sub(I, RES, task) 
end


##################regularity#######################

@doc raw"""
    cm_regularity(M::ModuleFP)

Given a finitely presented graded module `M` over a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the
Castelnuovo-Mumford regularity of `M`.

    cm_regularity(I::MPolyIdeal)

Given a (homogeneous) ideal `I` in a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the
Castelnuovo-Mumford regularity of `I`.

# Examples
```julia
julia> R, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(R, 1);

julia> M, _ = quo(F, [x^2*F[1], y^2*F[1], z^2*F[1]])
(Graded subquotient of submodule of F generated by
1 -> e[1]
by submodule of F generated by
1 -> x^2*e[1]
2 -> y^2*e[1]
3 -> z^2*e[1], F -> M
e[1] -> e[1]
Homogeneous module homomorphism)

julia> cm_regularity(M)
3

julia> minimal_betti_table(M)
       0 1 2 3 
--------------
0    : 1 - - - 
1    : - 3 - - 
2    : - - 3 - 
3    : - - - 1 
--------------
total: 1 3 3 1 
```
```julia
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [-x*z+y^2, x*y-w*z, x^2-w*y]);

julia> cm_regularity(I)
2

julia> A, _ = quo(R, I);

julia> minimal_betti_table(A)
       0 1 2 
------------
0    : 1 - - 
1    : - 3 2 
------------
total: 1 3 2 
```
"""
function cm_regularity(M::ModuleFP)
 error("Not implemented for the given type")
end

function cm_regularity(M::FreeMod)
   R = base_ring(M)
   @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
   @req is_standard_graded(R) "The base ring is not standard ZZ-graded"
   return 0
end

function cm_regularity(M::SubquoModule)
   R = base_ring(M)
   @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
   @req is_standard_graded(R) "The base ring is not standard ZZ-graded"
   B = minimal_betti_table(M)
   S = as_dictionary(B)
   V = [x[2][1] - x[1] for x in keys(S)] 
  return maximum(V)
end

#####affine algebras as modules#####

@doc raw"""
    quotient_ring_as_module(A::MPolyQuoRing)

Return `A` considered as an object of type `SubquoModule`.

    quotient_ring_as_module(I::MPolyIdeal)

As above, where `A` is the quotient of `base_ring(I)` modulo `I`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x^2, y^3])
ideal(x^2, y^3)

julia> quotient_ring_as_module(I)
Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 2 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
```
```jldoctest
julia> S, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(S, [x^2, y^3])
ideal(x^2, y^3)

julia> quotient_ring_as_module(I)
Graded subquotient of submodule of S^1 generated by
1 -> e[1]
by submodule of S^1 generated by
1 -> x^2*e[1]
2 -> y^3*e[1]

```
"""
function quotient_ring_as_module(A::MPolyQuoRing)
  return quotient_ring_as_module(modulus(A))
end

function quotient_ring_as_module(I::MPolyIdeal)
  R = base_ring(I)
  F = is_graded(R) ? graded_free_module(R, 1) : free_module(R, 1)
  e1 = F[1]
  return quo(F, [x * e1 for x = gens(I)], :module)
end

#####ideals as modules#####

@doc raw"""
    ideal_as_module(I::MPolyIdeal)

Return `I` considered as an object of type `SubquoModule`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x^2, y^3])
ideal(x^2, y^3)

julia> ideal_as_module(I)
Submodule with 2 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
represented as subquotient with no relations.
```
```jldoctest
julia> S, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(S, [x^2, y^3])
ideal(x^2, y^3)

julia> ideal_as_module(I)
Graded submodule of S^1
1 -> x^2*e[1]
2 -> y^3*e[1]
represented as subquotient with no relations
```
"""
function ideal_as_module(I::MPolyIdeal)
  R = base_ring(I)
  F = is_graded(R) ? graded_free_module(R, 1) : free_module(R, 1)
  e1 = F[1]
  return sub(F, [x * e1 for x = gens(I)], :module)
end



##########################################################################
##### Twists
##########################################################################

@doc raw"""
    twist(M::ModuleFP{T}, g::GrpAbFinGenElem) where {T<:MPolyDecRingElem}

Return the twisted module `M(g)`.

    twist(M::ModuleFP{T}, W::Vector{<:IntegerUnion}) where {T<:MPolyDecRingElem}

Given a module `M` over a $\mathbb Z^m$-graded polynomial ring and a vector `W` of $m$ integers, 
convert `W` into an element `g` of the grading group of the ring and proceed as above.

    twist(M::ModuleFP{T}, d::IntegerUnion) where {T<:MPolyDecRingElem}

Given a module `M` over a $\mathbb Z$-graded polynomial ring and an integer `d`, 
convert `d` into an element `g` of the grading group of the ring and proceed as above.

# Examples
```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [zero(R)])
ideal(0)

julia> M = quotient_ring_as_module(I)
Graded submodule of R^1
1 -> e[1]
represented as subquotient with no relations

julia> degree(gen(M, 1))
[0]

julia> N = twist(M, 2)
Graded submodule of R^1
1 -> e[1]
represented as subquotient with no relations

julia> degree(gen(N, 1))
[-2]

```
"""
function twist(M::ModuleFP{T}, g::GrpAbFinGenElem) where {T<:MPolyDecRingElem}
 error("Not implemented for the given type")
end

function twist(M::SubquoModule{T}, g::GrpAbFinGenElem) where {T<:MPolyDecRingElem}
 R = base_ring(M)
 @req parent(g) == grading_group(R) "Group element not contained in grading group of base ring"
 F = ambient_free_module(M)
 FN = twist(F, g)
 GN = free_module(R, ngens(M))
 HN = free_module(R, length(relations(M)))
 a = hom(GN, F, ambient_representatives_generators(M))
 b = hom(HN, F, relations(M))
 A = matrix(a)
 B = matrix(b)
 N = subquotient(FN, A, B)
 return N
end

function twist(F::FreeMod{T}, g::GrpAbFinGenElem) where {T<:MPolyDecRingElem}
 R = base_ring(F)
 @req parent(g) == grading_group(R) "Group element not contained in grading group of base ring"
 W = [x-g for x in F.d]
 G = graded_free_module(R, rank(F))
 G.d = W
 return G
end

function twist(M::ModuleFP{T}, W::Vector{<:IntegerUnion}) where {T<:MPolyDecRingElem}
  R = base_ring(M)
  @assert is_zm_graded(R)
  return twist(M, grading_group(R)(W))
end

function twist(M::ModuleFP{T}, d::IntegerUnion) where {T<:MPolyDecRingElem}
  R = base_ring(M)
  @assert is_z_graded(R)
  return twist(M, grading_group(R)([d]))
end
