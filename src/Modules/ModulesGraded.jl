

###############################################################################
# Graded Modules constructors
###############################################################################

@doc raw"""
    graded_free_module(R::AdmissibleSparseFPModuleRing, p::Int, W::Vector{FinGenAbGroupElem}=[grading_group(R)[0] for i in 1:p], name::String="e")

Given a graded ring `R` with grading group `G`, say,
and given a vector `W` with `p` elements of `G`, create the free module $R^p$ 
equipped with its basis of standard unit vectors, and assign weights to these 
vectors according to the entries of `W`. Return the resulting graded free module.

    graded_free_module(R::AdmissibleSparseFPModuleRing, W::Vector{FinGenAbGroupElem}, name::String="e")

As above, with `p = length(W)`.

!!! note
    The function applies to graded multivariate polynomial rings and their quotients.

The string `name` specifies how the basis vectors are printed. 

# Examples
```jldoctest
julia> R, (x,y) = graded_polynomial_ring(QQ, [:x, :y])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> graded_free_module(R,3)
Graded free module R^3([0]) of rank 3 over R

julia> G = grading_group(R)
Z

julia> graded_free_module(R, [G[1], 2*G[1]])
Graded free module R^1([-1]) + R^1([-2]) of rank 2 over R
```
"""
function graded_free_module(R::AdmissibleSparseFPModuleRing, p::Int, W::Vector{FinGenAbGroupElem}=[grading_group(R)[0] for i in 1:p], name::String="e")
  @assert length(W) == p
  @assert is_graded(R)
  all(x -> parent(x) == grading_group(R), W) || error("entries of W must be elements of the grading group of the base ring")
  M = FreeMod(R, p, name)
  M.d = W
  return M
end

function graded_free_module(R::AdmissibleSparseFPModuleRing, p::Int, W::Vector{Any}, name::String="e")
  @assert length(W) == p
  @assert is_graded(R)
  p == 0 || error("W should be either an empty array or a Vector{FinGenAbGroupElem}")
  W = FinGenAbGroupElem[]
  return graded_free_module(R, p, W, name)
end

function graded_free_module(R::AdmissibleSparseFPModuleRing, W::Vector{FinGenAbGroupElem}, name::String="e")
  p = length(W)
  return graded_free_module(R, p, W, name)
end

function graded_free_module(R::AdmissibleSparseFPModuleRing, W::Vector{Any}, name::String="e")
  p = length(W)
  @assert is_graded(R)
  p == 0 || error("W should be either an empty array or a Vector{FinGenAbGroupElem}")
  W = FinGenAbGroupElem[]
  return graded_free_module(R, p, W, name)
end

@doc raw"""
    graded_free_module(R::AdmissibleSparseFPModuleRing, W::Vector{<:Vector{<:IntegerUnion}}, name::String="e")

Given a graded ring `R` with grading group $G = \mathbb Z^m$, 
and given a vector `W` of integer vectors of the same size `p`, say, create the free 
module $R^p$ equipped with its basis of standard unit vectors, and assign weights to these 
vectors according to the entries of `W`, converted to elements of `G`. Return the 
resulting graded free module.

    graded_free_module(R::AdmissibleSparseFPModuleRing, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, name::String="e")

As above, converting the columns of `W`.

    graded_free_module(R::AdmissibleSparseFPModuleRing, W::Vector{<:IntegerUnion}, name::String="e")

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
julia> R, (x,y) = graded_polynomial_ring(QQ, [:x, :y]);

julia> F = graded_free_module(R, [1, 2])
Graded free module R^1([-1]) + R^1([-2]) of rank 2 over R
```

```jldoctest
julia> S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1 0 1; 0 1 1]);

julia> FF = graded_free_module(S, [[1, 2], [-1, 3]])
Graded free module S^1([-1 -2]) + S^1([1 -3]) of rank 2 over S

julia> FFF = graded_free_module(S, [1 -1; 2 3])
Graded free module S^1([-1 -2]) + S^1([1 -3]) of rank 2 over S

julia> FF == FFF
true
```
"""
function graded_free_module(R::AdmissibleSparseFPModuleRing, W::Vector{<:Vector{<:IntegerUnion}}, name::String="e")
  @assert is_zm_graded(R)
  n = length(W[1])
  @assert all(x->length(x) == n, W)
  A = grading_group(R)
  d = [A(w) for w = W]
  return graded_free_module(R, length(W), d, name)
end

function graded_free_module(R::AdmissibleSparseFPModuleRing, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}, name::String="e")
  @assert is_zm_graded(R)
  A = grading_group(R)
  d = [A(W[:, i]) for i = 1:size(W, 2)]
  return graded_free_module(R, size(W, 2), d, name)
end

function graded_free_module(R::AdmissibleSparseFPModuleRing, W::Vector{<:IntegerUnion}, name::String="e")
  @assert is_graded(R)
  A = grading_group(R)
  d = [W[i] * A[1] for i in 1:length(W)]
  return graded_free_module(R, length(W), d, name)
end

@doc raw"""
    grade(F::FreeMod, W::Vector{FinGenAbGroupElem})

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
julia> R, x, y = polynomial_ring(QQ, :x => 1:2, :y => 1:3);

julia> G = abelian_group([0, 0])
Z^2

julia> g = gens(G)
2-element Vector{FinGenAbGroupElem}:
 [1, 0]
 [0, 1]

julia> W = [g[1], g[1], g[2], g[2], g[2]];

julia> S, _ = grade(R, W)
(Graded multivariate polynomial ring in 5 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x[1], x[2], y[1], y[2], y[3]])

julia> F = free_module(S, 3)
Free module of rank 3 over S

julia> FF = grade(F)
Graded free module S^3([0, 0]) of rank 3 over S

julia> F
Free module of rank 3 over S
```
"""
function grade(F::FreeMod, W::Vector{FinGenAbGroupElem})
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
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z],  [1 0 1; 0 1 1])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> F = free_module(R, 2)
Free module of rank 2 over R

julia> FF = grade(F,  [[1, 0], [0, 1]])
Graded free module R^1([-1 0]) + R^1([0 -1]) of rank 2 over R

julia> FFF = grade(F,  [1 0; 0 1])
Graded free module R^1([-1 0]) + R^1([0 -1]) of rank 2 over R
```

```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> S, _ = quo(R, [x*y])
(Quotient of multivariate polynomial ring by ideal (x*y), Map: R -> S)

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
julia> R, (x,y) = graded_polynomial_ring(QQ, [:x, :y]);

julia> F = graded_free_module(R, 3)
Graded free module R^3([0]) of rank 3 over R

julia> grading_group(F)
Z
```
"""
function grading_group(M::FreeMod)
  return grading_group(base_ring(M))
end

# Forgetful functor for gradings 
function forget_grading(F::FreeMod)
  @assert is_graded(F) "module must be graded"
  R = base_ring(F)
  result = FreeMod(R, ngens(F))
  phi = hom(F, result, gens(result); check=false)
  psi = hom(result, F, gens(F); check=false)
  set_attribute!(phi, :inverse=>psi)
  set_attribute!(psi, :inverse=>phi)
  return result, phi
end

function forget_grading(I::SubModuleOfFreeModule; 
    ambient_forgetful_map::FreeModuleHom=begin
      R = base_ring(I)
      F = ambient_free_module(I)
      _, iso_F = forget_grading(F)
      iso_F
    end
  )
  g = gens(I)
  gg = ambient_forgetful_map.(g)
  FF = codomain(ambient_forgetful_map)
  result = SubModuleOfFreeModule(FF, gg)
  return result
end

function forget_grading(M::SubquoModule;
    ambient_forgetful_map::FreeModuleHom=begin
      R = base_ring(M)
      F = ambient_free_module(M)
      _, iso_F = forget_grading(F)
      iso_F
    end
  )
  @assert is_graded(M) "module must be graded"
  FF = codomain(ambient_forgetful_map)
  if isdefined(M, :sub) && isdefined(M, :quo)
    new_sub = forget_grading(M.sub; ambient_forgetful_map)
    new_quo = forget_grading(M.quo; ambient_forgetful_map)
    result = SubquoModule(new_sub, new_quo)
    phi = hom(M, result, gens(result); check=false)
    psi = hom(result, M, gens(M); check=false)
    set_attribute!(phi, :inverse=>psi)
    set_attribute!(psi, :inverse=>phi)
    return result, phi
  elseif isdefined(M, :sub)
    new_sub = forget_grading(M.sub; ambient_forgetful_map)
    result = SubquoModule(new_sub)
    phi = hom(M, result, gens(result); check=false)
    psi = hom(result, M, gens(M); check=false)
    set_attribute!(phi, :inverse=>psi)
    set_attribute!(psi, :inverse=>phi)
    return result, phi
  elseif isdefined(M, :quo)
    new_quo = forget_grading(M.quo; ambient_forgetful_map)
    pre_result = SubquoModule(new_quo)
    result, _ = quo(FF, pre_result)
    phi = hom(M, result, gens(result); check=false)
    psi = hom(result, M, gens(M); check=false)
    set_attribute!(phi, :inverse=>psi)
    set_attribute!(psi, :inverse=>phi)
    return result, phi
  end
end


# Dangerous: Only for internal use with care!!!
@doc raw"""
    set_grading!(F::FreeMod, W::Vector{FinGenAbGroupElem})

    set_grading!(F::FreeMod, W::Vector{<:Vector{<:IntegerUnion}})

    set_grading!(F::FreeMod, W::Union{ZZMatrix, Matrix{<:IntegerUnion}}) 

    set_grading!(F::FreeMod, W::Vector{<:IntegerUnion})

Assign weights to the generators of `F` according to the entries of `W`.

See the `grade` and `graded_free_module` functions.
```
"""
function set_grading!(M::FreeMod, W::Vector{FinGenAbGroupElem})
  @assert length(W) == ngens(M)
  @assert is_graded(base_ring(M))
  R = base_ring(M)
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

function degrees(M::FreeMod; check::Bool=true)
  @assert is_graded(M)
  return M.d::Vector{FinGenAbGroupElem}
end

@doc raw"""
    degrees_of_generators(F::FreeMod)

Return the degrees of the generators of `F`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(R, 2)
Graded free module R^2([0]) of rank 2 over R

julia> degrees_of_generators(F)
2-element Vector{FinGenAbGroupElem}:
 [0]
 [0]
```
"""
function degrees_of_generators(F::FreeMod; check::Bool=true)
  return degrees(F; check)
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
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z], [1, 2, 3]);

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
  !isnothing(el.d) && return true
  !is_graded(parent(el)) && error("the parent module is not graded")
  iszero(el) && return true
  el.d = isa(el.d, FinGenAbGroupElem) ? el.d : determine_degree_from_SR(coordinates(el), degrees(parent(el)))
  return isa(el.d, FinGenAbGroupElem)
end

const AnyGradedRingElem = Union{<:MPolyDecRingElem, <:MPolyQuoRingElem{<:MPolyDecRingElem},
                                <:MPolyLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing},
                                <:MPolyQuoLocRingElem{<:Ring, <:RingElem, <:MPolyDecRing}
                               }

@doc raw"""
    degree(f::FreeModElem{T}; check::Bool=true) where {T<:AnyGradedRingElem}

Given a homogeneous element `f` of a graded free module, return the degree of `f`.

    degree(::Type{Vector{Int}}, f::FreeModElem)

Given a homogeneous element `f` of a $\mathbb Z^m$-graded free module, return the degree of `f`, converted to a vector of integer numbers.

    degree(::Type{Int}, f::FreeModElem)

Given a homogeneous element `f` of a $\mathbb Z$-graded free module, return the degree of `f`, converted to an integer number.

If `check` is set to `false`, then there is no check for homegeneity. This should be called 
internally on provably sane input, as it speeds up computation significantly. 
# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> f = y^2*z − x^2*w
-w*x^2 + y^2*z

julia> degree(f)
[3]

julia> typeof(degree(f))
FinGenAbGroupElem

julia> degree(Int, f)
3

julia> typeof(degree(Int, f))
Int64
```
"""
function degree(f::FreeModElem{T}; check::Bool=true) where {T<:AnyGradedRingElem}
  !isnothing(f.d) && return f.d::FinGenAbGroupElem
  @check is_graded(parent(f)) "the parent module is not graded"
  @check is_homogeneous(f) "the element is not homogeneous"
  f.d = _degree_fast(f)
  return f.d::FinGenAbGroupElem
end

function _degree_of_parent_generator(f::FreeModElem, i::Int)
  return f.parent.d[i]::FinGenAbGroupElem
end

# TODO: This has the potential to be a "hot" function. 
# Should we store the information in the parent of `f` directly?
# Or is it enough that things are cached in the generators 
# of the `sub`?
function _degree_of_parent_generator(f::SubquoModuleElem, i::Int)
  return _degree_fast(gen(parent(f), i))::FinGenAbGroupElem
end

# Fast method only to be used on sane input; returns a `GrbAbFinGenElem`.
# This is exposed as an extra internal function so that `check=false` can be avoided. 
function _degree_fast(f::FreeModElem)
  iszero(f) && return zero(grading_group(base_ring(f)))
  for (i, c) in coordinates(f)
    !iszero(c) && return (_degree_fast(c) + _degree_of_parent_generator(f, i))::FinGenAbGroupElem
  end
  error("this line should never be reached")
end


function degree(::Type{Vector{Int}}, f::FreeModElem; check::Bool=true)
  @assert is_zm_graded(parent(f))
  d = degree(f; check)
  return Int[d[i] for i=1:ngens(parent(d))]
end

function degree(::Type{Int}, f::FreeModElem; check::Bool=true)
  @assert is_z_graded(parent(f))
  return Int(degree(f; check)[1])
end

# Checks for homogeneity and computes the degree. 
# If the input is not homogeneous, this returns nothing.
function determine_degree_from_SR(coords::SRow, unit_vector_degrees::Vector{FinGenAbGroupElem})
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

function graded_map(F::FreeMod{T}, A::MatrixElem{T}; check::Bool=true) where {T <: AdmissibleSparseFPModuleRingElem}
  R = base_ring(F)
  G = grading_group(R)
  source_degrees = Vector{eltype(G)}()
  for i in 1:nrows(A)
      for j in 1:ncols(A)
          if !is_zero(A[i, j])
              push!(source_degrees, degree(A[i, j]; check) + degree(F[j]; check))
              break
          end
      end
  end
  Fcdm = graded_free_module(R, source_degrees)
  phi = hom(Fcdm, F, A; check)
  return phi
end

function graded_map(F::FreeMod{T}, V::Vector{<:AbstractFreeModElem{T}}; check::Bool=true) where {T <: AdmissibleSparseFPModuleRingElem}
  R = base_ring(F)
  G = grading_group(R)
  nrows = length(V)
  ncols = rank(F)
  @check true # Trigger an error if checks are supposed to be disabled.
  
  source_degrees = Vector{eltype(G)}()
  for (i, v) in enumerate(V)
    if is_zero(v)
      push!(source_degrees, zero(G))
      continue
    end
    for (j, c) in coordinates(v)
      if !iszero(c)
        push!(source_degrees, degree(coordinates(V[i])[j]; check) + degree(F[j]; check))
        break
      end
    end
  end
  @assert length(source_degrees) == nrows
  Fcdm = graded_free_module(R, source_degrees)
  @assert ngens(Fcdm) == length(V) "number of generators must be equal to the number of images"
  phi = hom(Fcdm, F, V; check)
  return phi
end


function graded_map(F::SubquoModule{T}, V::Vector{<:SparseFPModuleElem{T}}; check::Bool=true) where {T <: AdmissibleSparseFPModuleRingElem}
  R = base_ring(F)
  G = grading_group(R)
  nrows = length(V)
  source_degrees = Vector{eltype(G)}()
  for (i, v) in enumerate(V)
    if is_zero(v)
      push!(source_degrees, zero(G))
      continue
    end
    for (j, c) in coordinates(v)
      if !iszero(c)
        push!(source_degrees, degree(coordinates(V[i])[j]; check) + degree(F[j]; check))
        break
      end
    end
  end
  Fcdm = graded_free_module(R, source_degrees)
  phi = hom(Fcdm, F, V; check)
  return phi
end

###############################################################################
# Graded Free Module homomorphisms functions
###############################################################################

function set_grading(f::FreeModuleHom{T1, T2}; check::Bool=true) where {T1 <: FreeMod, T2 <: Union{FreeMod, SubquoModule, Oscar.SubModuleOfFreeModule}}
  if !is_graded(domain(f)) || !is_graded(codomain(f))
      return f
  end
  f.d = degree(f; check)
  return f
end

function set_grading(f::FreeModuleHom{T1, T2}; check::Bool=true) where {T1 <: FreeMod_dec, T2 <: FreeMod_dec}
  return f
end
# for decorations: add SubquoModule_dec for codomain once it exists

@doc raw"""
    degree(a::FreeModuleHom; check::Bool=true)

If `a` is graded, return the degree of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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
Graded module homomorphism of degree [1]
  from F
  to G
defined by
  e[1] -> y*e[1]
  e[2] -> x*e[1] + y*e[2]
  e[3] -> z*e[2]

julia> degree(a)
[1]
```
"""
function degree(f::FreeModuleHom; check::Bool=true)
  # TODO: isdefined should not be necessary here. Can it be kicked?
  isdefined(f, :d) && isnothing(f.d) && return nothing # This stands for the map being not homogeneous
  isdefined(f, :d) && return f.d::FinGenAbGroupElem

  @check (is_graded(domain(f)) && is_graded(codomain(f))) "both domain and codomain must be graded"
  @check is_graded(f) "map is not graded"
  for i in 1:ngens(domain(f))
    if iszero(domain(f)[i]) || iszero(image_of_generator(f, i))
      continue
    end
    f.d = degree(image_of_generator(f, i); check) - degree(domain(f)[i]; check)
    return f.d::FinGenAbGroupElem
  end

  # If we got here, the map is the zero map. Return degree zero in this case
  return zero(grading_group(domain(f)))::FinGenAbGroupElem

  # Old code left for debugging
  return degree(image_of_generator(f, 1))
  domain_degrees = degrees(T1)
  df = nothing
  for i in 1:length(domain_degrees)
    image_vector = f(T1[i])
    if isempty(coordinates(image_vector)) || is_zero(image_vector)
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
    is_graded(a::SparseFPModuleHom)

Return `true` if `a` is graded, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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
Graded module homomorphism of degree [1]
  from F
  to G
defined by
  e[1] -> y*e[1]
  e[2] -> x*e[1] + y*e[2]
  e[3] -> z*e[2]

julia> is_graded(a)
true
```
"""
function is_graded(f::SparseFPModuleHom)
  isdefined(f, :d) && return true
  (is_graded(domain(f)) && is_graded(codomain(f))) || return false
  T1 = domain(f)
  T2 = codomain(f)
  domain_degrees = degrees_of_generators(T1)
  df = nothing
  for i in 1:length(domain_degrees)
    image_vector = f(T1[i])
    if isempty(coordinates(image_vector)) || is_zero(image_vector)
      continue
    end
    current_df = degree(image_vector) - domain_degrees[i]
    if df === nothing
      df = current_df
    elseif df != current_df
      return false
    end
  end
  if df === nothing
    R = base_ring(T1)
    G = grading_group(R)
    f.d = zero(G)
    return true
  end
  f.d = df
  return true
end

@doc raw"""
    grading_group(a::FreeModuleHom)

If `a` is graded, return the grading group of `a`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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
Graded module homomorphism of degree [1]
  from F
  to G
defined by
  e[1] -> y*e[1]
  e[2] -> x*e[1] + y*e[2]
  e[3] -> z*e[2]

julia> is_graded(a)
true

julia> grading_group(a)
Z
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
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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
Graded module homomorphism of degree [1]
  from F
  to G
defined by
  e[1] -> y*e[1]
  e[2] -> x*e[1] + y*e[2]
  e[3] -> z*e[2]

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


function degrees_of_generators(M::SubModuleOfFreeModule{T}; check::Bool=true) where T
  return map(gen -> degree(gen; check), gens(M))
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
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F1 = graded_free_module(Rg, [2,2,2]);

julia> F2 = graded_free_module(Rg, [2]);

julia> G = graded_free_module(Rg, [1,1]);

julia> V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]];

julia> V2 = [z*G[2]+y*G[1]];

julia> a1 = hom(F1, G, V1);

julia> a2 = hom(F2, G, V2);

julia> M = subquotient(a1,a2);

julia> grading_group(M)
Z
```
"""
function grading_group(M::SubquoModule)
  return grading_group(base_ring(M))
end

@doc raw"""
    degrees_of_generators(M::SubquoModule; check::Bool=true)

Return the degrees of the generators of `M`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F1 = graded_free_module(Rg, [2,2,2]);

julia> F2 = graded_free_module(Rg, [2]);

julia> G = graded_free_module(Rg, [1,1]);

julia> V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]];

julia> V2 = [z*G[2]+y*G[1]];

julia> a1 = hom(F1, G, V1);

julia> a2 = hom(F2, G, V2);

julia> M = subquotient(a1,a2);

julia> degrees_of_generators(M)
3-element Vector{FinGenAbGroupElem}:
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
function degrees_of_generators(M::SubquoModule{T}; check::Bool=true) where T
  isempty(gens(M)) ? FinGenAbGroupElem[] : map(gen -> degree(repres(gen); check), gens(M))
end

###############################################################################
# Graded subquotient elements
###############################################################################

 @doc raw"""
    is_homogeneous(m::SubquoModuleElem)

Return  `true` if `m` is homogeneous, `false` otherwise.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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
  el.is_reduced && return is_homogeneous(repres(el))
  return is_homogeneous(repres(simplify(el)))

  # The following call checks for homogeneity on the way and stores the degree thus determined.
  degree = determine_degree_from_SR(coordinates(el), degrees_of_generators(parent(el)))
  if degree === nothing
    reduced_el = simplify(el) # TODO: What do we expect `simplify` to do here generically???
    degree_reduced = determine_degree_from_SR(coordinates(reduced_el), degrees_of_generators(parent(reduced_el)))
    if degree_reduced === nothing
      el.d = nothing
      return false
    else
      el.d = degree_reduced
      return true
    end
  else
    el.d = degree
    return true
  end
end

@doc raw"""
    degree(m::SubquoModuleElem; check::Bool=true)

Given a homogeneous element `m` of a graded subquotient, return the degree of `m`.

    degree(::Type{Vector{Int}}, m::SubquoModuleElem)

Given a homogeneous element `m` of a $\mathbb Z^m$-graded subquotient, return the degree of `m`, converted to a vector of integer numbers.

    degree(::Type{Int}, m::SubquoModuleElem)

Given a homogeneous element `m` of a $\mathbb Z$-graded subquotient, return the degree of `m`, converted to an integer number.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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
function degree(el::SubquoModuleElem; check::Bool=true)
  # In general we can not assume that we have a groebner basis reduction available 
  # as a backend to bring the element to normal form. 
  # In particular, this may entail that `coordinates` produces non-homogeneous 
  # vectors via differently implemented liftings. 
  # Thus, the only thing we can do is to assume that the representative is 
  # homogeneous.
  return degree(repres(simplify!(el)); check)
end

# When there is a Groebner basis backend, we can reduce to normal form.
function degree(
    el::SubquoModuleElem{T};
    check::Bool=true
  ) where {T <:Union{<:MPolyRingElem{<:FieldElem}}}
  !el.is_reduced && return degree(simplify!(el); check)
  return degree(repres(el); check)
end

_degree_fast(el::SubquoModuleElem) = degree(el, check=false)

function degree(::Type{Vector{Int}}, el::SubquoModuleElem; check::Bool=true)
  @assert is_zm_graded(parent(el))
  d = degree(el; check)
  return Int[d[i] for i=1:ngens(parent(d))]
end

function degree(::Type{Int}, el::SubquoModuleElem; check::Bool=true)
  @assert is_z_graded(parent(el))
  return Int(degree(el; check)[1])
end


###############################################################################
# Graded subquotient homomorphisms functions
###############################################################################

function set_grading(f::SubQuoHom; check::Bool=true)
  if !is_graded(domain(f)) || !is_graded(codomain(f))
    return(f)
  end
  f.d = degree(f; check)
  return f
end

@doc raw"""
    degree(a::SubQuoHom; check::Bool=true)

If `a` is graded, return the degree of `a`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

julia> degree(a)
[2]
```
"""
function degree(f::SubQuoHom; check::Bool=true)
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
  @check is_graded(f) "homomorphism is not graded"
  for (i, v) in enumerate(gens(T1))
    (is_zero(v) || is_zero(image_of_generator(f, i))) && continue
    f.d = degree(image_of_generator(f, i); check) - degree(v; check)
    return f.d::FinGenAbGroupElem
  end

  # If we get here, we have the zero map
  f.d = zero(grading_group(T1))
  return f.d::FinGenAbGroupElem

  # Old code left here for debugging
  domain_degrees = degrees_of_generators(T1)
  df = nothing
  for i in 1:length(domain_degrees)
    image_vector = f(T1[i])
    if isempty(coordinates(image_vector)) || is_zero(image_vector)
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
    grading_group(a::SubQuoHom)

If `a` is graded, return the grading group of `a`.

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

julia> grading_group(a)
Z
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
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

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
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> I = ideal(R, [x*z, y*z, x*w^2, y*w^2])
Ideal generated by
  x*z
  y*z
  w^2*x
  w^2*y

julia> A, _= quo(R, I)
(Quotient of multivariate polynomial ring by ideal (x*z, y*z, w^2*x, w^2*y), Map: R -> A)

julia> FA  = free_resolution(A)
Free resolution of A
R^1 <---- R^4 <---- R^4 <---- R^1 <---- 0
0         1         2         3         4

julia> betti_table(FA)
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  2  1  -
     2: -  2  3  1
------------------
 total: 1  4  4  1

julia> R, (x, y) = graded_polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x, y, x+y]);

julia> M = quotient_ring_as_module(I);

julia> FM = free_resolution(M, algorithm = :nres);

julia> betti_table(FM)
degree: 0  1  2
---------------
     0: 1  2  1
---------------
 total: 1  2  1
```
"""
function betti_table(F::FreeResolution; project::Union{FinGenAbGroupElem, Nothing} = nothing, reverse_direction::Bool=false)
  @assert is_graded(F) "resolution must be graded"
  generator_count = Dict{Tuple{Int, Any}, Int}()
  C = F.C
  rng = Hecke.map_range(C)
  n = first(rng)
  for i in 0:n
    module_degrees = F[i].d
    for degree in module_degrees
      idx = (i, degree)
      generator_count[idx] = get(generator_count, idx, 0) + 1
    end
  end
  return BettiTable(generator_count, project = project, reverse_direction = reverse_direction)
end

function betti(b::FreeResolution; project::Union{FinGenAbGroupElem, Nothing} = nothing, reverse_direction::Bool = false)
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

#numbers are right aligned in the column
function Base.show(io::IO, b::BettiTable)
  if isempty(b.B)
    print(io, "Empty table")
    return
  end

  T = induce_shift(b.B)
  x = collect(keys(T))
  if isempty(x)
    print(io, "Empty table")
    return
  end
  step, min, maxv = b.reverse_direction ? (-1, maximum(first, x), minimum(first, x)) : (1, minimum(first, x), maximum(first, x))
  column_widths = Dict()
  for j in min:step:maxv
    sum_col = sum(getindex(T, x[m]) for m in 1:length(x) if x[m][1] == j)
    col_width_from_sum = ndigits(abs(sum_col))
    col_width_from_header = ndigits(abs(j))# + (j < 0 ? 1 : 0)
    column_widths[j] = max(col_width_from_sum, col_width_from_header) + 2
  end

  if b.project === nothing
    for i in 1:ngens(parent(x[1][2]))
      ngens(parent(x[1][2])) > 1 && println(io, "Betti Table for component ", i)

      # figure out width of first column
      L = sort(unique!(collect(x[k][2][i] for k in 1:length(x))))
      mi = minimum(L)
      mx = maximum(L)
      # 6 = length(degree); we take length of mi into account in case it is negative
      first_column_width = max(6, ndigits(mi), ndigits(mx))

      # header row
      print(io, lpad("degree", first_column_width), ":")
      total_space_count = first_column_width
      for j in min:step:maxv
        adjustment = j < 0 ? 1 : 0
        if j == min
          adjustment += 1
        end
        space_count = max(0, column_widths[j] - ndigits(j) - adjustment)
        print(io, " "^(space_count))
        print(io, j)
        total_space_count += space_count + ndigits(j) + adjustment
      end
      print(io, "\n")
      # separator row
      println(io, "-"^total_space_count)
      # print bulk of table
      for j in mi:mx
        print(io, lpad(j, first_column_width), ":")
        for h in min:step:maxv
          sum_current = sum([getindex(T, x[k]) for k in 1:length(x) if x[k][1] == h && x[k][2][i] == j])
          @assert column_widths[h] - ndigits(sum_current) >= 2
          print(io, " "^(column_widths[h] - ndigits(sum_current) - 2))
          print(io, " ", sum_current == 0 ? "-" : sum_current)
          if h!= maxv
            print(io, " ")
          end
        end
        print(io,"\n")
      end
      # separator row
      println(io, "-"^total_space_count)
      # footer row
      print(io, lpad("total", first_column_width), ":")
      for i_total in min:step:maxv
        sum_row = sum(getindex(T, x[j]) for j in 1:length(x) if x[j][1] == i_total)
        print(io, " ", sum_row)
        if i_total != maxv
          print(io, " "^(column_widths[i_total] - ndigits(sum_row)-1))
        end
      end
    end
  else
    parent(b.project) == parent(x[1][2]) || error("projection vector has wrong type")
    print(io, "Betti Table for scalar product of grading with ", coordinates(b.project), "\n")
    print(io, "  ")
    L = Vector{ZZRingElem}(undef,0)
    for i in 1:length(x)
      temp_sum = (coordinates(b.project) * transpose(coordinates(x[i][2])))[1]
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
              current_sum = (coordinates(b.project) * transpose(coordinates(x[i][2])))[1]
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
  row_ind = size(st.values, 1) - ind[1]
  col_ind = ind[2] - first(st.twist_range) + 1
  return st.values[row_ind, col_ind]
end

function Base.show(io::IO, table::sheafCohTable)
  chi = [any(==(-1), col) ? -1 : sum(col) for col in eachcol(table.values)]
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
    pushfirst!(print_rows[i], rpad("$(nrows-i): ", row_label_length, " "))
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
    sheaf_cohomology(M::SparseFPModule{T}, l::Int, h::Int; algorithm::Symbol = :bgg) where {T <: MPolyDecRingElem}

If `M` is a graded module over a standard graded multivariate polynomial ring with coefficients in a field `K`, 
say, and $\mathcal F = \widetilde{M}$ is the coherent sheaf associated to `M` on the corresponding projective 
space $\mathbb P^n(K)$, consider the cohomology groups $H^i(\mathbb P^n(K), \mathcal F(d))$ as vector spaces 
over $K$, and return their dimensions $h^i(\mathbb P^n(K), \mathcal F(d))$ in the range of twists $d$ 
indicated by `l` and `h`. The result is presented as a table, where '-' indicates that
$h^i(\mathbb P^n(K), \mathcal F(d)) = 0$. The line starting  with `chi` lists the Euler characteristic 
of each twist under consideration. The values in the table can be accessed as shown in the 
first example below. Note that this example addresses the cotangent bundle on projective 3-space, while the 
second example is concerned with the structure sheaf of projective 4-space.

The keyword `algorithm` can be set to
- `:bgg` (use the Tate resolution via the Bernstein-Gelfand-Gelfand correspondence),
- `:loccoh` (use local cohomology).

!!! note 
    Due to the shape of the Tate resolution, the algorithm addressed by `bgg` does not compute all values in the given range `l` $<$ `h`. The missing values are indicated by a `*`. To determine all values in the range `l` $<$ `h`, enter `sheaf_cohomology(M, l-ngens(base_ring(M)), h+ngens(base_ring(M)))`.

```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:4);

julia> S, _= grade(R);

julia> I = ideal(S, gens(S))
Ideal generated by
  x[1]
  x[2]
  x[3]
  x[4]

julia> FI = free_resolution(I)
Free resolution of I
S^4 <---- S^6 <---- S^4 <---- S^1 <---- 0
0         1         2         3         4

julia> M = cokernel(map(FI, 2));

julia> tbl = sheaf_cohomology(M, -6, 2)
twist:  -6  -5  -4  -3  -2  -1   0   1   2
------------------------------------------
3:      70  36  15   4   -   -   -   -   *
2:       *   -   -   -   -   -   -   -   -
1:       *   *   -   -   -   -   1   -   -
0:       *   *   *   -   -   -   -   -   6
------------------------------------------
chi:     *   *   *   4   -   -   1   -   *

julia> tbl[0, 2]
6

julia> tbl[1, 0]
1

julia> sheaf_cohomology(M, -9, 5)
twist:   -9   -8   -7   -6   -5   -4   -3   -2   -1    0    1    2    3    4    5
---------------------------------------------------------------------------------
3:      280  189  120   70   36   15    4    -    -    -    -    -    *    *    *
2:        *    -    -    -    -    -    -    -    -    -    -    -    -    *    *
1:        *    *    -    -    -    -    -    -    -    1    -    -    -    -    *
0:        *    *    *    -    -    -    -    -    -    -    -    6   20   45   84
---------------------------------------------------------------------------------
chi:      *    *    *   70   36   15    4    -    -    1    -    6    *    *    *
```

```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:5);

julia> S, _  = grade(R);

julia> F = graded_free_module(S, 1);

julia> sheaf_cohomology(F, -8, 3, algorithm = :loccoh)
twist:  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3
------------------------------------------------------
4:      35  15   5   1   -   -   -   -   -   -   -   -
3:       -   -   -   -   -   -   -   -   -   -   -   -
2:       -   -   -   -   -   -   -   -   -   -   -   -
1:       -   -   -   -   -   -   -   -   -   -   -   -
0:       -   -   -   -   -   -   -   -   1   5  15  35
------------------------------------------------------
chi:    35  15   5   1   -   -   -   -   1   5  15  35
```
"""
function sheaf_cohomology(M::SparseFPModule{T},
                          l::Int,
                          h::Int;
                          algorithm::Symbol = :bgg) where {T <: MPolyDecRingElem}
  if algorithm == :bgg
    return _sheaf_cohomology_bgg(M, l, h)
  elseif algorithm == :loccoh
    return _sheaf_cohomology_loccoh(M, l, h)
  else
    error("Algorithm not supported.")
  end
end

@doc raw"""
    _sheaf_cohomology_bgg(M::SparseFPModule{T}, l::Int, h::Int) where {T <: MPolyDecRingElem}

Compute the cohomology of twists of of the coherent sheaf on projective
space associated to `M`. The method used is based on the Bernstein-Gelfand-Gelfand correspondence. The range of twists is between `l` and `h`.
In the displayed result, '-' refers to a zero entry and '*' refers to a
negative entry (= dimension not yet determined). To determine all values
in the desired range between `l` and `h` use `_sheaf_cohomology_bgg(M, l-ngens(base_ring(M)), h+ngens(base_ring(M)))`.
The values of the returned table can be accessed by indexing it
with a cohomological index and a value between `l` and `h` as shown
in the example below.

```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:5);

julia> R, x = grade(R);

julia> F = graded_free_module(R, 1);

julia> Oscar._sheaf_cohomology_bgg(F, -7, 2)
twist:  -7  -6  -5  -4  -3  -2  -1   0   1   2
----------------------------------------------
4:      15   5   1   -   -   -   *   *   *   *
3:       *   -   -   -   -   -   -   *   *   *
2:       *   *   -   -   -   -   -   -   *   *
1:       *   *   *   -   -   -   -   -   -   *
0:       *   *   *   *   -   -   -   1   5  15
----------------------------------------------
chi:     *   *   *   *   -   -   *   *   *   *

julia> sheaf_cohomology(F, -11, 6)
twist:  -11  -10   -9   -8   -7   -6   -5   -4   -3   -2   -1    0    1    2    3    4    5    6
------------------------------------------------------------------------------------------------
4:      210  126   70   35   15    5    1    -    -    -    -    -    -    -    *    *    *    *
3:        *    -    -    -    -    -    -    -    -    -    -    -    -    -    -    *    *    *
2:        *    *    -    -    -    -    -    -    -    -    -    -    -    -    -    -    *    *
1:        *    *    *    -    -    -    -    -    -    -    -    -    -    -    -    -    -    *
0:        *    *    *    *    -    -    -    -    -    -    -    1    5   15   35   70  126  210
------------------------------------------------------------------------------------------------
chi:      *    *    *    *   15    5    1    -    -    -    -    1    5   15    *    *    *    *
```

```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:4);

julia> S, _= grade(R);

julia> I = ideal(S, gens(S))
Ideal generated by
  x[1]
  x[2]
  x[3]
  x[4]

julia> FI = free_resolution(I)
Free resolution of I
S^4 <---- S^6 <---- S^4 <---- S^1 <---- 0
0         1         2         3         4

julia> M = cokernel(map(FI, 2));

julia> tbl = sheaf_cohomology(M, -6, 2, algorithm = :loccoh)
twist:  -6  -5  -4  -3  -2  -1   0   1   2
------------------------------------------
3:      70  36  15   4   -   -   -   -   -
2:       -   -   -   -   -   -   -   -   -
1:       -   -   -   -   -   -   1   -   -
0:       -   -   -   -   -   -   -   -   6
------------------------------------------
chi:    70  36  15   4   -   -   1   -   6

julia> tbl[3, -6]
70

julia> tbl[1, 0]
1
```
"""
function _sheaf_cohomology_bgg(M::SparseFPModule{T},
                               l::Int,
                               h::Int) where {T <: MPolyDecRingElem}

  sing_mod, weights = _weights_and_sing_mod(M)
  reg = Int(cm_regularity(M; check=false))

  values = Singular.LibSheafcoh.sheafCohBGGregul_w(sing_mod,
                                                   l, h, reg,
                                                   weights)
  return sheafCohTable(l:h, values)
end

@doc raw"""
    _sheaf_cohomology_loccoh(M::SparseFPModule{T}, l::Int, h::Int) where {T <: MPolyDecRingElem}

Compute the cohomology of twists of of the coherent sheaf on projective
space associated to `M` The method used is based on local duality. The range of twists is between `l` and `h`.
In the displayed result, '-' refers to a zero entry and '*' refers to a
negative entry (= dimension not yet determined). To determine all values
in the desired range between `l` and `h` use `_sheaf_cohomology_loccoh(M, l-ngens(base_ring(M)), h+ngens(base_ring(M)))`.
The values of the returned table can be accessed by indexing it
with a cohomological index and a value between `l` and `h` as shown
in the example below.

```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:4);

julia> S, _= grade(R);

julia> I = ideal(S, gens(S))
Ideal generated by
  x[1]
  x[2]
  x[3]
  x[4]

julia> FI = free_resolution(I)
Free resolution of I
S^4 <---- S^6 <---- S^4 <---- S^1 <---- 0
0         1         2         3         4

julia> M = cokernel(map(FI, 2));

julia> tbl = Oscar._sheaf_cohomology_loccoh(M, -6, 2)
twist:  -6  -5  -4  -3  -2  -1   0   1   2
------------------------------------------
3:      70  36  15   4   -   -   -   -   -
2:       -   -   -   -   -   -   -   -   -
1:       -   -   -   -   -   -   1   -   -
0:       -   -   -   -   -   -   -   -   6
------------------------------------------
chi:    70  36  15   4   -   -   1   -   6

julia> tbl[3, -6]
70

julia> tbl[1, 0]
1

julia> R, x = polynomial_ring(QQ, :x => 1:5);

julia> R, x = grade(R);

julia> F = graded_free_module(R, 1);

julia> Oscar._sheaf_cohomology_loccoh(F, -7, 2)
twist:  -7  -6  -5  -4  -3  -2  -1   0   1   2
----------------------------------------------
4:      15   5   1   -   -   -   -   -   -   -
3:       -   -   -   -   -   -   -   -   -   -
2:       -   -   -   -   -   -   -   -   -   -
1:       -   -   -   -   -   -   -   -   -   -
0:       -   -   -   -   -   -   -   1   5  15
----------------------------------------------
chi:    15   5   1   -   -   -   -   1   5  15
```
"""
function _sheaf_cohomology_loccoh(M::SparseFPModule{T},
                                  l::Int,
                                  h::Int) where {T <: MPolyDecRingElem}

  sing_mod, weights = _weights_and_sing_mod(M)

  values = Singular.LibSheafcoh.sheafCoh_w(sing_mod,
                                           l, h,
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

function _weights_and_sing_mod(M::SparseFPModule{T}) where {T <: MPolyDecRingElem}

  CR = base_ring(base_ring(M))
  @assert isa(CR, AbstractAlgebra.Field) "Base ring of input module must be defined over a field."
  free_mod = ambient_free_module(M)
  @assert is_standard_graded(free_mod) "Only supported for the standard grading of polynomials."
  @assert is_graded(M) "Module must be graded."

  # get a cokernel presentation of M
  p = presentation(M)
  #cokern_repr = image(map(p, 1))[1] # Creating the inclusion map takes too long, see Issue #2999
  cokern_repr = SubquoModule(p[0], elem_type(p[0])[map(p, 1)(v) for v in gens(p[1])])
  cokern_gens = ambient_representatives_generators(cokern_repr)
  if isempty(cokern_gens)
    cokern_gens = [zero(ambient_free_module(cokern_repr))]
  end
  sing_mod = singular_generators(ModuleGens(cokern_gens))
  weights = [Int(d[1]) for d in degrees_of_generators(p[0])]
  return sing_mod, weights
end
  
##################################
### Tests on graded modules
##################################

function is_graded(M::FreeMod)
  return isa(M.d, Vector{FinGenAbGroupElem})
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
julia> R, (x,y) = graded_polynomial_ring(QQ, [:x, :y])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> free_module_dec(R,3)
Decorated free module of rank 3 over R

```
"""
free_module_dec(R::CRing_dec, n::Int, name::VarName = :e; cached::Bool = false) = FreeMod_dec(R, n, name, cached = cached)


@doc raw"""
    FreeMod_dec(R::CRing_dec, d::Vector{FinGenAbGroupElem}, name::VarName = :e; cached::Bool = false) 

Construct a decorated (graded or filtered) free module over the ring `R` 
with rank `n` where `n` is the length of `d`. `d` is the vector of degrees for the 
components, i.e. `d[i]` is the degree of `e[i]` where `e[i]` is the `i`th standard unit
vector of the free module.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for 
the standard unit vectors.
"""
function FreeMod_dec(R::CRing_dec, d::Vector{FinGenAbGroupElem}, name::VarName = :e; cached::Bool = false) 
  return FreeMod_dec{elem_type(R)}(R, [Symbol("$name[$i]") for i=1:length(d)],d)
end

@doc raw"""
    free_module_dec(R::CRing_dec, d::Vector{FinGenAbGroupElem}, name::VarName = :e; cached::Bool = false)

Create the decorated free module $R^n$ (`n` is the length of `d`)
equipped with its basis of standard unit vectors where the 
i-th standard unit vector has degree `d[i]`.

The string `name` specifies how the basis vectors are printed. 
"""
free_module_dec(R::CRing_dec, d::Vector{FinGenAbGroupElem}, name::VarName = :e; cached::Bool = false) = FreeMod_dec(R, d, name, cached = cached)


function FreeMod_dec(F::FreeMod, d::Vector{FinGenAbGroupElem})
  return FreeMod_dec{elem_type(base_ring(F))}(F, d)
end


function AbstractAlgebra.extra_name(F::FreeMod_dec)
  t = get_attribute(F, :twist)
  if t !== nothing
    n = AbstractAlgebra.get_name(t[1])
    if n !== nothing
      return "$n($(t[2]))"
    end
  end
  if length(Set(F.d)) == 1
    n = AbstractAlgebra.get_name(base_ring(forget_decoration(F)))
    if n !== nothing
      return "$n^$(ngens(F))($(-F.d[1]))"
    end
  end
  return nothing
end

function show(io::IO, F::FreeMod_dec)
  @show_name(io, F)
  @show_special(io, F)

  io = terse(io)
  io = pretty(io)
  print(io, "Decorated free module of rank $(rank(F)) over ")
  print(IOContext(io, :compact => true), Lowercase(), base_ring(F))
end

# Generic specialized show methods (formerly in Hecke)

function Hecke.show_hom(io::IO, G)
  D = get_attribute(G, :hom)
  D === nothing && error("only for hom")
  print(io, "hom of ")
  print(IOContext(io, :compact => true), D)
end

function Hecke.show_direct_product(io::IO, G)
  D = get_attribute(G, :direct_product)
  D === nothing && error("only for direct products")
  print(io, "direct product of ")
  show(IOContext(io, :compact => true), D)
end

function Hecke.show_direct_sum(io::IO, G)
  D = get_attribute(G, :direct_product)
  D === nothing && error("only for direct sums")
  print(io, "direct sum of ")
  show(IOContext(io, :compact => true), D)
end

function Hecke.show_tensor_product(io::IO, G)
  D = get_attribute(G, :tensor_product)
  D === nothing && error("only for tensor products")
  print(io, "tensor product of ")
  show(IOContext(io, :compact => true), D)
end

function forget_decoration(F::FreeMod_dec)
  return F.F
end

@doc raw"""
    base_ring(F::FreeMod_dec)

Return the underlying ring of `F`.
"""
base_ring(F::FreeMod_dec) = forget_decoration(F).R

base_ring_type(::Type{FreeMod_dec{T}}) where {T} = base_ring_type(FreeMod{T})

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
  res = Dict{FinGenAbGroupElem, FreeModElem_dec}()
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
    homogeneous_component(a::FreeModElem_dec, g::FinGenAbGroupElem)

Return the homogeneous component of `a` which has degree `g`.
"""
function homogeneous_component(a::FreeModElem_dec, g::FinGenAbGroupElem)
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
  return F, MapFromFunc(Hecke.TupleParent(Tuple([zero(g) for g = G])), F, pure, inv_pure)
end


###############################################################################
# FreeModuleHom_dec constructors
###############################################################################

FreeModuleHom_dec(F::FreeMod_dec{T}, G::SparseFPModule_dec, a::Vector) where {T} = FreeModuleHom_dec{T}(F, G, a)

FreeModuleHom_dec(F::FreeMod_dec{T}, G::SparseFPModule_dec, mat::MatElem{T}) where {T} = FreeModuleHom{T}(F, G, mat)

function forget_decoration_on_morphism(f::FreeModuleHom_dec)
  return f.f
end

function forget_decoration(f::FreeModuleHom_dec)
  F = forget_decoration(domain(f))
  G = forget_decoration(codomain(f))
  return hom(F, G, [forget_decoration(f(v)) for v in gens(domain(f))]; check=false)
end

function matrix(a::FreeModuleHom_dec)
  return matrix(forget_decoration_on_morphism(a))
end

(h::FreeModuleHom_dec)(a::FreeModElem_dec) = image(h, a)

hom(F::FreeMod_dec{T}, G::SparseFPModule_dec{T}, V::Vector{<:FreeModElem_dec}) where T = FreeModuleHom_dec(F, G, V) 
hom(F::FreeMod_dec{T}, G::SparseFPModule_dec{T}, A::MatElem{T}) where T = FreeModuleHom_dec(F, G, A)


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
    minimal_betti_table(M::SparseFPModule{T}) where {T<:MPolyDecRingElem}
    minimal_betti_table(A::MPolyQuoRing{T}) where {T<:MPolyDecRingElem}
    minimal_betti_table(I::MPolyIdeal{T}) where {T<:MPolyDecRingElem} 

Given a finitely presented graded module `M` over a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the Betti Table
of the minimal free resolution of `M`. Similarly for `A` and `I`.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> I = ideal(R, [w^2-x*z, w*x-y*z, x^2-w*y, x*y-z^2, y^2-w*z]);

julia> A, _ = quo(R, I)
(Quotient of multivariate polynomial ring by ideal (w^2 - x*z, w*x - y*z, -w*y + x^2, x*y - z^2, -w*z + y^2), Map: R -> A)

julia> minimal_betti_table(A)
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  5  5  -
     2: -  -  -  1
------------------
 total: 1  5  5  1
```
"""
function minimal_betti_table(M::SparseFPModule{T}; check::Bool=true) where {T<:MPolyDecRingElem}
 error("Not implemented for the given type")
end

function minimal_betti_table(M::SubquoModule{T}; check::Bool=true) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(M); check)
end

function minimal_betti_table(F::FreeMod{T}; check::Bool=true) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(M); check)
end

function minimal_betti_table(A::MPolyQuoRing{T}; check::Bool=true) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(A); check)
end

function minimal_betti_table(I::MPolyIdeal{T}; check::Bool=true) where {T<:MPolyDecRingElem}
  return minimal_betti_table(free_resolution(I); check)
end



@doc raw"""
    minimal_betti_table(F::FreeResolution{T}; check::Bool=true) where {T<:SparseFPModule}

Given a graded free resolution `F` over a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the
Betti table of the minimal free resolution arising from `F`.

!!! note
    The algorithm proceeds without actually minimizing the resolution.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> I = ideal(R, [w^2-x*z, w*x-y*z, x^2-w*y, x*y-z^2, y^2-w*z]);

julia> A, _ = quo(R, I)
(Quotient of multivariate polynomial ring by ideal (w^2 - x*z, w*x - y*z, -w*y + x^2, x*y - z^2, -w*z + y^2), Map: R -> A)

julia> FA = free_resolution(A)
Free resolution of A
R^1 <---- R^5 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> betti_table(FA)
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  5  5  1
     2: -  -  1  1
------------------
 total: 1  5  6  2

julia> minimal_betti_table(FA)
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  5  5  -
     2: -  -  -  1
------------------
 total: 1  5  5  1
```
"""
function minimal_betti_table(res::FreeResolution{T}; check::Bool=true) where {T<:SparseFPModule}
  @assert is_standard_graded(base_ring(res)) "resolution must be defined over a standard graded ring"
  @assert is_graded(res) "resolution must be graded"
  C = complex(res)
  @assert is_complete(res) "resolution must be complete"
  rng = range(C)
  # The following needs the resolution to be complete to be true
  res_length = first(rng)-1
  offsets = Dict{FinGenAbGroupElem, Int}()
  betti_hash_table = Dict{Tuple{Int, Any}, Int}()
  for i in 1:res_length+1
    phi = map(C, i)
    F = domain(phi)
    G = codomain(phi)
    dom_degs = unique!([degree(g; check) for g in gens(F)])
    cod_degs = unique!([degree(g; check) for g in gens(G)])
    for d in cod_degs
      d::FinGenAbGroupElem
      if d in dom_degs
        _, _, sub_mat = _constant_sub_matrix(phi, d; check)
        r = rank(sub_mat)
        c = ncols(sub_mat) - r - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
        offsets[d] = r
      else
        c = length(_indices_of_generators_of_degree(G, d; check)) - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
      end
    end
  end
  return BettiTable(betti_hash_table)
end

function hash_table(B::BettiTable) 
  return B.B
end

# TODO: Where is this called??? Adjust the use of `check` there!
function generators_of_degree(
    C::FreeResolution{T},
    i::Int,
    d::FinGenAbGroupElem;
    check::Bool=true
  ) where {T<:SparseFPModule}
  F = C[i]
  return [g for g in gens(F) if degree(g) == d]
end

function _indices_of_generators_of_degree(
    F::FreeMod{T}, d::FinGenAbGroupElem; check::Bool=true
  ) where {T<:Union{MPolyDecRingElem{<:FieldElem}, 
                    MPolyQuoRingElem{<:MPolyDecRingElem{<:FieldElem}}}}
  return Int[i for (i, g) in enumerate(gens(F)) if degree(g; check) == d]
end

function _constant_sub_matrix(
    phi::FreeModuleHom{T, T},
    d::FinGenAbGroupElem;
    check::Bool=true
  ) where {RET<:Union{MPolyDecRingElem{<:FieldElem}, 
                      MPolyQuoRingElem{<:MPolyDecRingElem{<:FieldElem}}
                     }, T<:FreeMod{RET}}
  S = base_ring(domain(phi))::Union{MPolyDecRing, MPolyQuoRing{<:MPolyDecRingElem}}
  kk = coefficient_ring(S)::Field
  F = domain(phi)
  G = codomain(phi)
  ind_dom = _indices_of_generators_of_degree(F, d; check)
  ind_cod = _indices_of_generators_of_degree(G, d; check)
  m = length(ind_dom)
  n = length(ind_cod)
  result = zero_matrix(kk, m, n)
  img_gens = images_of_generators(phi)
  for (i, l) in enumerate(ind_dom)
    v = coordinates(img_gens[l])
    for (j, k) in enumerate(ind_cod)
      success, c = _has_index(v, k)
      !success && continue
      result[i, j] = _constant_coeff(c)
    end
  end
  return ind_dom, ind_cod, result
end

_constant_coeff(f::MPolyDecRingElem) = first(coefficients(f))
_constant_coeff(f::MPolyQuoRingElem) = first(coefficients(lift(f)))

# TODO: This will be provided soon from different sources.
function complex(F::FreeResolution) 
  return F.C
end

function base_ring(res::FreeResolution{T}) where {T<:SparseFPModule}
  return base_ring(res[-1])
end
                                                                                                                                    
#############truncation#############

@doc raw"""
    truncate(I::SparseFPModule, g::FinGenAbGroupElem, task::Symbol=:with_morphism)

Given a finitely presented graded module `M` over a $\mathbb Z$-graded multivariate 
polynomial ring with positive weights, return the truncation of `M` at degree `g`.

Put more precisely, return the truncation as an object of type `SubquoModule`. 

Additionally, if `N` denotes this object,

- return the inclusion map `N` $\to$ `M` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `M` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.

    truncate(M::SparseFPModule, d::Int, task::Symbol = :with_morphism)

Given a module `M` as above, and given an integer `d`, convert `d` into an element `g`
of the grading group of `base_ring(I)` and proceed as above.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(R, 1)
Graded free module R^1([0]) of rank 1 over R

julia> V = [x*F[1]; y^4*F[1]; z^5*F[1]];

julia> M, _ = quo(F, V);

julia> M[1]
e[1]

julia> MT = truncate(M, 3);

julia> MT[1]
Graded subquotient of graded submodule of F with 10 generators
  1: z^3*e[1]
  2: y*z^2*e[1]
  3: y^2*z*e[1]
  4: y^3*e[1]
  5: x*z^2*e[1]
  6: x*y*z*e[1]
  7: x*y^2*e[1]
  8: x^2*z*e[1]
  9: x^2*y*e[1]
  10: x^3*e[1]
by graded submodule of F with 3 generators
  1: x*e[1]
  2: y^4*e[1]
  3: z^5*e[1]
```
"""
function truncate(I::SparseFPModule, g::FinGenAbGroupElem, task::Symbol=:with_morphism)
  return truncate(I, Int(g[1]), task)
end

function truncate(I::SparseFPModule, d::Int, task::Symbol=:with_morphism; check::Bool=true)
  @req I isa FreeMod || I isa SubquoModule "Not implemented for the given type"
  R = base_ring(I)
  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  @req is_z_graded(R) "The base ring must be ZZ-graded"
  W = R.d
  W = [Int(W[i][1]) for i = 1:ngens(R)]
  @req minimum(W) > 0 "The weights must be positive"
  if is_zero(I)
     return _return_wrt_task((I, identity_map(I)), task)
  end
  dmin = minimum(degree(Int, x; check) for x in gens(I))
  if  d <= dmin
     return _return_wrt_task((I, identity_map(I)), task)
  end
  V = sort(gens(I); by=a -> degree(Int, a; check))
  RES = elem_type(I)[]
  s = dmin
  B = monomial_basis(R, d-s)
  for i = 1:length(V)
    if degree(Int, V[i]; check) < d
      if degree(Int, V[i]; check) > s
        s = degree(Int, V[i]; check)
        B = monomial_basis(R, d-s)
      end
      append!(RES, [x*V[i] for x in B])
    else
      push!(RES, V[i])
    end
  end
  return _return_wrt_task(sub(I, RES), task)
end

function _return_wrt_task(result, task)
  if task == :with_morphism
    return result
  elseif task == :cache_morphism
    register_morphism!(result[2])
    return result[2]
  elseif task == :only_morphism
    return result[2]
  elseif task == :none
    return result[1]
  else
    error("task not recognized")
  end
end
  


##################regularity#######################

@doc raw"""
    cm_regularity(M::SparseFPModule; check::Bool=true)

Given a finitely presented graded module `M` over a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the
Castelnuovo-Mumford regularity of `M`.

    cm_regularity(I::MPolyIdeal)

Given a (homogeneous) ideal `I` in a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the
Castelnuovo-Mumford regularity of `I`.

# Examples
```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(R, 1);

julia> M, _ = quo(F, [x^2*F[1], y^2*F[1], z^2*F[1]])
(Graded subquotient of graded submodule of F with 1 generator
  1: e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^2*e[1]
  3: z^2*e[1], Hom: F -> M)

julia> cm_regularity(M)
3

julia> minimal_betti_table(M)
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  3  -  -
     2: -  -  3  -
     3: -  -  -  1
------------------
 total: 1  3  3  1 
```
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> I = ideal(R, [-x*z+y^2, x*y-w*z, x^2-w*y]);

julia> cm_regularity(I)
2

julia> A, _ = quo(R, I);

julia> minimal_betti_table(A)
degree: 0  1  2
---------------
     0: 1  -  -
     1: -  3  2
---------------
 total: 1  3  2 
```
"""
function cm_regularity(M::SparseFPModule; check::Bool=true)
 error("Not implemented for the given type")
end

function cm_regularity(M::FreeMod{T}; check::Bool=true) where {T<:MPolyRingElem{<:FieldElem}}
   R = base_ring(M)
   @check is_standard_graded(R) "The base ring is not standard ZZ-graded"
   return 0
end

function cm_regularity(M::SubquoModule{T}; check::Bool=true) where {T<:MPolyRingElem{<:FieldElem}}
   R = base_ring(M)
   @check coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
   @check is_standard_graded(R) "The base ring is not standard ZZ-graded"
   B = minimal_betti_table(M; check)
   S = as_dictionary(B)
   V = [x[2][1] - x[1] for x in keys(S)] 
  return maximum(V)
end

#####affine algebras as modules#####

@doc raw"""
    quotient_ring_as_module(A::MPolyQuoRing)

Return `A` considered as an object of type `SubquoModule`.

    quotient_ring_as_module(I::Union{MPolyIdeal, MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal})

As above, where `A` is the quotient of `base_ring(I)` modulo `I`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> IR = ideal(R, [x^2, y^3]);

julia> quotient_ring_as_module(IR)
Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 2 generators
  1: x^2*e[1]
  2: y^3*e[1]

julia> base_ring(ans)
Multivariate polynomial ring in 2 variables x, y
  over rational field

julia> A, _ = quo(R, ideal(R,[x*y]));

julia> AI = ideal(A, [x^2, y^3]);

julia> quotient_ring_as_module(AI)
Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 2 generators
  1: x^2*e[1]
  2: y^3*e[1]

julia> base_ring(ans)
Quotient
  of multivariate polynomial ring in 2 variables x, y
    over rational field
  by ideal (x*y)

```
```jldoctest
julia> S, (x, y) = graded_polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(S, [x^2, y^3])
Ideal generated by
  x^2
  y^3

julia> quotient_ring_as_module(I)
Graded subquotient of graded submodule of S^1 with 1 generator
  1: e[1]
by graded submodule of S^1 with 2 generators
  1: x^2*e[1]
  2: y^3*e[1]

```
"""
function quotient_ring_as_module(A::MPolyQuoRing)
  return quotient_ring_as_module(modulus(A))
end

function quotient_ring_as_module(I::Union{MPolyIdeal, MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal})
  R = base_ring(I)
  F = is_graded(R) ? graded_free_module(R, 1) : free_module(R, 1)
  e1 = F[1]
  return quo_object(F, [x * e1 for x = gens(I)])
end

#####ideals as modules#####

@doc raw"""
    ideal_as_module(I::Union{MPolyIdeal, MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal})

Return `I` considered as an object of type `SubquoModule`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x^2, y^3])
Ideal generated by
  x^2
  y^3

julia> ideal_as_module(I)
Submodule with 2 generators
  1: x^2*e[1]
  2: y^3*e[1]
represented as subquotient with no relations
```
```jldoctest
julia> S, (x, y) = graded_polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(S, [x^2, y^3])
Ideal generated by
  x^2
  y^3

julia> ideal_as_module(I)
Graded submodule of S^1 with 2 generators
  1: x^2*e[1]
  2: y^3*e[1]
represented as subquotient with no relations
```
"""
function ideal_as_module(I::Union{MPolyIdeal, MPolyQuoIdeal, MPolyLocalizedIdeal, MPolyQuoLocalizedIdeal})
  R = base_ring(I)
  F = is_graded(R) ? graded_free_module(R, 1) : free_module(R, 1)
  e1 = F[1]
  return sub_object(F, [x * e1 for x = gens(I)])
end



##########################################################################
##### Twists
##########################################################################

@doc raw"""
    twist(M::SparseFPModule{T}, g::FinGenAbGroupElem) where {T<:MPolyDecRingElem}

Return the twisted module `M(g)`.

    twist(M::SparseFPModule{T}, W::Vector{<:IntegerUnion}) where {T<:MPolyDecRingElem}

Given a module `M` over a $\mathbb Z^m$-graded polynomial ring and a vector `W` of $m$ integers, 
convert `W` into an element `g` of the grading group of the ring and proceed as above.

    twist(M::SparseFPModule{T}, d::IntegerUnion) where {T<:MPolyDecRingElem}

Given a module `M` over a $\mathbb Z$-graded polynomial ring and an integer `d`, 
convert `d` into an element `g` of the grading group of the ring and proceed as above.

# Examples
```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [zero(R)])
Ideal generated by
  0

julia> M = quotient_ring_as_module(I)
Graded submodule of R^1 with 1 generator
  1: e[1]
represented as subquotient with no relations

julia> degree(gen(M, 1))
[0]

julia> N = twist(M, 2)
Graded submodule of R^1 with 1 generator
  1: e[1]
represented as subquotient with no relations

julia> degree(gen(N, 1))
[-2]

```
"""
function twist(M::SparseFPModule{T}, g::FinGenAbGroupElem) where {T<:Union{MPolyDecRingElem, MPolyQuoRingElem{<:MPolyDecRingElem}}}
 error("Not implemented for the given type")
end

function twist(M::SubquoModule{T}, g::FinGenAbGroupElem) where {T<:Union{MPolyDecRingElem, MPolyQuoRingElem{<:MPolyDecRingElem}}}
 R = base_ring(M)
 @req parent(g) == grading_group(R) "Group element not contained in grading group of base ring"
 F = ambient_free_module(M)
 FN = twist(F, g)
 GN = free_module(R, ngens(M))
 HN = free_module(R, length(relations(M)))
 a = hom(GN, F, ambient_representatives_generators(M); check=false)
 b = hom(HN, F, relations(M); check=false)
 A = matrix(a)
 B = matrix(b)
 N = subquotient(FN, A, B)
 return N
end

function twist(F::FreeMod{T}, g::FinGenAbGroupElem) where {T<:Union{MPolyDecRingElem, MPolyQuoRingElem{<:MPolyDecRingElem}}}
 R = base_ring(F)
 @req parent(g) == grading_group(R) "Group element not contained in grading group of base ring"
 W = [x-g for x in F.d]
 G = graded_free_module(R, rank(F))
 G.d = W
 return G
end

function twist(M::SparseFPModule{T}, W::Vector{<:IntegerUnion}) where {T<:MPolyDecRingElem}
  R = base_ring(M)
  @assert is_zm_graded(R)
  return twist(M, grading_group(R)(W))
end

function twist(M::SparseFPModule{T}, d::IntegerUnion) where {T<:MPolyDecRingElem}
  R = base_ring(M)
  @assert is_z_graded(R)
  return twist(M, grading_group(R)([d]))
end

# TODO: implement further base changes. But for this we need more functionality 
# for potential change of grading groups. See the open pr #2677.
function change_base_ring(
    f::RingFlattening{DomType, CodType}, F::FreeMod
  ) where {DomType<:MPolyDecRing, CodType<:Ring}
  domain(f) === base_ring(F) || error("ring map not compatible with the module")
  S = codomain(f)
  r = ngens(F)
  FS = grade(FreeMod(S, F.S), degree.(gens(F)))
  map = hom(F, FS, gens(FS), f; check=false)
  return FS, map
end

function change_base_ring(f::RingFlattening{DomType, CodType}, M::SubquoModule) where {DomType<:MPolyDecRing, CodType<:Ring}
  domain(f) == base_ring(M) || error("ring map not compatible with the module")
  S = codomain(f)
  F = ambient_free_module(M)
  R = base_ring(M)
  FS, mapF, iso_inv = f(F) # Makes sure ambient modules are preserved due to caching in flattenings
  g = ambient_representatives_generators(M)
  rels = relations(M)
  MS = SubquoModule(FS, mapF.(g), mapF.(rels))
  map = SubQuoHom(M, MS, gens(MS), f; check=false)
  return MS, map
end

function _regularity_bound(M::SubquoModule)
  @assert is_graded(M) "module must be graded"
  S = base_ring(M)
  @assert is_z_graded(S) "base ring must be ZZ-graded"
  @assert all(x->degree(Int, x; check=false) >= 0, gens(S)) "base ring variables must be non-negatively graded"
  res = free_resolution(M)
  result = maximum((x->degree(Int, x; check=false)).(gens(res[0])))
  for i in 0:first(chain_range(res))
    result = maximum(push!((x->degree(Int, x; check=false)).(gens(res[i])), result))
  end
  return result
end


###############################################################################
# Random elements
###############################################################################

function rand_homogeneous(R::MPolyRing, degree::Int)
  K = base_ring(R)
  if !is_standard_graded(R)
      throw(ArgumentError("Base ring is not standard graded"))
  end
  if !is_finite(K)
      throw(ArgumentError("Base ring is not finite"))
  end
  n = nvars(R)
  M = MPolyBuildCtx(R)
  for p in weak_compositions(degree, n)
    push_term!(M, rand(K), data(p))
  end
  return finish(M)
end


function rand_homogeneous(V::SparseFPModule, d::Int)
  R = base_ring(V)
  random_element = zero(V)
  for gen in gens(V)
      gen_degree = (degree(gen).coeff)[1,1]
      if gen_degree <= d
          rand_poly = rand_homogeneous(R, Int(d - gen_degree))
          random_element += rand_poly*gen
      end
  end
  return random_element
end
