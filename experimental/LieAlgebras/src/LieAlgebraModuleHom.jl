@attributes mutable struct LieAlgebraModuleHom{T1<:LieAlgebraModule,T2<:LieAlgebraModule} <:
                           Map{T1,T2,Hecke.HeckeMap,LieAlgebraModuleHom}
  header::MapHeader{T1,T2}
  matrix::MatElem

  inverse_isomorphism::LieAlgebraModuleHom{T2,T1}

  function LieAlgebraModuleHom(
    V1::LieAlgebraModule,
    V2::LieAlgebraModule,
    imgs::Vector{<:LieAlgebraModuleElem};
    check::Bool=true,
  )
    @req base_lie_algebra(V1) === base_lie_algebra(V2) "Lie algebras must be the same" # for now at least
    @req all(x -> parent(x) === V2, imgs) "Images must lie in the codomain"
    @req length(imgs) == dim(V1) "Number of images must match dimension of domain"

    mat = zero_matrix(coefficient_ring(V2), dim(V1), dim(V2))
    for (i, img) in enumerate(imgs)
      mat[i, :] = _matrix(img)
    end
    return LieAlgebraModuleHom(V1, V2, mat; check)
  end

  function LieAlgebraModuleHom(
    V1::LieAlgebraModule, V2::LieAlgebraModule, mat::MatElem; check::Bool=true
  )
    @req base_lie_algebra(V1) === base_lie_algebra(V2) "Lie algebras must be the same" # for now at least
    @req size(mat) == (dim(V1), dim(V2)) "Matrix size must match dimensions of domain and codomain"
    h = new{typeof(V1),typeof(V2)}()
    h.matrix = mat::dense_matrix_type(coefficient_ring(V2))
    h.header = MapHeader(V1, V2)
    if check
      for x in basis(base_lie_algebra(V1)), v in basis(V1)
        @req x * h(v) == h(x * v) "Not a homomorphism"
      end
    end
    return h
  end
end

###############################################################################
#
#   Basic properties
#
###############################################################################

@doc raw"""
    matrix(h::LieAlgebraModuleHom) -> MatElem
    
Return the transformation matrix of `h` w.r.t. the bases of the domain and codomain.

Note: The matrix operates on the coefficient vectors from the right.
"""
function matrix(
  h::LieAlgebraModuleHom{<:LieAlgebraModule,<:LieAlgebraModule{C2}}
) where {C2<:RingElement}
  return (h.matrix)::dense_matrix_type(C2)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", h::LieAlgebraModuleHom)
  io = pretty(io)
  println(io, LowercaseOff(), "Lie algebra module morphism")
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(h))
  print(io, "to   ", Lowercase(), codomain(h))
  print(io, Dedent())
end

function Base.show(io::IO, h::LieAlgebraModuleHom)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, LowercaseOff(), "Lie algebra module morphism")
  else
    print(io, LowercaseOff(), "Lie algebra module morphism: ")
    print(io, Lowercase(), domain(h), " -> ", Lowercase(), codomain(h))
  end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(h1::LieAlgebraModuleHom, h2::LieAlgebraModuleHom)
  return domain(h1) === domain(h2) &&
         codomain(h1) === codomain(h2) &&
         matrix(h1) == matrix(h2)
end

function Base.hash(f::LieAlgebraModuleHom, h::UInt)
  b = 0x7b887214521c25b7 % UInt
  h = hash(objectid(domain(f)), h)
  h = hash(objectid(codomain(f)), h)
  h = hash(matrix(f), h)
  return xor(h, b)
end

###############################################################################
#
#   Image and kernel
#
###############################################################################

@doc raw"""
    image(h::LieAlgebraModuleHom, v::LieAlgebraModuleElem) -> LieAlgebraModuleElem

Return the image of `v` under `h`.
"""
function image(
  h::LieAlgebraModuleHom{T1,T2}, v::LieAlgebraModuleElem
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule}
  @req parent(v) === domain(h) "Domain mismatch"
  return codomain(h)(_matrix(v) * matrix(h))
end

# TODO: image and kernel, once submodules are implemented

###############################################################################
#
#   Map operations
#
###############################################################################

@doc raw"""
    compose(f::LieAlgebraModuleHom, g::LieAlgebraModuleHom) -> LieAlgebraModuleHom

Return the composition of `f` and `g`, i.e. the homomorphism `h` such that
`h(x) = g(f(x))` for all `x` in the domain of `f`.
The codomain of `f` must be identical to the domain of `g`.
"""
function compose(
  f::LieAlgebraModuleHom{T1,T2}, g::LieAlgebraModuleHom{T2,T3}
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule,T3<:LieAlgebraModule}
  @req codomain(f) === domain(g) "Composition: Maps are not compatible"
  h = LieAlgebraModuleHom(domain(f), codomain(g), matrix(f) * matrix(g); check=false)
  if isdefined(f, :inverse_isomorphism) && isdefined(g, :inverse_isomorphism)
    h.inverse_isomorphism = LieAlgebraModuleHom(
      codomain(g),
      domain(f),
      matrix(g.inverse_isomorphism) * matrix(f.inverse_isomorphism);
      check=false,
    )
    h.inverse_isomorphism.inverse_isomorphism = h
  end
  return h
end

@doc raw"""
    inv(h::LieAlgebraModuleHom) -> LieAlgebraModuleHom

Return the inverse of `h`.
Requires `h` to be an isomorphism.
"""
function inv(h::LieAlgebraModuleHom)
  @req is_isomorphism(h) "Homomorphism must be invertible"
  return h.inverse_isomorphism
end

@doc raw"""
    is_isomorphism(h::LieAlgebraModuleHom) -> Bool

Return `true` if `h` is an isomorphism.
This function tries to invert the transformation matrix of `h` and caches the result.
The inverse isomorphism can be cheaply accessed via `inv(h)` after calling this function.
"""
@attr Bool function is_isomorphism(h::LieAlgebraModuleHom)
  isdefined(h, :inverse_isomorphism) && return true
  fl, invmat = is_invertible_with_inverse(matrix(h))
  fl || return false
  h.inverse_isomorphism = LieAlgebraModuleHom(codomain(h), domain(h), invmat; check=false)
  h.inverse_isomorphism.inverse_isomorphism = h
  return true
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    hom(V1::LieAlgebraModule, V2::LieAlgebraModule, imgs::Vector{<:LieAlgebraModuleElem}; check::Bool=true) -> LieAlgebraModuleHom

Construct the homomorphism from `V1` to `V2` by sending the `i`-th basis element of `V1` 
to `imgs[i]` and extending linearly.
All elements of `imgs` must lie in `V2`.
Currently, `V1` and `V2` must be modules over the same Lie algebra.

By setting `check=false`, the linear map is not checked to be compatible with the module action.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 2);

julia> V1 = standard_module(L);

julia> V3 = trivial_module(L, 3);

julia> V2 = direct_sum(V1, V3);

julia> h = hom(V1, V2, [V2([v, zero(V3)]) for v in basis(V1)])
Lie algebra module morphism
  from standard module of dimension 2 over sl_2
  to   direct sum module of dimension 5 over sl_2

julia> [(v, h(v)) for v in basis(V1)]
2-element Vector{Tuple{LieAlgebraModuleElem{QQFieldElem}, LieAlgebraModuleElem{QQFieldElem}}}:
 (v_1, v_1^(1))
 (v_2, v_2^(1))
```
"""
function hom(
  V1::LieAlgebraModule{C},
  V2::LieAlgebraModule{C},
  imgs::Vector{<:LieAlgebraModuleElem{C}};
  check::Bool=true,
) where {C<:RingElement}
  return LieAlgebraModuleHom(V1, V2, imgs; check)
end

@doc raw"""
    hom(V1::LieAlgebraModule, V2::LieAlgebraModule, mat::MatElem; check::Bool=true) -> LieAlgebraModuleHom

Construct the homomorphism from `V1` to `V2` by acting with the matrix `mat`
from the right on the coefficient vector w.r.t. the basis of `V1`.
`mat` must be a matrix of size `dim(V1) \times dim(V2)` over `coefficient_ring(V2)`.
Currently, `V1` and `V2` must be modules over the same Lie algebra.

By setting `check=false`, the linear map is not checked to be compatible with the module action.

# Examples
```jldoctest
julia> L = general_linear_lie_algebra(QQ, 3);

julia> V1 = standard_module(L);

julia> V2 = trivial_module(L);

julia> h = hom(V1, V2, matrix(QQ, 3, 1, [0, 0, 0]))
Lie algebra module morphism
  from standard module of dimension 3 over gl_3
  to   abstract Lie algebra module of dimension 1 over gl_3
  
julia> [(v, h(v)) for v in basis(V1)]
3-element Vector{Tuple{LieAlgebraModuleElem{QQFieldElem}, LieAlgebraModuleElem{QQFieldElem}}}:
 (v_1, 0)
 (v_2, 0)
 (v_3, 0)
```
"""
function hom(
  V1::LieAlgebraModule{C}, V2::LieAlgebraModule{C}, mat::MatElem{C}; check::Bool=true
) where {C<:RingElement}
  return LieAlgebraModuleHom(V1, V2, mat; check)
end

@doc raw"""
    identity_map(V::LieAlgebraModule) -> LieAlgebraModuleHom

Construct the identity map on `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = standard_module(L)
Standard module
  of dimension 3
over special linear Lie algebra of degree 3 over QQ

julia> identity_map(V)
Lie algebra module morphism
  from standard module of dimension 3 over sl_3
  to   standard module of dimension 3 over sl_3
```
"""
function identity_map(V::LieAlgebraModule)
  return hom(V, V, basis(V); check=false)
end
