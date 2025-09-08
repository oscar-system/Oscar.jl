###############################################################################
#
#   Basic properties
#
###############################################################################

@doc raw"""
    matrix(h::LieAlgebraHom) -> MatElem
    
Return the transformation matrix of `h` w.r.t. the bases of the domain and codomain.

Note: The matrix operates on the coefficient vectors from the right.
"""
function matrix(h::LieAlgebraHom{<:LieAlgebra,<:LieAlgebra,MatT}) where {MatT<:MatElem}
  return (h.matrix)::MatT
end

@doc raw"""
    is_welldefined(h::LieAlgebraHom) -> Bool

Return `true` if `h` is a well-defined homomorphism of Lie algebras.
This function is used internally when calling `hom` with `check=true`.
"""
function is_welldefined(h::LieAlgebraHom)
  L1 = domain(h)
  for x1 in basis(L1), x2 in basis(L1)
    h(x1) * h(x2) == h(x1 * x2) || return false
  end
  return true
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", h::LieAlgebraHom)
  @show_name(io, h)
  @show_special(io, mime, h)
  io = pretty(io)
  println(io, LowercaseOff(), "Lie algebra morphism")
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(h))
  print(io, "to ", Lowercase(), codomain(h))
  print(io, Dedent())
end

function Base.show(io::IO, h::LieAlgebraHom)
  @show_name(io, h)
  @show_special(io, h)
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Lie algebra morphism")
  else
    print(io, LowercaseOff(), "Lie algebra morphism: ")
    print(io, Lowercase(), domain(h), " -> ", Lowercase(), codomain(h))
  end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(h1::LieAlgebraHom, h2::LieAlgebraHom)
  return domain(h1) === domain(h2) &&
         codomain(h1) === codomain(h2) &&
         matrix(h1) == matrix(h2)
end

function Base.hash(f::LieAlgebraHom, h::UInt)
  b = 0xc023adc432e006be % UInt
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
    image(h::LieAlgebraHom, x::LieAlgebraElem) -> LieAlgebraElem

Return the image of `x` under `h`.
"""
function image(h::LieAlgebraHom, x::LieAlgebraElem)
  @req parent(x) === domain(h) "Domain mismatch"
  return codomain(h)(_matrix(x) * matrix(h))
end

@doc raw"""
    image(h::LieAlgebraHom) -> LieSubalgebra

Return the image of `h` as a Lie subalgebra of the codomain.
"""
function image(h::LieAlgebraHom)
  return sub(codomain(h), [image(h, x) for x in basis(domain(h))])
end

@doc raw"""
    image(h::LieAlgebraHom, I::LieAlgebraIdeal) -> LieSubalgebra

Return the image of `I` under `h` as a Lie subalgebra of the codomain.
"""
function image(h::LieAlgebraHom, I::LieAlgebraIdeal)
  @req base_lie_algebra(I) === domain(h) "Domain mismatch"
  return sub(codomain(h), [image(h, x) for x in basis(I)])
end

@doc raw"""
    image(h::LieAlgebraHom, S::LieSubalgebra) -> LieSubalgebra

Return the image of `S` under `h` as a Lie subalgebra of the codomain.
"""
function image(h::LieAlgebraHom, S::LieSubalgebra)
  @req base_lie_algebra(S) === domain(h) "Domain mismatch"
  return sub(codomain(h), [image(h, x) for x in basis(S)])
end

@doc raw"""
    kernel(h::LieAlgebraHom) -> LieAlgebraIdeal

Return the kernel of `h` as an ideal of the domain.
"""
function kernel(h::LieAlgebraHom)
  ker_b = kernel(matrix(h); side=:left)
  ker_dim = nrows(ker_b)

  return ideal(domain(h), [domain(h)(ker_b[i, :]) for i in 1:ker_dim])
end

###############################################################################
#
#   Map operations
#
###############################################################################

@doc raw"""
    compose(f::LieAlgebraHom, g::LieAlgebraHom) -> LieAlgebraHom

Return the composition of `f` and `g`, i.e. the homomorphism `h` such that
`h(x) = g(f(x))` for all `x` in the domain of `f`.
The codomain of `f` must be identical to the domain of `g`.
"""
function compose(
  f::LieAlgebraHom{T1,T2}, g::LieAlgebraHom{T2,T3}
) where {T1<:LieAlgebra,T2<:LieAlgebra,T3<:LieAlgebra}
  @req codomain(f) === domain(g) "Composition: Maps are not compatible"
  h = LieAlgebraHom(domain(f), codomain(g), matrix(f) * matrix(g); check=false)
  if isdefined(f, :inverse_isomorphism) && isdefined(g, :inverse_isomorphism)
    h.inverse_isomorphism = LieAlgebraHom(
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
    inv(h::LieAlgebraHom) -> LieAlgebraHom

Return the inverse of `h`.
Requires `h` to be an isomorphism.
"""
function inv(h::LieAlgebraHom)
  @req is_isomorphism(h) "Homomorphism must be invertible"
  return h.inverse_isomorphism
end

@doc raw"""
    is_isomorphism(h::LieAlgebraHom) -> Bool

Return `true` if `h` is an isomorphism.
This function tries to invert the transformation matrix of `h` and caches the result.
The inverse isomorphism can be cheaply accessed via `inv(h)` after calling this function.
"""
@attr Bool function is_isomorphism(h::LieAlgebraHom)
  isdefined(h, :inverse_isomorphism) && return true
  fl, invmat = is_invertible_with_inverse(matrix(h))
  fl || return false
  h.inverse_isomorphism = LieAlgebraHom(codomain(h), domain(h), invmat; check=false)
  h.inverse_isomorphism.inverse_isomorphism = h
  return true
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    hom(L1::LieAlgebra, L2::LieAlgebra, imgs::Vector{<:LieAlgebraElem}; check::Bool=true) -> LieAlgebraHom

Construct the homomorphism from `L1` to `L2` by sending the `i`-th basis element of `L1` 
to `imgs[i]` and extending linearly.
All elements of `imgs` must lie in `L2`.

By setting `check=false`, the linear map is not checked to be compatible with the Lie bracket.

# Examples
```jldoctest
julia> L1 = special_linear_lie_algebra(QQ, 2);

julia> L2 = special_linear_lie_algebra(QQ, 3);

julia> h = hom(L1, L2, [basis(L2, 1), basis(L2, 4), basis(L2, 7)]) # embed sl_2 into sl_3
Lie algebra morphism
  from special linear Lie algebra of degree 2 over QQ
  to special linear Lie algebra of degree 3 over QQ
  
julia> [(x, h(x)) for x in basis(L1)]
3-element Vector{Tuple{LinearLieAlgebraElem{QQFieldElem}, LinearLieAlgebraElem{QQFieldElem}}}:
 (e_1_2, e_1_2)
 (f_1_2, f_1_2)
 (h_1, h_1)
```
"""
function hom(
  L1::LieAlgebra{C}, L2::LieAlgebra{C}, imgs::Vector{<:LieAlgebraElem{C}}; check::Bool=true
) where {C<:FieldElem}
  return LieAlgebraHom(L1, L2, imgs; check)
end

@doc raw"""
    hom(L1::LieAlgebra, L2::LieAlgebra, mat::MatElem; check::Bool=true) -> LieAlgebraHom

Construct the homomorphism from `L1` to `L2` by acting with the matrix `mat`
from the right on the coefficient vector w.r.t. the basis of `L1`.
`mat` must be a matrix of size `dim(L1) \times dim(L2)` over `coefficient_ring(L2)`.

By setting `check=false`, the linear map is not checked to be compatible with the Lie bracket.

# Examples
```jldoctest
julia> L1 = special_linear_lie_algebra(QQ, 2);

julia> L2 = general_linear_lie_algebra(QQ, 2);

julia> h = hom(L1, L2, matrix(QQ, [0 1 0 0; 0 0 1 0; 1 0 0 -1]))
Lie algebra morphism
  from special linear Lie algebra of degree 2 over QQ
  to general linear Lie algebra of degree 2 over QQ

julia> [(x, h(x)) for x in basis(L1)]
3-element Vector{Tuple{LinearLieAlgebraElem{QQFieldElem}, LinearLieAlgebraElem{QQFieldElem}}}:
 (e_1_2, x_1_2)
 (f_1_2, x_2_1)
 (h_1, x_1_1 - x_2_2)
```
"""
function hom(
  L1::LieAlgebra{C}, L2::LieAlgebra{C}, mat::MatElem{C}; check::Bool=true
) where {C<:FieldElem}
  return LieAlgebraHom(L1, L2, mat; check)
end

@doc raw"""
    id_hom(L::LieAlgebra) -> LieAlgebraHom

Construct the identity morphism on `L`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3)
Special linear Lie algebra of degree 3
  of dimension 8
over rational field

julia> id_hom(L)
Lie algebra morphism
  from special linear Lie algebra of degree 3 over QQ
  to special linear Lie algebra of degree 3 over QQ
```
"""
function id_hom(L::LieAlgebra)
  return hom(L, L, basis(L); check=false)
end

@doc raw"""
    zero_map(L1::LieAlgebra, L2::LieAlgebra) -> LieAlgebraHom
    zero_map(L::LieAlgebra) -> LieAlgebraHom

Construct the zero map from `L1` to `L2` or from `L` to `L`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3)
Special linear Lie algebra of degree 3
  of dimension 8
over rational field

julia> zero_map(L)
Lie algebra morphism
  from special linear Lie algebra of degree 3 over QQ
  to special linear Lie algebra of degree 3 over QQ
```
"""
function zero_map(L1::LieAlgebra{C}, L2::LieAlgebra{C}) where {C<:FieldElem}
  return hom(L1, L2, zero_matrix(coefficient_ring(L2), dim(L1), dim(L2)); check=false)
end

function zero_map(L::LieAlgebra)
  return zero_map(L, L)
end
