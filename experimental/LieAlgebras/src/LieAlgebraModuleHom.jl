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
) where {C2<:FieldElem}
  return (h.matrix)::dense_matrix_type(C2)
end

@doc raw"""
    is_welldefined(h::LieAlgebraModuleHom) -> Bool

Return `true` if `h` is a well-defined homomorphism of Lie algebra modules.
This function is used internally when calling `hom` with `check=true`.
"""
function is_welldefined(h::LieAlgebraModuleHom)
  V1 = domain(h)
  for x in basis(base_lie_algebra(V1)), v in basis(V1)
    x * h(v) == h(x * v) || return false
  end
  return true
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", h::LieAlgebraModuleHom)
  @show_name(io, h)
  @show_special(io, mime, h)
  io = pretty(io)
  println(io, LowercaseOff(), "Lie algebra module morphism")
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(h))
  print(io, "to ", Lowercase(), codomain(h))
  print(io, Dedent())
end

function Base.show(io::IO, h::LieAlgebraModuleHom)
  @show_name(io, h)
  @show_special(io, h)
  io = pretty(io)
  if is_terse(io)
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
  from standard module of dimension 2 over L
  to direct sum module of dimension 5 over L

julia> [(v, h(v)) for v in basis(V1)]
2-element Vector{Tuple{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}, LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}}:
 (v_1, v_1^(1))
 (v_2, v_2^(1))
```
"""
function hom(
  V1::LieAlgebraModule{C},
  V2::LieAlgebraModule{C},
  imgs::Vector{<:LieAlgebraModuleElem{C}};
  check::Bool=true,
) where {C<:FieldElem}
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
  from standard module of dimension 3 over L
  to abstract Lie algebra module of dimension 1 over L

julia> [(v, h(v)) for v in basis(V1)]
3-element Vector{Tuple{LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}, LieAlgebraModuleElem{QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}}:
 (v_1, 0)
 (v_2, 0)
 (v_3, 0)
```
"""
function hom(
  V1::LieAlgebraModule{C}, V2::LieAlgebraModule{C}, mat::MatElem{C}; check::Bool=true
) where {C<:FieldElem}
  return LieAlgebraModuleHom(V1, V2, mat; check)
end

@doc raw"""
    id_hom(V::LieAlgebraModule) -> LieAlgebraModuleHom

Construct the identity map on `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = standard_module(L)
Standard module
  of dimension 3
over special linear Lie algebra of degree 3 over QQ

julia> id_hom(V)
Lie algebra module morphism
  from standard module of dimension 3 over L
  to standard module of dimension 3 over L
```
"""
function id_hom(V::LieAlgebraModule)
  return hom(V, V, basis(V); check=false)
end

@doc raw"""
    zero_map(V1::LieAlgebraModule, V2::LieAlgebraModule) -> LieAlgebraModuleHom
    zero_map(V::LieAlgebraModule) -> LieAlgebraModuleHom

Construct the zero map from `V1` to `V2` or from `V` to `V`.

# Examples
```jldoctest
julia> L = special_linear_lie_algebra(QQ, 3);

julia> V = standard_module(L)
Standard module
  of dimension 3
over special linear Lie algebra of degree 3 over QQ

julia> zero_map(V)
Lie algebra module morphism
  from standard module of dimension 3 over L
  to standard module of dimension 3 over L
```
"""
function zero_map(V1::LieAlgebraModule{C}, V2::LieAlgebraModule{C}) where {C<:FieldElem}
  return hom(V1, V2, zero_matrix(coefficient_ring(V2), dim(V1), dim(V2)); check=false)
end

function zero_map(V::LieAlgebraModule)
  return zero_map(V, V)
end

###############################################################################
#
#   Hom constructions
#
###############################################################################

function Base.:-(h::LieAlgebraModuleHom)
  return hom(domain(h), codomain(h), -matrix(h); check=false)
end

function Base.:+(
  h1::LieAlgebraModuleHom{T1,T2}, h2::LieAlgebraModuleHom{T1,T2}
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule}
  @req domain(h1) === domain(h2) "Maps must have the same domain"
  @req codomain(h1) === codomain(h2) "Maps must have the same codomain"
  return hom(domain(h1), codomain(h1), matrix(h1) + matrix(h2); check=false)
end

function Base.:-(
  h1::LieAlgebraModuleHom{T1,T2}, h2::LieAlgebraModuleHom{T1,T2}
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule}
  @req domain(h1) === domain(h2) "Maps must have the same domain"
  @req codomain(h1) === codomain(h2) "Maps must have the same codomain"
  return hom(domain(h1), codomain(h1), matrix(h1) - matrix(h2); check=false)
end

@doc raw"""
    canonical_injections(V::LieAlgebraModule) -> Vector{LieAlgebraModuleHom}

Return the canonical injections from all components into $V$
where $V$ has been constructed as $V_1 \oplus \cdot \oplus V_n$.
"""
function canonical_injections(V::LieAlgebraModule)
  fl, Vs = Oscar._is_direct_sum(V)
  @req fl "Module must be a direct sum"
  return [canonical_injection(V, i) for i in 1:length(Vs)]
end

@doc raw"""
    canonical_injection(V::LieAlgebraModule, i::Int) -> LieAlgebraModuleHom

Return the canonical injection $V_i \to V$
where $V$ has been constructed as $V_1 \oplus \cdot \oplus V_n$.
"""
function canonical_injection(V::LieAlgebraModule, i::Int)
  fl, Vs = Oscar._is_direct_sum(V)
  @req fl "Module must be a direct sum"
  @req 1 <= i <= length(Vs) "Index out of bound"
  j = sum(dim(Vs[l]) for l in 1:(i - 1); init=0)
  emb = hom(Vs[i], V, [basis(V, l + j) for l in 1:dim(Vs[i])]; check=false)
  return emb
end

@doc raw"""
    canonical_projections(V::LieAlgebraModule) -> Vector{LieAlgebraModuleHom}

Return the canonical projections from $V$ to all components
where $V$ has been constructed as $V_1 \oplus \cdot \oplus V_n$.
"""
function canonical_projections(V::LieAlgebraModule)
  fl, Vs = Oscar._is_direct_sum(V)
  @req fl "Module must be a direct sum"
  return [canonical_projection(V, i) for i in 1:length(Vs)]
end

@doc raw"""
    canonical_projection(V::LieAlgebraModule, i::Int) -> LieAlgebraModuleHom

Return the canonical projection $V \to V_i$
where $V$ has been constructed as $V_1 \oplus \cdot \oplus V_n$.
"""
function canonical_projection(V::LieAlgebraModule, i::Int)
  fl, Vs = Oscar._is_direct_sum(V)
  @req fl "Module must be a direct sum"
  @req 1 <= i <= length(Vs) "Index out of bound"
  j = sum(dim(Vs[l]) for l in 1:(i - 1); init=0)
  proj = hom(
    V,
    Vs[i],
    [
      [zero(Vs[i]) for l in 1:j]
      basis(Vs[i])
      [zero(Vs[i]) for l in (j + dim(Vs[i]) + 1):dim(V)]
    ];
    check=false,
  )
  return proj
end

@doc raw"""
    hom_direct_sum(V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, hs::Matrix{<:LieAlgebraModuleHom}) -> LieAlgebraModuleHom
    hom_direct_sum(V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, hs::Vector{<:LieAlgebraModuleHom}) -> LieAlgebraModuleHom

Given modules `V` and `W` which are direct sums with `r` respective `s` summands,
say $M = M_1 \oplus \cdots \oplus M_r$, $N = N_1 \oplus \cdots \oplus N_s$, and given a $r \times s$ matrix
`hs` of homomorphisms $h_{ij} : V_i \to W_j$, return the homomorphism
$V \to W$ with $ij$-components $h_{ij}$.

If `hs` is a vector, then it is interpreted as a diagonal matrix.
"""
function hom_direct_sum(
  V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, hs::Matrix{<:LieAlgebraModuleHom}
) where {C<:FieldElem}
  fl, Vs = Oscar._is_direct_sum(V)
  @req fl "First module must be a direct sum"
  fl, Ws = Oscar._is_direct_sum(W)
  @req fl "Second module must be a direct sum"
  @req length(Vs) == size(hs, 1) "Length mismatch"
  @req length(Ws) == size(hs, 2) "Length mismatch"
  @req all(
    domain(hs[i, j]) === Vs[i] && codomain(hs[i, j]) === Ws[j] for i in 1:size(hs, 1),
    j in 1:size(hs, 2)
  ) "Domain/codomain mismatch"

  Winjs = canonical_injections(W)
  Vprojs = canonical_projections(V)
  function map_basis(v)
    return sum(
      Winjs[j](sum(hs[i, j](Vprojs[i](v)) for i in 1:length(Vs); init=zero(Ws[j]))) for
      j in 1:length(Ws);
      init=zero(W),
    )
  end
  return hom(V, W, map(map_basis, basis(V)); check=false)
end

function hom_direct_sum(
  V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, hs::Vector{<:LieAlgebraModuleHom}
) where {C<:FieldElem}
  fl, Vs = Oscar._is_direct_sum(V)
  @req fl "First module must be a direct sum"
  fl, Ws = Oscar._is_direct_sum(W)
  @req fl "Second module must be a direct sum"
  @req length(Vs) == length(Ws) == length(hs) "Length mismatch"
  @req all(i -> domain(hs[i]) === Vs[i] && codomain(hs[i]) === Ws[i], 1:length(hs)) "Domain/codomain mismatch"

  return hom(V, W, diagonal_matrix(matrix.(hs)); check=false)
end

@doc raw"""
    hom_tensor(V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, hs::Vector{<:LieAlgebraModuleHom}) -> LieAlgebraModuleHom

Given modules `V` and `W` which are tensor products with the same number of factors,
say $V = V_1 \otimes \cdots \otimes V_r$, $W = W_1 \otimes \cdots \otimes W_r$,
and given a vector `hs` of homomorphisms $a_i : V_i \to W_i$, return
$a_1 \otimes \cdots \otimes a_r$.

This works for $r$th tensor powers as well.
"""
function hom_tensor(
  V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, hs::Vector{<:LieAlgebraModuleHom}
) where {C<:FieldElem} # TODO: cleanup after refactoring tensor_product
  if ((fl, Vs) = _is_tensor_product(V); fl)
    # nothing to do
  elseif ((fl, Vb, k) = _is_tensor_power(V); fl)
    Vs = [Vb for _ in 1:k]
  else
    throw(ArgumentError("First module must be a tensor product or power"))
  end
  if ((fl, Ws) = _is_tensor_product(W); fl)
    # nothing to do
  elseif ((fl, Wb, k) = _is_tensor_power(W); fl)
    Ws = [Wb for _ in 1:k]
  else
    throw(ArgumentError("Second module must be a tensor product or power"))
  end

  @req length(Vs) == length(Ws) == length(hs) "Length mismatch"
  @req all(i -> domain(hs[i]) === Vs[i] && codomain(hs[i]) === Ws[i], 1:length(hs)) "Domain/codomain mismatch"

  mat = reduce(
    kronecker_product,
    [matrix(hi) for hi in hs];
    init=identity_matrix(coefficient_ring(W), 1),
  )
  return hom(V, W, mat; check=false)
end

@doc raw"""
    hom(V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, h::LieAlgebraModuleHom) -> LieAlgebraModuleHom

Given modules `V` and `W` which are exterior/symmetric/tensor powers of the same kind with the same exponent,
say, e.g., $V = S^k V'$, $W = S^k W'$, and given a homomorphism $h : V' \to W'$, return
$S^k h: V \to W$ (analogous for other types of powers).
"""
function hom(
  V::LieAlgebraModule{C}, W::LieAlgebraModule{C}, h::LieAlgebraModuleHom
) where {C<:FieldElem}
  if ((fl, _, k) = _is_exterior_power(V); fl)
    return induced_map_on_exterior_power(h, k; domain=V, codomain=W)
  elseif ((fl, _, k) = _is_symmetric_power(V); fl)
    return induced_map_on_symmetric_power(h, k; domain=V, codomain=W)
  elseif ((fl, _, k) = _is_tensor_power(V); fl)
    return induced_map_on_tensor_power(h, k; domain=V, codomain=W)
  else
    throw(ArgumentError("First module must be a power module"))
  end
end

function _induced_map_on_power(
  D::LieAlgebraModule, C::LieAlgebraModule, h::LieAlgebraModuleHom, power::Int, type::Symbol
)
  TD = type == :tensor ? D : get_attribute(D, :embedding_tensor_power)::typeof(D)
  TC = type == :tensor ? C : get_attribute(C, :embedding_tensor_power)::typeof(C)

  mat = reduce(
    kronecker_product,
    [matrix(h) for _ in 1:power];
    init=identity_matrix(coefficient_ring(C), 1),
  )
  TD_to_TC = hom(TD, TC, mat; check=false)

  if type == :tensor
    return TD_to_TC
  else
    D_to_TD = get_attribute(
      D, :embedding_tensor_power_embedding
    )::LieAlgebraModuleHom{typeof(D),typeof(TD)}
    TC_to_C = get_attribute(
      C, :embedding_tensor_power_projection
    )::LieAlgebraModuleHom{typeof(TC),typeof(C)}
    return D_to_TD * TD_to_TC * TC_to_C
  end
end

function induced_map_on_exterior_power(
  h::LieAlgebraModuleHom,
  k::Int;
  domain::LieAlgebraModule{C}=exterior_power(Oscar.domain(h), k)[1],
  codomain::LieAlgebraModule{C}=exterior_power(Oscar.codomain(h), k)[1],
) where {C<:FieldElem}
  (domain_fl, domain_base, domain_k) = _is_exterior_power(domain)
  (codomain_fl, codomain_base, codomain_k) = _is_exterior_power(codomain)
  @req domain_fl "Domain must be an exterior power"
  @req codomain_fl "Codomain must be an exterior power"
  @req k == domain_k == codomain_k "Exponent mismatch"
  @req Oscar.domain(h) === domain_base && Oscar.codomain(h) === codomain_base "Domain/codomain mismatch"

  return _induced_map_on_power(domain, codomain, h, k, :ext)
end

function induced_map_on_symmetric_power(
  h::LieAlgebraModuleHom,
  k::Int;
  domain::LieAlgebraModule{C}=symmetric_power(Oscar.domain(h), k)[1],
  codomain::LieAlgebraModule{C}=symmetric_power(Oscar.codomain(h), k)[1],
) where {C<:FieldElem}
  (domain_fl, domain_base, domain_k) = _is_symmetric_power(domain)
  (codomain_fl, codomain_base, codomain_k) = _is_symmetric_power(codomain)
  @req domain_fl "Domain must be an symmetric power"
  @req codomain_fl "Codomain must be an symmetric power"
  @req k == domain_k == codomain_k "Exponent mismatch"
  @req Oscar.domain(h) === domain_base && Oscar.codomain(h) === codomain_base "Domain/codomain mismatch"

  return _induced_map_on_power(domain, codomain, h, k, :sym)
end

function induced_map_on_tensor_power(
  h::LieAlgebraModuleHom,
  k::Int;
  domain::LieAlgebraModule{C}=tensor_power(Oscar.domain(h), k)[1],
  codomain::LieAlgebraModule{C}=tensor_power(Oscar.codomain(h), k)[1],
) where {C<:FieldElem}
  (domain_fl, domain_base, domain_k) = _is_tensor_power(domain)
  (codomain_fl, codomain_base, codomain_k) = _is_tensor_power(codomain)
  @req domain_fl "Domain must be an tensor power"
  @req codomain_fl "Codomain must be an tensor power"
  @req domain_k == codomain_k "Exponent mismatch"
  @req Oscar.domain(h) === domain_base && Oscar.codomain(h) === codomain_base "Domain/codomain mismatch"

  return _induced_map_on_power(domain, codomain, h, k, :tensor)
end
