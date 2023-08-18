mutable struct LieAlgebraHom{T1<:LieAlgebra,T2<:LieAlgebra} <:
               Map{T1,T2,Hecke.HeckeMap,LieAlgebraHom}
  header::MapHeader
  matrix::MatElem

  inverse_isomorphism::LieAlgebraHom{T2,T1}

  function LieAlgebraHom(
    L1::LieAlgebra, L2::LieAlgebra, imgs::Vector{<:LieAlgebraElem}; check::Bool=true
  )
    @req all(x -> parent(x) === L2, imgs) "Images must lie in the codomain"
    @req length(imgs) == dim(L1) "Number of images must match dimension of domain"

    mat = zero_matrix(coefficient_ring(L2), dim(L1), dim(L2))
    for (i, img) in enumerate(imgs)
      mat[i, :] = _matrix(img)
    end
    return LieAlgebraHom(L1, L2, mat; check)
  end

  function LieAlgebraHom(L1::LieAlgebra, L2::LieAlgebra, mat::MatElem; check::Bool=true)
    h = new{typeof(L1),typeof(L2)}()
    h.matrix = mat::dense_matrix_type(coefficient_ring(L2))
    h.header = MapHeader(L1, L2)
    if check
      for x1 in basis(L1), x2 in basis(L1)
        @req h(x1) * h(x2) == h(x1 * x2) "Not a homomorphism"
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

function matrix(h::LieAlgebraHom)
  return h.matrix
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", h::LieAlgebraHom)
  io = pretty(io)
  println(io, LowercaseOff(), "Lie algebra morphism")
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(h))
  print(io, "to   ", Lowercase(), domain(h))
  print(io, Dedent())
end

function Base.show(io::IO, h::LieAlgebraHom)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, LowercaseOff(), "Lie algebra morphism")
  else
    print(io, LowercaseOff(), "Lie algebra morphism: ")
    print(io, Lowercase(), domain(h), " -> ", Lowercase(), domain(h))
  end
end

###############################################################################
#
#   Image and kernel
#
###############################################################################

function image(h::LieAlgebraHom, x::LieAlgebraElem)
  @req parent(x) === domain(h) "Domain mismatch"
  return codomain(h)(_matrix(x) * h.matrix)
end

function image(h::LieAlgebraHom)
  return sub(codomain(h), [image(h, x) for x in basis(domain(h))])
end

function image(h::LieAlgebraHom, I::LieAlgebraIdeal)
  @req base_lie_algebra(I) === domain(h) "Domain mismatch"
  return sub(codomain(h), [image(h, x) for x in basis(I)])
end

function image(h::LieAlgebraHom, S::LieSubalgebra)
  @req base_lie_algebra(S) === domain(h) "Domain mismatch"
  return sub(codomain(h), [image(h, x) for x in basis(S)])
end

function kernel(h::LieAlgebraHom)
  ker_dim, ker_b = left_kernel(matrix(h))
  return ideal(domain(h), [domain(h)(ker_b[i, :]) for i in 1:ker_dim])
end

###############################################################################
#
#   Map operations
#
###############################################################################

function compose(
  f::LieAlgebraHom{T1,T2}, g::LieAlgebraHom{T2,T3}
) where {T1<:LieAlgebra,T2<:LieAlgebra,T3<:LieAlgebra}
  @req codomain(f) === domain(g) "Composition: Maps are not compatible"
  h = LieAlgebraHom(domain(f), codomain(g), f.matrix * g.matrix; check=false)
  if isdefined(f, :inverse_isomorphism) && isdefined(g, :inverse_isomorphism)
    h.inverse_isomorphism = LieAlgebraHom(
      codomain(g),
      domain(f),
      g.inverse_isomorphism.matrix * f.inverse_isomorphism.matrix;
      check=false,
    )
  end
  return h
end

function inv(h::LieAlgebraHom{T1,T2}) where {T1<:LieAlgebra,T2<:LieAlgebra}
  @req is_isomorphism(h) "Homomorphism must be invertible"
  return h.inverse_isomorphism
end

function is_isomorphism(h::LieAlgebraHom{T1,T2}) where {T1<:LieAlgebra,T2<:LieAlgebra}
  isdefined(h, :inverse_isomorphism) && return true
  fl, invmat = is_invertible_with_inverse(h.matrix)
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

function hom(
  L1::LieAlgebra{C}, L2::LieAlgebra{C}, imgs::Vector{<:LieAlgebraElem{C}}; check::Bool=true
) where {C<:RingElement}
  return LieAlgebraHom(L1, L2, imgs; check)
end

function hom(
  L1::LieAlgebra{C}, L2::LieAlgebra{C}, mat::MatElem{C}; check::Bool=true
) where {C<:RingElement}
  return LieAlgebraHom(L1, L2, mat; check)
end

function identity_map(L::LieAlgebra)
  return hom(L, L, basis(L); check=false)
end
