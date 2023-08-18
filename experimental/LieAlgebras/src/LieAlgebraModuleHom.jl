mutable struct LieAlgebraModuleHom{T1<:LieAlgebraModule,T2<:LieAlgebraModule} <:
               Map{T1,T2,Hecke.HeckeMap,LieAlgebraModuleHom}
  header::MapHeader
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

function matrix(h::LieAlgebraModuleHom)
  return h.matrix
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
  print(io, "to   ", Lowercase(), domain(h))
  print(io, Dedent())
end

function Base.show(io::IO, h::LieAlgebraModuleHom)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, LowercaseOff(), "Lie algebra module morphism")
  else
    print(io, LowercaseOff(), "Lie algebra module morphism: ")
    print(io, Lowercase(), domain(h), " -> ", Lowercase(), domain(h))
  end
end

###############################################################################
#
#   Image and kernel
#
###############################################################################

function image(
  h::LieAlgebraModuleHom{T1,T2}, v::LieAlgebraModuleElem
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule}
  @req parent(v) === domain(h) "Domain mismatch"
  return codomain(h)(_matrix(v) * h.matrix)
end

# TODO: image and kernel, once submodules are implemented

###############################################################################
#
#   Map operations
#
###############################################################################

function compose(
  f::LieAlgebraModuleHom{T1,T2}, g::LieAlgebraModuleHom{T2,T3}
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule,T3<:LieAlgebraModule}
  @req codomain(f) === domain(g) "Composition: Maps are not compatible"
  h = LieAlgebraModuleHom(domain(f), codomain(g), f.matrix * g.matrix; check=false)
  if isdefined(f, :inverse_isomorphism) && isdefined(g, :inverse_isomorphism)
    h.inverse_isomorphism = LieAlgebraModuleHom(
      codomain(g),
      domain(f),
      g.inverse_isomorphism.matrix * f.inverse_isomorphism.matrix;
      check=false,
    )
  end
  return h
end

function inv(
  h::LieAlgebraModuleHom{T1,T2}
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule}
  @req is_isomorphism(h) "Homomorphism must be invertible"
  return h.inverse_isomorphism
end

function is_isomorphism(
  h::LieAlgebraModuleHom{T1,T2}
) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule}
  isdefined(h, :inverse_isomorphism) && return true
  fl, invmat = is_invertible_with_inverse(h.matrix)
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

function hom(
  V1::LieAlgebraModule{C},
  V2::LieAlgebraModule{C},
  imgs::Vector{<:LieAlgebraModuleElem{C}};
  check::Bool=true,
) where {C<:RingElement}
  return LieAlgebraModuleHom(V1, V2, imgs; check)
end

function hom(
  V1::LieAlgebraModule{C}, V2::LieAlgebraModule{C}, mat::MatElem{C}; check::Bool=true
) where {C<:RingElement}
  return LieAlgebraModuleHom(V1, V2, mat; check)
end

function identity_map(V::LieAlgebraModule)
  return hom(V, V, basis(V); check=false)
end
