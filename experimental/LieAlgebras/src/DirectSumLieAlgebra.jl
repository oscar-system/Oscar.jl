###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{DirectSumLieAlgebraElem{C}}) where {C<:FieldElem} =
  DirectSumLieAlgebra{C}

elem_type(::Type{DirectSumLieAlgebra{C}}) where {C<:FieldElem} = DirectSumLieAlgebraElem{C}

parent(x::DirectSumLieAlgebraElem) = x.parent

coefficient_ring(L::DirectSumLieAlgebra{C}) where {C<:FieldElem} = L.R::parent_type(C)

dim(L::DirectSumLieAlgebra) = L.dim

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, mime::MIME"text/plain", L::DirectSumLieAlgebra)
  @show_name(io, L)
  @show_special(io, mime, L)
  io = pretty(io)
  println(io, "Direct sum Lie algebra")
  println(io, Indent(), "of dimension $(dim(L))", Dedent())
  println(io, "with summands", Indent())
  for S in L.summands
    println(terse(io), Lowercase(), S)
  end
  print(io, Dedent(), "over ")
  print(io, Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::DirectSumLieAlgebra)
  @show_name(io, L)
  @show_special(io, L)
  if is_terse(io)
    print(io, "Direct sum Lie algebra")
  else
    io = pretty(io)
    print(
      io,
      "Direct sum of $(_number_of_direct_product_factors(L)) Lie algebras over ",
      Lowercase(),
    )
    print(terse(io), coefficient_ring(L))
  end
end

function symbols(L::DirectSumLieAlgebra)
  return L.s
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function bracket(
  x::DirectSumLieAlgebraElem{C}, y::DirectSumLieAlgebraElem{C}
) where {C<:FieldElem}
  check_parent(x, y)
  D = parent(x)
  vec = zero(_matrix(x))
  index = 0
  for S in D.summands
    offset = index
    dimS = dim(S)
    index += dimS

    xS_vec = _matrix(x)[:, (offset + 1):(offset + dimS)]
    iszero(xS_vec) && continue
    xS = S(xS_vec)
    yS_vec = _matrix(y)[:, (offset + 1):(offset + dimS)]
    iszero(yS_vec) && continue
    yS = S(_matrix(y)[:, (offset + 1):(offset + dimS)])
    vec[:, (offset + 1):(offset + dimS)] = _matrix(bracket(xS, yS))
  end
  return D(vec)
end

###############################################################################
#
#   Properties
#
###############################################################################

@attr Bool function is_abelian(L::DirectSumLieAlgebra)
  return all(is_abelian, L.summands)
end

###############################################################################
#
#   Direct sum functionality
#
###############################################################################

AbstractAlgebra._number_of_direct_product_factors(D::DirectSumLieAlgebra) =
  length(D.summands)

function canonical_injection(D::DirectSumLieAlgebra, i::Int)
  @boundscheck @req 1 <= i <= _number_of_direct_product_factors(D) "Invalid index."
  S = D.summands[i]
  mat = Generic.inj_proj_mat(
    coefficient_ring(D), dim(S), dim(D), sum(dim, D.summands[1:(i - 1)]; init=0) + 1
  )
  return hom(S, D, mat; check=false)
end

function canonical_projection(D::DirectSumLieAlgebra, i::Int)
  @boundscheck @req 1 <= i <= _number_of_direct_product_factors(D) "Invalid index."
  S = D.summands[i]
  mat = Generic.inj_proj_mat(
    coefficient_ring(D), dim(D), dim(S), sum(dim, D.summands[1:(i - 1)]; init=0) + 1
  )
  return hom(D, S, mat; check=false)
end

###############################################################################
#
#   Root system getters
#
###############################################################################

has_root_system(D::DirectSumLieAlgebra) = isdefined(D, :root_system)

function root_system(D::DirectSumLieAlgebra)
  assure_root_system(D)
  return D.root_system
end

function chevalley_basis(D::DirectSumLieAlgebra)
  assure_root_system(D)
  return D.chevalley_basis::NTuple{3,Vector{elem_type(D)}}
end

function set_root_system_and_chevalley_basis!(
  D::DirectSumLieAlgebra{C},
  R::RootSystem,
  chev::NTuple{3,Vector{DirectSumLieAlgebraElem{C}}},
) where {C<:FieldElem}
  D.root_system = R
  D.chevalley_basis = chev
end

###############################################################################
#
#   Constructor
#
###############################################################################

function direct_sum(Ss::Vector{<:LieAlgebra{C}}) where {C<:FieldElem}
  R = coefficient_ring(Ss[1])
  return direct_sum(R, Ss)
end

function direct_sum(R::Field, Ss::Vector{<:LieAlgebra{C}}) where {C<:FieldElem}
  @req elem_type(R) == C "Coefficient ring mismatch."
  return DirectSumLieAlgebra{C}(R, Ss)
end

function direct_sum(S::LieAlgebra{C}, Ss::LieAlgebra{C}...) where {C<:FieldElem}
  return direct_sum(LieAlgebra{C}[S; Ss...])
end
