@attributes mutable struct AbstractLieAlgebra{C<:RingElement} <: LieAlgebra{C}
  R::Ring
  dim::Int
  struct_consts::Matrix{SRow{C}}
  s::Vector{Symbol}

  function AbstractLieAlgebra{C}(
    R::Ring,
    struct_consts::Matrix{SRow{C}},
    s::Vector{Symbol};
    cached::Bool=true,
    check::Bool=true,
  ) where {C<:RingElement}
    return get_cached!(
      AbstractLieAlgebraDict, (R, struct_consts, s), cached
    ) do
      (n1, n2) = size(struct_consts)
      @req n1 == n2 "Invalid structure constants dimensions."
      dimL = n1
      @req length(s) == dimL "Invalid number of basis element names."
      if check
        @req all(
          r -> all(e -> parent(last(e)) === R, r), struct_consts
        ) "Invalid structure constants."
        @req all(
          iszero, struct_consts[i, i][k] for i in 1:dimL, k in 1:dimL
        ) "Not anti-symmetric."
        @req all(
          iszero,
          struct_consts[i, j][k] + struct_consts[j, i][k] for i in 1:dimL, j in 1:dimL,
          k in 1:dimL
        ) "Not anti-symmetric."
        @req all(
          iszero,
          sum(
            struct_consts[i, j][k] * struct_consts[k, l][m] +
            struct_consts[j, l][k] * struct_consts[k, i][m] +
            struct_consts[l, i][k] * struct_consts[k, j][m] for k in 1:dimL
          ) for i in 1:dimL, j in 1:dimL, l in 1:dimL, m in 1:dimL
        ) "Jacobi identity does not hold."
      end
      new{C}(R, dimL, struct_consts, s)
    end::AbstractLieAlgebra{C}
  end
end

const AbstractLieAlgebraDict = CacheDictType{
  Tuple{Ring,Matrix{SRow},Vector{Symbol}},AbstractLieAlgebra
}()

struct AbstractLieAlgebraElem{C<:RingElement} <: LieAlgebraElem{C}
  parent::AbstractLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{AbstractLieAlgebraElem{C}}) where {C<:RingElement} =
  AbstractLieAlgebra{C}

elem_type(::Type{AbstractLieAlgebra{C}}) where {C<:RingElement} = AbstractLieAlgebraElem{C}

parent(x::AbstractLieAlgebraElem) = x.parent

coefficient_ring(L::AbstractLieAlgebra{C}) where {C<:RingElement} = L.R::parent_type(C)

dim(L::AbstractLieAlgebra) = L.dim

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", L::AbstractLieAlgebra)
  io = pretty(io)
  println(io, "Abstract Lie algebra")
  println(io, Indent(), "of dimension $(dim(L))", Dedent())
  print(io, "over ")
  print(io, Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::AbstractLieAlgebra)
  if get(io, :supercompact, false)
    print(io, "Abstract Lie algebra")
  else
    io = pretty(io)
    print(io, "Abstract Lie algebra over ", Lowercase())
    print(IOContext(io, :supercompact => true), coefficient_ring(L))
  end
end

function symbols(L::AbstractLieAlgebra)
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
  x::AbstractLieAlgebraElem{C}, y::AbstractLieAlgebraElem{C}
) where {C<:RingElement}
  check_parent(x, y)
  L = parent(x)
  mat = sum(
    cxi * cyj * L.struct_consts[i, j] for (i, cxi) in enumerate(coefficients(x)),
    (j, cyj) in enumerate(coefficients(y))
  )
  return L(mat)
end

###############################################################################
#
#   Properties
#
###############################################################################

function is_abelian(L::AbstractLieAlgebra)
  return all(e -> iszero(length(e)), L.struct_consts)
end

###############################################################################
#
#   Constructor
#
###############################################################################

@doc raw"""
    lie_algebra(R::Ring, struct_consts::Matrix{SRow{elem_type(R)}}, s::Vector{<:VarName}; cached::Bool, check::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the Lie algebra over the ring `R` with structure constants `struct_consts`
and with basis element names `s`.

The Lie bracket on the newly constructed Lie algebra `L` is determined by the structure
constants in `struct_consts` as follows: let $x_i$ denote the $i$-th standard basis vector
of `L`. Then the entry `struct_consts[i,j][k]` is a scalar $a_{i,j,k}$
such that $[x_i, x_j] = \sum_k a_{i,j,k} x_k$.

* `s`: A vector of basis element names. This is 
  `[Symbol("x_$i") for i in 1:size(struct_consts, 1)]` by default.
* `cached`: If `true`, cache the result. This is `true` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.
"""
function lie_algebra(
  R::Ring,
  struct_consts::Matrix{SRow{C}},
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:size(struct_consts, 1)];
  cached::Bool=true,
  check::Bool=true,
) where {C<:RingElement}
  @req C == elem_type(R) "Invalid coefficient type."
  return AbstractLieAlgebra{elem_type(R)}(R, struct_consts, Symbol.(s); cached, check)
end

@doc raw"""
    lie_algebra(R::Ring, struct_consts::Array{elem_type(R),3}, s::Vector{<:VarName}; cached::Bool, check::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the Lie algebra over the ring `R` with structure constants `struct_consts`
and with basis element names `s`.

The Lie bracket on the newly constructed Lie algebra `L` is determined by the structure
constants in `struct_consts` as follows: let $x_i$ denote the $i$-th standard basis vector
of `L`. Then the entry `struct_consts[i,j,k]` is a scalar $a_{i,j,k}$
such that $[x_i, x_j] = \sum_k a_{i,j,k} x_k$.

* `s`: A vector of basis element names. This is
  `[Symbol("x_$i") for i in 1:size(struct_consts, 1)]` by default.
* `cached`: If `true`, cache the result. This is `true` by default.
* `check`: If `true`, check that the structure constants are anti-symmetric and
  satisfy the Jacobi identity. This is `true` by default.

# Examples
```jldoctest
julia> struct_consts = zeros(QQ, 3, 3, 3);

julia> struct_consts[1, 2, 3] = QQ(1);

julia> struct_consts[2, 1, 3] = QQ(-1);

julia> struct_consts[3, 1, 1] = QQ(2);

julia> struct_consts[1, 3, 1] = QQ(-2);

julia> struct_consts[3, 2, 2] = QQ(-2);

julia> struct_consts[2, 3, 2] = QQ(2);

julia> sl2 = lie_algebra(QQ, struct_consts, ["e", "f", "h"])
Abstract Lie algebra
  of dimension 3
over rational field

julia> e, f, h = basis(sl2)
3-element Vector{AbstractLieAlgebraElem{QQFieldElem}}:
 e
 f
 h

julia> e * f
h

julia> h * e
2*e

julia> h * f
-2*f
```
"""
function lie_algebra(
  R::Ring,
  struct_consts::Array{C,3},
  s::Vector{<:VarName}=[Symbol("x_$i") for i in 1:size(struct_consts, 1)];
  cached::Bool=true,
  check::Bool=true,
) where {C<:RingElement}
  @req C == elem_type(R) "Invalid coefficient type."
  struct_consts2 = Matrix{SRow{elem_type(R)}}(
    undef, size(struct_consts, 1), size(struct_consts, 2)
  )
  for i in axes(struct_consts, 1), j in axes(struct_consts, 2)
    struct_consts2[i, j] = sparse_row(
      R, collect(axes(struct_consts, 3)), struct_consts[i, j, :]
    )
  end

  return AbstractLieAlgebra{elem_type(R)}(R, struct_consts2, Symbol.(s); cached, check)
end

@doc raw"""
    lie_algebra(R::Ring, dynkin::Tuple{Char,Int}; cached::Bool) -> AbstractLieAlgebra{elem_type(R)}

Construct the simple Lie algebra over the ring `R` with Dynkin type given by `dynkin`.
The actual construction is done in GAP.

If `cached` is `true`, the constructed Lie algebra is cached.
"""
function lie_algebra(R::Ring, dynkin::Tuple{Char,Int}; cached::Bool=true)
  @req dynkin[1] in 'A':'G' "Unknown Dynkin type"

  coeffs_iso = inv(Oscar.iso_oscar_gap(R))
  LG = GAP.Globals.SimpleLieAlgebra(
    GAP.Obj(string(dynkin[1])), dynkin[2], domain(coeffs_iso)
  )
  s = [Symbol("x_$i") for i in 1:GAPWrap.Dimension(LG)]
  LO = codomain(
    _iso_gap_oscar_abstract_lie_algebra(LG, s; coeffs_iso, cached)
  )::AbstractLieAlgebra{elem_type(R)}

  return LO
end
