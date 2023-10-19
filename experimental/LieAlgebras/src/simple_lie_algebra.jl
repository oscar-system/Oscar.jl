#Construct a simple Lie algebra over a given field with a given root system
@attributes mutable struct SimpleLieAlgebra{C<:FieldElem} <: LieAlgebra{C}
  R::Field
  root_system::RootSystem
  dim::Int
  s::Vector{Symbol}
  root_system_type::Tuple{Symbol,Int}
  struct_consts::Matrix{SRow{C}}
  function SimpleLieAlgebra{C}(
    R::Field, S::Symbol, n::Int; cached::Bool=true
  ) where {C<:FieldElem}
    RS = root_system(S, n)
    return get_cached!(SimpleLieAlgebraDict, (R, RS), cached) do
      dimL = number_of_roots(RS) + length(RS.simple_roots)
      s = [Symbol("e_$i") for i in 1:dimL]
      st = root_system_type(RS)
      #get the structure constants of the Lie algebra L 
      #note that it is enough to do this over QQ, as we can later coerce the constants
      #into the field R
      coeffs_iso = inv(Oscar.iso_oscar_gap(QQ))
      LG = GAP.Globals.SimpleLieAlgebra(GAP.Obj(S), n, domain(coeffs_iso))
      sc_table_G =
        (
          entry -> (entry[1], Vector{elem_type(QQ)}(map(coeffs_iso, entry[2])))
        ).(
          Matrix{Tuple{Vector{Int},Vector{GAP.Obj}}}(
            (GAP.Globals.StructureConstantsTable(GAPWrap.Basis(LG)))[1:dimL]
          )
        )
      struct_consts = Matrix{SRow{elem_type(R)}}(undef, dimL, dimL)
      for i in 1:dimL, j in 1:dimL
        struct_consts[i, j] = sparse_row(
          R, Tuple{Int,elem_type(R)}[(k, R(c)) for (k, c) in zip(sc_table_G[i, j]...)]
        )
      end
      new{C}(R, RS, dimL, s, st, struct_consts)
    end::SimpleLieAlgebra{C}
  end
end

const SimpleLieAlgebraDict = CacheDictType{Tuple{Field,RootSystem},SimpleLieAlgebra}()

mutable struct SimpleLieAlgebraElem{C<:FieldElem} <: LieAlgebraElem{C}
  parent::SimpleLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{SimpleLieAlgebraElem{C}}) where {C<:FieldElem} = SimpleLieAlgebra{C}

elem_type(::Type{SimpleLieAlgebra{C}}) where {C<:FieldElem} = SimpleLieAlgebraElem{C}

parent(x::SimpleLieAlgebraElem) = x.parent

coefficient_ring(L::SimpleLieAlgebra{C}) where {C<:FieldElem} = L.R::parent_type(C)

dim(L::SimpleLieAlgebra) = L.dim

root_system(L::SimpleLieAlgebra) = L.root_system

root_system_type(L::SimpleLieAlgebra) = L.root_system_type

function symbols(L::SimpleLieAlgebra)
  return L.s
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", L::SimpleLieAlgebra)
  io = pretty(io)
  println(io, "Simple Lie algebra")
  println(io, Indent(), "of type $(root_system_type_string(root_system(L)))", Dedent())
  println(io, Indent(), "of dimension $(dim(L))", Dedent())
  print(io, "over ")
  print(io, Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::SimpleLieAlgebra)
  if get(io, :supercompact, false)
    print(io, "Simple Lie algebra $(root_system_type_string(root_system(L)))")
  else
    io = pretty(io)
    print(
      io, "Simple Lie algebra $(root_system_type_string(root_system(L))) over ", Lowercase()
    )
    print(IOContext(io, :supercompact => true), coefficient_ring(L))
  end
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

# Binary operations

function bracket(
  x::SimpleLieAlgebraElem{C}, y::SimpleLieAlgebraElem{C}
) where {C<:FieldElem}
  check_parent(x, y)
  L = parent(x)
  mat = sum(
    cxi * cyj * L.struct_consts[i, j] for (i, cxi) in enumerate(coefficients(x)),
    (j, cyj) in enumerate(coefficients(y));
    init=sparse_row(coefficient_ring(L)),
  )
  return L(mat)
end
###############################################################################
#
#   Constructors
#
###############################################################################

@doc raw"""
    lie_algebra(R::Field, S::Symbol, n::Int; cached::Bool=true) -> SimpleLieAlgebra{elem_type(R)}

Construct the simple Lie algebra over the field `R` with root system of type `Sn` (see `root_system(S::Symbol, n::Int)` ref).
The internally used basis of this Lie algebra is the Chevalley basis.
"""
function lie_algebra(R::Field, S::Symbol, n::Int; cached::Bool=true)
  return SimpleLieAlgebra{elem_type(R)}(R, S, n; cached)
end

###############################################################################
#
#   Chevalley basis
#
###############################################################################

@doc raw"""
    chevalley_basis(L::SimpleLieAlgebra{T}) -> NTuple{3,Vector{SimpleLieAlgebraElem{T}}}

Return the Chevalley basis of the simple Lie algebra `L` in three vectors, stating first the positive root vectors, 
then the negative root vectors, and finally the basis of the Cartan subalgebra. The order of root vectors corresponds
to the order of the roots in the root system.
"""
function chevalley_basis(L::SimpleLieAlgebra)
  RS = root_system(L)
  n = length(RS.positive_roots)
  B = basis(L)
  # root vectors
  r_plus = B[1:n]
  r_minus = B[(n + 1):(2 * n)]
  # basis for cartan algebra
  h = B[(2 * n + 1):dim(L)]
  return (r_plus, r_minus, h)
end
