@attributes mutable struct LieAlgebraStructure
  lie_type::Symbol
  rank::Int
  lie_algebra::AbstractLieAlgebra{QQFieldElem}
  chevalley_basis_gap::NTuple{3,Vector{GAP.Obj}}

  function LieAlgebraStructure(lie_type::Symbol, rank::Int)
    lie_algebra = Oscar.lie_algebra(QQ, lie_type, rank)
    chevalley_basis_gap = NTuple{3,Vector{GAP.Obj}}(GAP.Globals.ChevalleyBasis(codomain(Oscar.iso_oscar_gap(lie_algebra))))
    return new(lie_type, rank, lie_algebra, chevalley_basis_gap)
  end
end

rank(L::LieAlgebraStructure) = L.rank

function cartan_matrix(L::LieAlgebraStructure)
  return cartan_matrix(L.lie_algebra)
end

function cartan_matrix_inv(L::LieAlgebraStructure)
  return cartan_matrix_inv(L.lie_algebra)
end

function Base.show(io::IO, L::LieAlgebraStructure)
  io = pretty(io)
  print(io, LowercaseOff(), "Lie algebra of type $(L.lie_type)$(L.rank)")
end

function lie_algebra(type::Symbol, rk::Int)
  return LieAlgebraStructure(type, rk)
end

function chevalley_basis_gap(L::LieAlgebraStructure)
  return L.chevalley_basis_gap
end

function cartan_sub_basis_gap(L::LieAlgebraStructure)
  return L.chevalley_basis_gap[3]
end

function root_system(L::LieAlgebraStructure)
  return root_system(L.lie_algebra)
end

function dim_of_simple_module(T::Type, L::LieAlgebraStructure, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(T, L.lie_algebra, hw)
end

function dim_of_simple_module(L::LieAlgebraStructure, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(Int, L, hw)
end

function matrices_of_operators_gap(
  L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}, operators::Vector{GAP.Obj}
)
  # used in tensor_matrices_of_operators
  M = GAP.Globals.HighestWeightModule(codomain(Oscar.iso_oscar_gap(L.lie_algebra)), GAP.Obj(Int.(highest_weight)))
  matrices_of_operators = [
    sparse_matrix(matrix(QQ, GAP.Globals.MatrixOfAction(GAPWrap.Basis(M), o))) for
    o in operators
  ]
  denominators = map(y -> denominator(y[2]), union(union(matrices_of_operators...)...))
  common_denominator = lcm(denominators)# // 1
  matrices_of_operators =
    (A -> change_base_ring(ZZ, common_denominator * A)).(matrices_of_operators)
  return matrices_of_operators
end

@doc raw"""
    weight(L::LieAlgebraStructure, operator::GAP.Obj) -> Vector{ZZRingElem}

Calculate the weight of `operator` w.r.t. the fundamental weights w_i.
"""
function weight(L::LieAlgebraStructure, operator::GAP.Obj)
  @req !iszero(operator) "Operators should be non-zero"
  basis = GAPWrap.Basis(codomain(Oscar.iso_oscar_gap(L.lie_algebra)))
  basis_ind = GAP.Globals.Position(basis, operator)
  denom = GAPWrap.Coefficients(basis, operator)[basis_ind]
  return [
    ZZ(GAPWrap.Coefficients(basis, h * operator)[basis_ind]//denom) for
    h in cartan_sub_basis_gap(L)
  ]
end
