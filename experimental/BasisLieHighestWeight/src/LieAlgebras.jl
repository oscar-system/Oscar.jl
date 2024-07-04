@attributes mutable struct LieAlgebraStructure
  lie_type::Symbol
  rank::Int
  lie_algebra_gap::GAP.Obj
  chevalley_basis::NTuple{3,Vector{GAP.Obj}}

  function LieAlgebraStructure(lie_type::Symbol, rank::Int)
    lie_algebra_gap = GAP.Globals.SimpleLieAlgebra(
      GAP.Obj(lie_type), rank, GAP.Globals.Rationals
    )
    chevalley_basis = NTuple{3,Vector{GAP.Obj}}(GAP.Globals.ChevalleyBasis(lie_algebra_gap))
    return new(lie_type, rank, lie_algebra_gap, chevalley_basis)
  end
end

rank(L::LieAlgebraStructure) = L.rank

@attr QQMatrix function cartan_matrix(L::LieAlgebraStructure)
  R = GAPWrap.RootSystem(L.lie_algebra_gap)
  C = matrix(QQ, GAP.Globals.CartanMatrix(R))
  return C
end

@attr QQMatrix function inv_cartan_matrix(L::LieAlgebraStructure)
  return inv(cartan_matrix(L))
end

function Base.show(io::IO, L::LieAlgebraStructure)
  io = pretty(io)
  print(io, LowercaseOff(), "Lie algebra of type $(L.lie_type)$(L.rank)")
end

function lie_algebra(type::Symbol, rk::Int)
  return LieAlgebraStructure(type, rk)
end

function chevalley_basis_gap(L::LieAlgebraStructure)
  return L.chevalley_basis
end

function cartan_sub_basis(L::LieAlgebraStructure)
  return L.chevalley_basis[3]
end

function root_system_gap(L::LieAlgebraStructure)
  return GAPWrap.RootSystem(L.lie_algebra_gap)
end

function number_of_positive_roots(L::LieAlgebraStructure)
  return length(GAP.Globals.PositiveRoots(root_system_gap(L)))
end

function dim_of_simple_module(T::Type, L::LieAlgebraStructure, hw::Vector{<:IntegerUnion})
  hw_ = Int.(hw)
  @req Oscar.LieAlgebras.is_dominant_weight(hw_) "Not a dominant weight."
  return T(GAPWrap.DimensionOfHighestWeightModule(L.lie_algebra_gap, GAP.Obj(Int.(hw_))))
end

function dim_of_simple_module(L::LieAlgebraStructure, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(Int, L, hw)
end

function matrices_of_operators_gap(
  L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}, operators::Vector{GAP.Obj}
)
  # used in tensor_matrices_of_operators
  M = GAP.Globals.HighestWeightModule(L.lie_algebra_gap, GAP.Obj(Int.(highest_weight)))
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
  basis = GAPWrap.Basis(L.lie_algebra_gap)
  basis_ind = GAP.Globals.Position(basis, operator)
  denom = GAPWrap.Coefficients(basis, operator)[basis_ind]
  return [
    ZZ(GAPWrap.Coefficients(basis, h * operator)[basis_ind]//denom) for
    h in cartan_sub_basis(L)
  ]
end
