@attributes mutable struct LieAlgebraStructure
  lie_type::Symbol
  rank::Int
  lie_algebra_gap::GAP.Obj

  function LieAlgebraStructure(lie_type::Symbol, rank::Int)
    lie_algebra_gap = GAP.Globals.SimpleLieAlgebra(
      GAP.Obj(lie_type), rank, GAP.Globals.Rationals
    )
    return new(lie_type, rank, lie_algebra_gap)
  end
end

rank(L::LieAlgebraStructure) = L.rank

@attr QQMatrix function cartan_matrix(L::LieAlgebraStructure)
  R = GAP.Globals.RootSystem(L.lie_algebra_gap)
  C = matrix(QQ, GAP.Globals.CartanMatrix(R))
  return C
end

@attr QQMatrix function inv_cartan_matrix(L::LieAlgebraStructure)
  return inv(cartan_matrix(L))
end

function Base.show(io::IO, lie_algebra::LieAlgebraStructure)
  print(io, "Lie-Algebra of type ", lie_algebra.lie_type, " and rank ", lie_algebra.rank)
end

function lie_algebra_with_basis(type::Symbol, rk::Int)
  lie_algebra = LieAlgebraStructure(type, rk)
  chevalley_basis = NTuple{3,Vector{GAP.Obj}}(
    GAP.Globals.ChevalleyBasis(lie_algebra.lie_algebra_gap)
  )
  return lie_algebra, chevalley_basis
end

function matricesForOperators(
  lie_algebra::GAP.Obj, highest_weight::Vector{ZZRingElem}, operators::Vector{GAP.Obj}
)::Vector{SMat{ZZRingElem}}
  """
  used to create tensorMatricesForOperators
  """
  M = GAP.Globals.HighestWeightModule(lie_algebra, GAP.Obj(Vector{Int}(highest_weight)))
  matrices_of_operators = [
    sparse_matrix(transpose(matrix(QQ, GAP.Globals.MatrixOfAction(GAPWrap.Basis(M), o))))
    for o in operators
  ]
  denominators = map(y -> denominator(y[2]), union(union(matrices_of_operators...)...))
  common_denominator = lcm(denominators)# // 1
  matrices_of_operators =
    (A -> change_base_ring(ZZ, common_denominator * A)).(matrices_of_operators)
  return matrices_of_operators
end

function weights_for_operators(
  lie_algebra::GAP.Obj, cartan_sub::Vector{GAP.Obj}, operators::Vector{GAP.Obj}
)::Vector{Vector{ZZRingElem}}
  """
  Calculates the weight weights[i] in w_i for each operator operators[i]
  """
  """cartan = [Vector{Int}(x) for x in GAP.Globals.ExtRepOfObj.(cartan)]
  operators = [Vector{Int}(x) for x in GAP.Globals.ExtRepOfObj.(operators)]
  if any(iszero.(operators))
      error("ops should be non-zero")
  end
  println([findfirst(v .!= 0) for v in operators])

  return [
      [(dot(h, v))[findfirst(v .!= 0)] / (v)[findfirst(v .!= 0)] for h in cartan] for v in operators
  ]
  """
  if any(iszero, operators)
    error("ops should be non-zero")
  end
  basis = GAP.Globals.Basis(lie_algebra)
  return [
    begin
      ind = GAP.Globals.Position(basis, v)
      denom = GAP.Globals.Coefficients(basis, v)[ind]
      [ZZ(GAP.Globals.Coefficients(basis, h * v)[ind]//denom) for h in cartan_sub]
    end for v in operators
  ]
end
