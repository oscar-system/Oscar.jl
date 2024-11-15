@attributes mutable struct LieAlgebraStructure
  lie_type::Symbol
  lie_algebra::AbstractLieAlgebra{QQFieldElem}

  function LieAlgebraStructure(lie_type::Symbol, rank::Int)
    lie_algebra = Oscar.lie_algebra(QQ, lie_type, rank)
    return new(lie_type, lie_algebra)
  end
end

rank(L::LieAlgebraStructure) = rank(root_system(L))

function cartan_matrix(L::LieAlgebraStructure)
  return cartan_matrix(L.lie_algebra)
end

function cartan_matrix_inv(L::LieAlgebraStructure)
  return cartan_matrix_inv(L.lie_algebra)
end

function dim(L::LieAlgebraStructure)
  return dim(L.lie_algebra)
end

function coefficient_ring(L::LieAlgebraStructure)
  return coefficient_ring(L.lie_algebra)
end

function Base.show(io::IO, L::LieAlgebraStructure)
  io = pretty(io)
  print(io, LowercaseOff(), "Lie algebra of type $(L.lie_type)$(rank(L))")
end

function lie_algebra(type::Symbol, rk::Int)
  return LieAlgebraStructure(type, rk)
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

function dim_of_simple_module(T::Type, L::LieAlgebraStructure, hw::WeightLatticeElem)
  return dim_of_simple_module(T, L.lie_algebra, hw)
end

function dim_of_simple_module(L::LieAlgebraStructure, hw::WeightLatticeElem)
  return dim_of_simple_module(Int, L, hw)
end

function matrices_of_operators(
  L::LieAlgebraStructure, highest_weight::WeightLatticeElem, operators::Vector{RootSpaceElem}
)
  # used in tensor_matrices_of_operators
  R = root_system(L)
  struct_consts = lie_algebra_simple_module_struct_consts_gap(L.lie_algebra, highest_weight)
  dimV = size(struct_consts, 2)
  transformation_matrices = [sparse_matrix(coefficient_ring(L), dimV, dimV) for _ in 1:number_of_positive_roots(R)]
  for i in 1:number_of_positive_roots(R), j in 1:dimV
    transformation_matrices[i][j] = struct_consts[i + number_of_positive_roots(R), j] # take f_alpha for positive root alpha
  end

  matrices_of_operators = map(operators) do op
    fl, i = is_positive_root_with_index(op)
    @assert fl
    change_base_ring(ZZ, transformation_matrices[i])
  end

  return matrices_of_operators
end
