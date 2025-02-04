###############################################################################
#
#   Representation theory for semisimple Lie algebras
#
###############################################################################

# TODO: add semisimplicity check once that is available

###############################################################################
#
#   Simple modules (via highest weight) of semisimple Lie algebras
#
###############################################################################

# TODO: move to RootSystem.jl
function is_dominant_weight(hw::Vector{<:IntegerUnion})
  return all(>=(0), hw)
end

@doc raw"""
    simple_module(L::LieAlgebra{C}, hw::Vector{Int}) -> LieAlgebraModule{C}

Construct the simple module of the Lie algebra `L` with highest weight `hw`.
"""
function simple_module(L::LieAlgebra, hw::Vector{Int})
  @req is_dominant_weight(hw) "Not a dominant weight."
  struct_consts = lie_algebra_simple_module_struct_consts_gap(L, hw)
  dimV = size(struct_consts, 2)
  V = abstract_module(L, dimV, struct_consts; check=false)
  # TODO: set appropriate attributes
  return V
end

@doc raw"""
    dim_of_simple_module([T = Int], L::LieAlgebra{C}, hw::WeightLatticeElem) -> T
    dim_of_simple_module([T = Int], L::LieAlgebra{C}, hw::Vector{<:IntegerUnion}) -> T

Compute the dimension of the simple module of the Lie algebra `L` with highest weight `hw`
using Weyl's dimension formula.
The return value is of type `T`.

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :A, 3);

julia> dim_of_simple_module(L, [1, 1, 1])
64
```
"""
function dim_of_simple_module(L::LieAlgebra, hw::WeightLatticeElem)
  return dim_of_simple_module(root_system(L), hw)
end

function dim_of_simple_module(T::Type, L::LieAlgebra, hw::WeightLatticeElem)
  return dim_of_simple_module(T, root_system(L), hw)
end

function dim_of_simple_module(L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(root_system(L), hw)
end

function dim_of_simple_module(T::Type, L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return dim_of_simple_module(T, root_system(L), hw)
end

@doc raw"""
    dominant_weights(L::LieAlgebra{C}, hw::WeightLatticeElem) -> Vector{WeightLatticeElem}
    dominant_weights(L::LieAlgebra{C}, hw::Vector{<:IntegerUnion}) -> Vector{WeightLatticeElem}

Computes the dominant weights occurring in the simple module of the Lie algebra `L` with highest weight `hw`,
sorted ascendingly by the total height of roots needed to reach them from `hw`.

See [MP82](@cite) for details and the implemented algorithm.

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :B, 3);

julia> dominant_weights(L, [1, 0, 3])
7-element Vector{WeightLatticeElem}:
 w_1 + 3*w_3
 w_1 + w_2 + w_3
 2*w_1 + w_3
 3*w_3
 w_2 + w_3
 w_1 + w_3
 w_3
```
"""
function dominant_weights(L::LieAlgebra, hw::WeightLatticeElem)
  return dominant_weights(root_system(L), hw)
end

function dominant_weights(L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return dominant_weights(root_system(L), hw)
end

@doc raw"""
    dominant_character([T = Int], L::LieAlgebra{C}, hw::WeightLatticeElem) -> Dict{WeightLatticeElem, T}
    dominant_character([T = Int], L::LieAlgebra{C}, hw::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}

Computes the dominant weights occurring in the simple module of the Lie algebra `L` with highest weight `hw`,
together with their multiplicities.

This function uses an optimized version of the Freudenthal formula, see [MP82](@cite) for details.

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :A, 3);

julia> dominant_character(L, [2, 1, 0])
Dict{WeightLatticeElem, Int64} with 4 entries:
  0           => 3
  2*w_2       => 1
  w_1 + w_3   => 2
  2*w_1 + w_2 => 1
```
"""
function dominant_character(L::LieAlgebra, hw::WeightLatticeElem)
  return dominant_character(root_system(L), hw)
end

function dominant_character(T::DataType, L::LieAlgebra, hw::WeightLatticeElem)
  return dominant_character(T, root_system(L), hw)
end

function dominant_character(L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return dominant_character(root_system(L), hw)
end

function dominant_character(T::DataType, L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return dominant_character(T, root_system(L), hw)
end

@doc raw"""
    character([T = Int], L::LieAlgebra{C}, hw::WeightLatticeElem) -> Dict{WeightLatticeElem, T}
    character([T = Int], L::LieAlgebra{C}, hw::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}

Computes all weights occurring in the simple module of the Lie algebra `L` with highest weight `hw`,
together with their multiplicities.
This is achieved by acting with the Weyl group on the [`dominant_character`](@ref dominant_character(::LieAlgebra, ::Vector{<:IntegerUnion})).

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :A, 3);

julia> character(L, [2, 0, 0])
Dict{WeightLatticeElem, Int64} with 10 entries:
  -2*w_3           => 1
  -2*w_2 + 2*w_3   => 1
  2*w_1            => 1
  -2*w_1 + 2*w_2   => 1
  -w_1 + w_3       => 1
  w_2              => 1
  w_1 - w_2 + w_3  => 1
  w_1 - w_3        => 1
  -w_1 + w_2 - w_3 => 1
  -w_2             => 1
```
"""
function character(L::LieAlgebra, hw::WeightLatticeElem)
  return character(root_system(L), hw)
end

function character(T::DataType, L::LieAlgebra, hw::WeightLatticeElem)
  return character(T, root_system(L), hw)
end

function character(L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return character(root_system(L), hw)
end

function character(T::DataType, L::LieAlgebra, hw::Vector{<:IntegerUnion})
  return character(T, root_system(L), hw)
end

@doc raw"""
    tensor_product_decomposition(L::LieAlgebra, hw1::WeightLatticeElem, hw2::WeightLatticeElem) -> MSet{Vector{Int}}
    tensor_product_decomposition(L::LieAlgebra, hw1::Vector{<:IntegerUnion}, hw2::Vector{<:IntegerUnion}) -> MSet{Vector{Int}}

Computes the decomposition of the tensor product of the simple modules of the Lie algebra `L` with highest weights `hw1` and `hw2`
into simple modules with their multiplicities.
This function uses Klimyk's formula (see [Hum72; Exercise 24.9](@cite)).

The return type may change in the future.

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :A, 2);

julia> tensor_product_decomposition(L, [1, 0], [0, 1])
MSet{Vector{Int64}} with 2 elements:
  [0, 0]
  [1, 1]

julia> tensor_product_decomposition(L, [1, 1], [1, 1])
MSet{Vector{Int64}} with 6 elements:
  [0, 0]
  [1, 1] : 2
  [2, 2]
  [3, 0]
  [0, 3]
```
"""
function tensor_product_decomposition(
  L::LieAlgebra, hw1::Vector{<:IntegerUnion}, hw2::Vector{<:IntegerUnion}
)
  return tensor_product_decomposition(root_system(L), hw1, hw2)
end

function tensor_product_decomposition(
  L::LieAlgebra, hw1::WeightLatticeElem, hw2::WeightLatticeElem
)
  return tensor_product_decomposition(root_system(L), hw1, hw2)
end


###############################################################################
#
#   Demazure modules (via highest weight) of semisimple Lie algebras
#
###############################################################################

@doc raw"""
    demazure_character([T = Int], L::LieAlgebra, w::WeightLatticeElem, x::WeylGroupElem) -> Dict{WeightLatticeElem, T}
    demazure_character([T = Int], L::LieAlgebra, w::Vector{<:IntegerUnion}, x::WeylGroupElem) -> Dict{WeightLatticeElem, T}
    demazure_character([T = Int], L::LieAlgebra, w::WeightLatticeElem, reduced_expr::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}
    demazure_character([T = Int], L::LieAlgebra, w::Vector{<:IntegerUnion}, reduced_expr::Vector{<:IntegerUnion}) -> Dict{WeightLatticeElem, T}

Computes all weights occurring in the Demazure module of the Lie algebra `L``
with extremal weight `w * x`, together with their multiplicities.

Instead of a Weyl group element `x`, a reduced expression for `x` can be supplied.
This function may return arbitrary results if the provided expression is not reduced.

# Examples
```jldoctest
julia> L = lie_algebra(QQ, :A, 2);

julia> demazure_character(L, [1, 1], [2, 1])
Dict{WeightLatticeElem, Int64} with 5 entries:
  2*w_1 - w_2  => 1
  w_1 + w_2    => 1
  0            => 1
  -w_1 + 2*w_2 => 1
  -2*w_1 + w_2 => 1
```
"""
function demazure_character(L::LieAlgebra, w::WeightLatticeElem, x::WeylGroupElem)
  return demazure_character(root_system(L), w, x)
end

function demazure_character(
  T::DataType, L::LieAlgebra, w::WeightLatticeElem, x::WeylGroupElem
)
  return demazure_character(T, root_system(L), w, x)
end

function demazure_character(
  L::LieAlgebra, w::WeightLatticeElem, reduced_expression::Vector{<:IntegerUnion}
)
  return demazure_character(root_system(L), w, reduced_expression)
end

function demazure_character(
  T::DataType,
  L::LieAlgebra,
  w::WeightLatticeElem,
  reduced_expression::Vector{<:IntegerUnion},
)
  return demazure_character(T, root_system(L), w, reduced_expression)
end

function demazure_character(L::LieAlgebra, w::Vector{<:IntegerUnion}, x::WeylGroupElem)
  return demazure_character(root_system(L), w, x)
end

function demazure_character(
  T::DataType, L::LieAlgebra, w::Vector{<:IntegerUnion}, x::WeylGroupElem
)
  return demazure_character(T, root_system(L), w, x)
end

function demazure_character(
  L::LieAlgebra, w::Vector{<:IntegerUnion}, reduced_expression::Vector{<:IntegerUnion}
)
  return demazure_character(root_system(L), w, reduced_expression)
end

function demazure_character(
  T::DataType,
  L::LieAlgebra,
  w::Vector{<:IntegerUnion},
  reduced_expression::Vector{<:IntegerUnion},
)
  return demazure_character(T, root_system(L), w, reduced_expression)
end
