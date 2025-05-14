abstract type ModuleData{C <: FieldElem, LieT <: LieAlgebraElem{C}} end

# To be implemented by subtypes:
# Mandatory:
#   base_lie_algebra(V::MyModuleData) -> LieAlgebra
#   highest_weight(V::MyModuleData) -> WeightLatticeElem
#   dim(V::MyModuleData) -> ZZRingElem
#   character(V::MyModuleData) -> Dict{WeightLatticeElem,ZZRingElem}


mutable struct SimpleModuleData{C <: FieldElem, LieT <: LieAlgebraElem{C}} <: ModuleData{C, LieT}
    L::LieAlgebra{C} # parent_type(LieT)
    highest_weight::WeightLatticeElem

    # The following fields are not set by default, just for caching
    dim::ZZRingElem
    character::Dict{WeightLatticeElem, ZZRingElem}

    function SimpleModuleData(L::LieAlgebra{C}, highest_weight::WeightLatticeElem) where {C <: FieldElem}
        return new{C, elem_type(L)}(L, highest_weight)
    end
end

function SimpleModuleData(L::LieAlgebra, highest_weight::Vector{Int})
    return SimpleModuleData(L, WeightLatticeElem(root_system(L), highest_weight))
end

function base_lie_algebra(V::SimpleModuleData{C, LieT}) where {C <: FieldElem, LieT <: LieAlgebraElem{C}}
    return V.L::parent_type(LieT)
end

function highest_weight(V::SimpleModuleData)
    return V.highest_weight
end

function dim(V::SimpleModuleData)
    if !isdefined(V, :dim)
        V.dim = dim_of_simple_module(ZZRingElem, base_lie_algebra(V), highest_weight(V))
    end
    return V.dim
end

function character(V::SimpleModuleData)
    if !isdefined(V, :character)
        V.character = character(ZZRingElem, base_lie_algebra(V), highest_weight(V))
    end
    return V.character
end
