abstract type ModuleData end

# To be implemented by subtypes:
# Mandatory:
#   base_lie_algebra(V::MyModuleData) -> LieAlgebra
#   highest_weight(V::MyModuleData) -> WeightLatticeElem
#   dim(V::MyModuleData) -> ZZRingElem
#   character(V::MyModuleData) -> Dict{WeightLatticeElem,ZZRingElem}


mutable struct SimpleModuleData <: ModuleData
    L::LieAlgebra
    highest_weight::WeightLatticeElem

    # The following fields are not set by default, just for caching
    dim::ZZRingElem
    character::Dict{WeightLatticeElem, ZZRingElem}

    function SimpleModuleData(L::LieAlgebra, highest_weight::WeightLatticeElem)
        new(L, highest_weight)
    end
end

function base_lie_algebra(V::SimpleModuleData)
    return V.L
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

mutable struct DemazureModuleData <: ModuleData
    L::LieAlgebra
    highest_weight::WeightLatticeElem
    weyl_group_elem::WeylGroupElem

    # The following fields are not set by default, just for caching
    dim::ZZRingElem
    character::Dict{WeightLatticeElem, ZZRingElem}

    function DemazureModuleData(L::LieAlgebra, highest_weight::WeightLatticeElem, weyl_group_elem::WeylGroupElem)
        new(L, highest_weight, weyl_group_elem)
    end
end

function base_lie_algebra(V::DemazureModuleData)
    return V.L
end

function highest_weight(V::DemazureModuleData)
    return V.highest_weight
end

function character(V::DemazureModuleData)
    if !isdefined(V, :character)
        V.character = demazure_character(ZZRingElem, base_lie_algebra(V), highest_weight(V), V.weyl_group_elem)
    end
    return V.character
end

function dim(V::DemazureModuleData)
    if !isdefined(V, :dim)
        V.dim = sum(values(character(V)); init=zero(ZZ))
    end
    return V.dim
end
