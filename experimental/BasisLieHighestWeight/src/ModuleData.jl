abstract type ModuleData end

# To be implemented by subtypes:
# Mandatory:
#   dim(M::MyModuleData) -> ZZRingElem
#   character(M::MyModuleData) -> Dict{WeightLatticeElem,<:IntegerUnion}


mutable struct SimpleModuleData <: ModuleData
    L::LieAlgebra
    highest_weight::WeightLatticeElem

    # The following fields are not set by default just set for caching
    dim::ZZRingElem
    character::Dict{WeightLatticeElem,<:IntegerUnion}

    function SimpleModuleData(L::LieAlgebra, highest_weight::WeightLatticeElem)
        new(L, highest_weight)
    end
end

function L(M::SimpleModuleData)
    return M.L
end

function highest_weight(M::SimpleModuleData)
    return M.highest_weight
end

function dim(M::SimpleModuleData)
    if !isdefined(M, :dim)
        M.dim = dim_of_simple_module(ZZRingElem, M.L, M.highest_weight)
    end
    return M.dim
end

function character(M::SimpleModuleData)
    if !isdefined(M, :character)
        M.character = LieAlgebraModule.character(ZZRingElem, M.L, M.highest_weight)
    end
    return M.character
end

mutable struct DemazureModuleData <: ModuleData
    L::LieAlgebra
    highest_weight::WeightLatticeElem
    weyl_group_elem::WeylGroupElem

    # The following fields are not set by default just set for caching
    dim::ZZRingElem
    character::Dict{WeightLatticeElem,<:IntegerUnion}

    function DemazureModuleData(L::LieAlgebra, highest_weight::WeightLatticeElem, weyl_group_elem::WeylGroupElem)
        new(L, highest_weight, weyl_group_elem)
    end
end

function L(M::SimDemazuModuleData)
    return M.L
end

function highest_weight(M::DemazureModuleData)
    return M.highest_weight
end

function character(M::DemazureModuleData)
    if !isdefined(M, :character)
        M.character = demazure_character(ZZRingElem, M.L, M.highest_weight, M.weyl_group_elem)
    end
    return M.character
end

function dim(M::DemazureModuleData)
    if !isdefined(M, :dim)
        if !isdefined(M, :character)
            character(M)
        else
            dim = 0
            for (w, d) in M.character
                dim += d
            end
            M.dim = ZZRingElem(dim)
        end
    end
    return M.dim
end