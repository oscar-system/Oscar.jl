@attributes mutable struct LieAlgebraStdModule{C <: RingElement} <: LieAlgebraModule{C}
    L::LinearLieAlgebra{C}

    function LieAlgebraStdModule{C}(L::LinearLieAlgebra{C}; cached::Bool=true) where {C <: RingElement}
        return get_cached!(LieAlgebraStdModuleDict, L, cached) do
            new{C}(L)
        end::LieAlgebraStdModule{C}
    end
end

const LieAlgebraStdModuleDict = CacheDictType{LinearLieAlgebra, LieAlgebraStdModule}()

struct LieAlgebraStdModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraStdModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraStdModuleElem{C}}) where {C <: RingElement} = LieAlgebraStdModule{C}

elem_type(::Type{LieAlgebraStdModule{C}}) where {C <: RingElement} = LieAlgebraStdModuleElem{C}

parent(v::LieAlgebraStdModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraStdModule{C}) where {C <: RingElement} = base_ring(base_liealgebra(V))

base_liealgebra(V::LieAlgebraStdModule{C}) where {C <: RingElement} = V.L

@attr dim(V::LieAlgebraStdModule{C}) where {C <: RingElement} = base_liealgebra(V).n


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraStdModule{C}) where {C <: RingElement}
    print(io, "StdModule of ")
    print(IOContext(io, :compact => true), base_liealgebra(V))
end

function symbols(V::LieAlgebraStdModule{C}) where {C <: RingElement}
    return [Symbol("v_$(i)") for i in 1:dim(V)]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_by_basisindex(V::LieAlgebraStdModule{C}, i::Int) where {C <: RingElement}
    return matrix_repr_basis(base_liealgebra(V), i)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function standard_module(L::LinearLieAlgebra{C}; cached::Bool=true) where {C <: RingElement}
    return LieAlgebraStdModule{elem_type(base_ring(L))}(L; cached)
end
