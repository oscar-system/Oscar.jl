@attributes mutable struct LieAlgebraAbstractModule{C <: RingElement} <: LieAlgebraModule{C}
    L::LieAlgebra{C}
    dim::Int
    transformation_matrices::Vector{MatElem{C}}
    s::Vector{Symbol}

    function LieAlgebraAbstractModule{C}(
        L::LieAlgebra{C},
        dimV::Int,
        transformation_matrices::Vector{<:MatElem{C}},
        s::Vector{Symbol};
        cached::Bool=true,
        check::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(
            LieAlgebraAbstractModuleDict,
            (L, dimV, transformation_matrices, s),
            cached,
        ) do
            all(m -> size(m) == (dimV, dimV), transformation_matrices) ||
                error("Invalid transformation matrix dimensions.")
            dimV == length(s) || error("Invalid number of basis element names.")

            V = new{C}(L, dimV, transformation_matrices, s)
            if check
                for xi in gens(L), xj in gens(L), v in gens(V)
                    bracket(xi, xj) * v == xi * (xj * v) - xj * (xi * v) ||
                        error("Structure constants do not define a module.")
                end
            end
            V
        end::LieAlgebraAbstractModule{C}
    end

end

const LieAlgebraAbstractModuleDict =
    CacheDictType{Tuple{LieAlgebra, Int, Vector{MatElem}, Vector{Symbol}}, LieAlgebraAbstractModule}()

struct LieAlgebraAbstractModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraAbstractModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraAbstractModuleElem{C}}) where {C <: RingElement} = LieAlgebraAbstractModule{C}

elem_type(::Type{LieAlgebraAbstractModule{C}}) where {C <: RingElement} = LieAlgebraAbstractModuleElem{C}

parent(v::LieAlgebraAbstractModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraAbstractModule{C}) where {C <: RingElement} = base_ring(base_liealgebra(V))

base_liealgebra(V::LieAlgebraAbstractModule{C}) where {C <: RingElement} = V.L

dim(V::LieAlgebraAbstractModule{C}) where {C <: RingElement} = V.dim


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraAbstractModule{C}) where {C <: RingElement}
    print(io, "AbstractModule of ")
    print(IOContext(io, :compact => true), base_liealgebra(V))
end

function symbols(V::LieAlgebraAbstractModule{C}) where {C <: RingElement}
    return V.s
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

function transformation_matrix_by_basisindex(V::LieAlgebraAbstractModule{C}, i::Int) where {C <: RingElement}
    return (V.transformation_matrices[i])::dense_matrix_type(C)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function abstract_module(
    L::LieAlgebra{C},
    dimV::Int,
    transformation_matrices::Vector{<:MatElem{C}},
    s::Vector{<:Union{AbstractString, Char, Symbol}}=[Symbol("v_$i") for i in 1:dimV];
    cached::Bool=true,
    check::Bool=true,
) where {C <: RingElement}
    return LieAlgebraAbstractModule{C}(L, dimV, transformation_matrices, Symbol.(s); cached, check=check)
end

function abstract_module(
    L::LieAlgebra{C},
    dimV::Int,
    struct_consts::Matrix{SRow{C}},
    s::Vector{<:Union{AbstractString, Char, Symbol}}=[Symbol("v_$i") for i in 1:dimV];
    cached::Bool=true,
    check::Bool=true,
) where {C <: RingElement}
    dim(L) == size(struct_consts, 1) || error("Invalid structure constants dimensions.")
    dimV == size(struct_consts, 2) || error("Invalid structure constants dimensions.")
    dimV == length(s) || error("Invalid number of basis element names.")

    transformation_matrices = [zero_matrix(base_ring(L), dimV, dimV) for _ in 1:dim(L)]
    for i in 1:dim(L), j in 1:dimV
        transformation_matrices[i][:, j] = transpose(dense_row(struct_consts[i, j], dimV))
    end

    return LieAlgebraAbstractModule{C}(L, dimV, transformation_matrices, Symbol.(s); cached, check=check)
end


function highest_weight_module(L::LieAlgebra{C}, weight::Vector{Int}; cached::Bool=true) where {C <: RingElement}
    struct_consts = liealgebra_highest_weight_module_struct_consts_gap(L, weight)
    dimV = size(struct_consts, 2)
    V = abstract_module(L, dimV, struct_consts; cached, check=false)
    set_attribute!(V, :highest_weight, weight)
    return V
end