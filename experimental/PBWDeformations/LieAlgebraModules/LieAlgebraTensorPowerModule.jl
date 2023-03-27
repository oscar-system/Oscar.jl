@attributes mutable struct LieAlgebraTensorPowerModule{C <: RingElement} <: LieAlgebraModule{C}
    inner_mod::LieAlgebraModule{C}
    power::Int
    ind_map::Vector{Vector{Int}}
    transformation_matrix_cache::Vector{Union{Nothing, <:MatElem{C}}}

    function LieAlgebraTensorPowerModule{C}(
        inner_mod::LieAlgebraModule{C},
        power::Int;
        cached::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(
            LieAlgebraTensorPowerModuleDict,
            (inner_mod, power),
            cached,
        ) do
            ind_map = collect(ProductIterator(1:dim(inner_mod), power)) .|> reverse
            transformation_matrix_cache = Vector{Union{Nothing, <:MatElem{C}}}(
                nothing,
                dim(base_liealgebra(inner_mod)),
            )
            new{C}(inner_mod, power, ind_map, transformation_matrix_cache)
        end::LieAlgebraTensorPowerModule{C}
    end
end

const LieAlgebraTensorPowerModuleDict = CacheDictType{Tuple{LieAlgebraModule, Int}, LieAlgebraTensorPowerModule}()
struct LieAlgebraTensorPowerModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraTensorPowerModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraTensorPowerModuleElem{C}}) where {C <: RingElement} = LieAlgebraTensorPowerModule{C}

elem_type(::Type{LieAlgebraTensorPowerModule{C}}) where {C <: RingElement} = LieAlgebraTensorPowerModuleElem{C}

parent(v::LieAlgebraTensorPowerModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraTensorPowerModule{C}) where {C <: RingElement} = base_ring(V.inner_mod)

base_liealgebra(V::LieAlgebraTensorPowerModule{C}) where {C <: RingElement} = base_liealgebra(V.inner_mod)

@attr dim(V::LieAlgebraTensorPowerModule{C}) where {C <: RingElement} = dim(V.inner_mod)^V.power


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraTensorPowerModule{C}) where {C <: RingElement}
    print(io, "$(V.power)-th tensor power of ")
    print(IOContext(io, :compact => true), V.inner_mod)
end

function symbols(V::LieAlgebraTensorPowerModule{C}) where {C <: RingElement}
    if V.power == 1
        return symbols(V.inner_mod)
    end
    if isa(V.inner_mod, LieAlgebraStdModule)
        parentheses = identity
    else
        parentheses = x -> "($x)"
    end

    return [Symbol(join(s .|> parentheses, " ⊗ ")) for s in ProductIterator(symbols(V.inner_mod), V.power) .|> reverse]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraTensorPowerModule{C})(a::Vector{T}) where {T <: LieAlgebraModuleElem{C}} where {C <: RingElement}
    length(a) == V.power || error("Length of vector does not match tensor power.")
    all(x -> parent(x) == V.inner_mod, a) || error("Incompatible modules.")
    mat = zero_matrix(base_ring(V), 1, dim(V))
    for (i, inds) in enumerate(V.ind_map)
        mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
    return LieAlgebraTensorPowerModuleElem{C}(V, mat)
end


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_by_basisindex(V::LieAlgebraTensorPowerModule{C}, i::Int) where {C <: RingElement}
    if V.transformation_matrix_cache[i] === nothing
        y = transformation_matrix_by_basisindex(V.inner_mod, i)::dense_matrix_type(C)
        V.transformation_matrix_cache[i] =
            sum(reduce(kronecker_product, (j == i ? y : one(y) for j in 1:V.power)) for i in 1:V.power)
    end
    return (V.transformation_matrix_cache[i])::dense_matrix_type(C)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function tensor_power(V::LieAlgebraModule{C}, k::Int; cached::Bool=true) where {C <: RingElement}
    return LieAlgebraTensorPowerModule{C}(V, k; cached)
end
