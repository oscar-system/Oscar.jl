@attributes mutable struct LieAlgebraExteriorPowerModule{C <: RingElement} <: LieAlgebraModule{C}
    inner_mod::LieAlgebraModule{C}
    power::Int
    ind_map::Vector{Vector{Int}}
    transformation_matrix_cache::Vector{Union{Nothing, <:MatElem{C}}}

    function LieAlgebraExteriorPowerModule{C}(
        inner_mod::LieAlgebraModule{C},
        power::Int;
        cached::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(
            LieAlgebraExteriorPowerModuleDict,
            (inner_mod, power),
            cached,
        ) do
            ind_map = collect(Combinatorics.combinations(1:dim(inner_mod), power))
            transformation_matrix_cache = Vector{Union{Nothing, <:MatElem{C}}}(
                nothing,
                dim(base_liealgebra(inner_mod)),
            )
            new{C}(inner_mod, power, ind_map, transformation_matrix_cache)
        end::LieAlgebraExteriorPowerModule{C}
    end
end

const LieAlgebraExteriorPowerModuleDict = CacheDictType{Tuple{LieAlgebraModule, Int}, LieAlgebraExteriorPowerModule}()

struct LieAlgebraExteriorPowerModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraExteriorPowerModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraExteriorPowerModuleElem{C}}) where {C <: RingElement} = LieAlgebraExteriorPowerModule{C}

elem_type(::Type{LieAlgebraExteriorPowerModule{C}}) where {C <: RingElement} = LieAlgebraExteriorPowerModuleElem{C}

parent(v::LieAlgebraExteriorPowerModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement} = base_ring(V.inner_mod)

base_liealgebra(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement} = base_liealgebra(V.inner_mod)

@attr dim(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement} = binomial(dim(V.inner_mod), V.power)


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement}
    print(io, "$(V.power)-th exterior power of ")
    print(IOContext(io, :compact => true), V.inner_mod)
end

function symbols(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement}
    if V.power == 1
        return symbols(V.inner_mod)
    end
    if isa(V.inner_mod, LieAlgebraStdModule)
        parentheses = identity
    else
        parentheses = x -> "($x)"
    end

    return [Symbol(join(s .|> parentheses, " ∧ ")) for s in Combinatorics.combinations(symbols(V.inner_mod), V.power)]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraExteriorPowerModule{C})(
    a::Vector{T},
) where {T <: LieAlgebraModuleElem{C}} where {C <: RingElement}
    length(a) == V.power || error("Length of vector does not match tensor power.")
    all(x -> parent(x) == V.inner_mod, a) || error("Incompatible modules.")
    mat = zero_matrix(base_ring(V), 1, dim(V))
    for (i, _inds) in enumerate(V.ind_map), inds in Combinatorics.permutations(_inds)
        sgn = levicivita(sortperm(inds))
        mat[1, i] += sgn * prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
    return LieAlgebraExteriorPowerModuleElem{C}(V, mat)
end


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_by_basisindex(V::LieAlgebraExteriorPowerModule{C}, i::Int) where {C <: RingElement}
    if V.transformation_matrix_cache[i] === nothing
        T = tensor_power(V.inner_mod, V.power)
        xT = transformation_matrix_by_basisindex(T, i)

        basis_change_E2T = zero(xT, dim(T), dim(V))
        basis_change_T2E = zero(xT, dim(V), dim(T))

        for (i, _inds) in enumerate(V.ind_map), inds in Combinatorics.permutations(_inds)
            sgn = levicivita(sortperm(inds))
            j = findfirst(==(inds), T.ind_map)
            basis_change_E2T[j, i] = sgn // factorial(V.power)
            basis_change_T2E[i, j] = sgn
        end

        V.transformation_matrix_cache[i] = basis_change_T2E * xT * basis_change_E2T
    end
    return (V.transformation_matrix_cache[i])::dense_matrix_type(C)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function exterior_power(V::LieAlgebraModule{C}, k::Int; cached::Bool=true) where {C <: RingElement}
    return LieAlgebraExteriorPowerModule{C}(V, k; cached)
end
