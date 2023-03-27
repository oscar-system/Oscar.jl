@attributes mutable struct LieAlgebraSymmetricPowerModule{C <: RingElement} <: LieAlgebraModule{C}
    inner_mod::LieAlgebraModule{C}
    power::Int
    ind_map::Vector{Vector{Int}}
    transformation_matrix_cache::Vector{Union{Nothing, <:MatElem{C}}}

    function LieAlgebraSymmetricPowerModule{C}(
        inner_mod::LieAlgebraModule{C},
        power::Int;
        cached::Bool=true,
    ) where {C <: RingElement}
        return get_cached!(
            LieAlgebraSymmetricPowerModuleDict,
            (inner_mod, power),
            cached,
        ) do
            ind_map = collect(
                Combinatorics.with_replacement_combinations(1:dim(inner_mod), power),
            )
            transformation_matrix_cache = Vector{Union{Nothing, <:MatElem{C}}}(
                nothing,
                dim(base_liealgebra(inner_mod)),
            )
            new{C}(inner_mod, power, ind_map, transformation_matrix_cache)
        end::LieAlgebraSymmetricPowerModule{C}
    end
end

const LieAlgebraSymmetricPowerModuleDict = CacheDictType{Tuple{LieAlgebraModule, Int}, LieAlgebraSymmetricPowerModule}()

struct LieAlgebraSymmetricPowerModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraSymmetricPowerModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraSymmetricPowerModuleElem{C}}) where {C <: RingElement} = LieAlgebraSymmetricPowerModule{C}

elem_type(::Type{LieAlgebraSymmetricPowerModule{C}}) where {C <: RingElement} = LieAlgebraSymmetricPowerModuleElem{C}

parent(v::LieAlgebraSymmetricPowerModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement} = base_ring(V.inner_mod)

base_liealgebra(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement} = base_liealgebra(V.inner_mod)

@attr dim(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement} =
    binomial(dim(V.inner_mod) + V.power - 1, V.power)


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement}
    print(io, "$(V.power)-th symmetric power of ")
    print(IOContext(io, :compact => true), V.inner_mod)
end

function symbols(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement}
    if V.power == 1
        return symbols(V.inner_mod)
    end
    if isa(V.inner_mod, LieAlgebraStdModule)
        parentheses = identity
    else
        parentheses = x -> "($x)"
    end

    return [
        Symbol(
            join((
                begin
                    e = count(==(i), inds)
                    if e == 1
                        s |> parentheses
                    else
                        "$(s |> parentheses)^$e"
                    end
                end for (i, s) in enumerate(symbols(V.inner_mod)) if in(i, inds)
            ), "*"),
        ) for inds in Combinatorics.with_replacement_combinations(1:dim(V.inner_mod), V.power)
    ]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraSymmetricPowerModule{C})(
    a::Vector{T},
) where {T <: LieAlgebraModuleElem{C}} where {C <: RingElement}
    length(a) == V.power || error("Length of vector does not match tensor power.")
    all(x -> parent(x) == V.inner_mod, a) || error("Incompatible modules.")
    mat = zero_matrix(base_ring(V), 1, dim(V))
    for (i, _inds) in enumerate(V.ind_map), inds in unique(Combinatorics.permutations(_inds))
        mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
    return LieAlgebraSymmetricPowerModuleElem{C}(V, mat)
end


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_by_basisindex(V::LieAlgebraSymmetricPowerModule{C}, i::Int) where {C <: RingElement}
    if V.transformation_matrix_cache[i] === nothing
        T = tensor_power(V.inner_mod, V.power)
        xT = transformation_matrix_by_basisindex(T, i)

        basis_change_S2T = zero(xT, dim(T), dim(V))
        basis_change_T2S = zero(xT, dim(V), dim(T))
        for (i, _inds) in enumerate(V.ind_map), inds in Combinatorics.permutations(_inds)
            j = findfirst(==(inds), T.ind_map)
            basis_change_S2T[j, i] += 1 // factorial(V.power)
            basis_change_T2S[i, j] = 1
        end

        V.transformation_matrix_cache[i] = basis_change_T2S * xT * basis_change_S2T
    end
    return (V.transformation_matrix_cache[i])::dense_matrix_type(C)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function symmetric_power(V::LieAlgebraModule{C}, k::Int; cached::Bool=true) where {C <: RingElement}
    return LieAlgebraSymmetricPowerModule{C}(V, k; cached)
end
