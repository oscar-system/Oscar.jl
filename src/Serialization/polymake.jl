# Loading polymake objects into the appropriate types.


##############################################################################
"""
    load_from_polymake(::Type{Graphs.Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Graphs.Directed, Graphs.Undirected}}

Load a graph stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{Graphs.Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Graphs.Directed, Graphs.Undirected}}
    polymake_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    return Graphs.Graph{T}(polymake_object)
end


const polymake2OscarTypes = Dict{String, Type}([
    "polytope::Cone<Rational>" => Cone{fmpq},
    "polytope::Polytope<Rational>" => Polyhedron{fmpq},
    "fan::PolyhedralFan<Rational>" => PolyhedralFan{fmpq},
    "fan::PolyhedralComplex<Rational>" => PolyhedralComplex{fmpq},
    "fan::SubdivisionOfPoints<Rational>" => SubdivisionOfPoints{fmpq},
    "topaz::SimplicialComplex" => SimplicialComplex,
    "common::GraphAdjacency<Undirected>" => Graphs.Graph{Graphs.Undirected},
    "common::GraphAdjacency<Directed>" => Graphs.Graph{Graphs.Directed},
])


function load_from_polymake(::Type{T}, jsondict::Dict{Symbol, Any}) where {
        S<:scalar_types,
        T<:Union{Cone{S}, Polyhedron{S}, PolyhedralFan{S}, 
            PolyhedralComplex{S}, SubdivisionOfPoints{S}, SimplicialComplex}}
    
    inner_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    return T(inner_object)
end


# Distinguish between the various polymake datatypes.
function load_from_polymake(jsondict::Dict{Symbol, Any})
    typename = jsondict[:_type]
    if haskey(polymake2OscarTypes, typename)
        oscar_type = polymake2OscarTypes[typename]
        return load_from_polymake(oscar_type, jsondict)
    else 
        # We just try to default to something from Polymake.jl
        deserialized = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))

        try
            return convert_to_oscar(deserialized)
        catch e
            if e isa MethodError
                @warn "No function for converting the deserialized Polymake type to Oscar type"
                return deserialized
            else
                throw(e)
            end
        end
    end
end

function convert_to_oscar(p::Polymake.PolynomialAllocated{Polymake.Rational, Int64};
                          R::Union{MPolyRing, Nothing} = nothing)
    coeff_vec = convert(Vector{fmpq}, Polymake.coefficients_as_vector(p))
    monomials = Matrix{Int}(Polymake.monomials_as_matrix(p))
    n_vars = length(monomials[:, 1])
    # not sure if the numbering is the best choice but it matches Polymake
    if isnothing(R)
        R, _ = PolynomialRing(QQ, "x" => 0:n_vars - 1, cached=false)
    end

    return R(coeff_vec, [monomials[:, i] for i in 1:ncols(monomials)])
end

function convert_to_oscar(O::Polymake.BigObjectAllocated)
    big_object_name = Polymake.type_name(O)

    if "Ideal" == big_object_name
        n_vars = O.N_VARIABLES
        R, _ = PolynomialRing(QQ, "x" => 0:n_vars - 1, cached=false)
        converted_generators = map(p -> convert_to_oscar(p, R=R), O.GENERATORS)
        
        return ideal(R, converted_generators)
    else
        throw(MethodError)
    end
end
