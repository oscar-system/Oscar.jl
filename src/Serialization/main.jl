using JSON
using UUIDs

# struct which tracks state for (de)serialization
mutable struct SerializerState
    # dict to track already serialized objects
    objmap::IdDict{Any, UUID}
    depth::Int

    # TODO: if we don't want to produce intermediate dictionaries (which is a lot of overhead), we would probably store an `IO` object here
    # io::IO
end

function SerializerState()
    return SerializerState(IdDict{Any, UUID}(), 0)
end

struct DeserializerState
    objs::Dict{UUID, Any}  # or perhaps Dict{Int,Any} to be resilient against corrupts/malicious files using huge ids
end

function DeserializerState()
    return DeserializerState(Dict{UUID, Any}())
end


const backref_sym = Symbol("#backref")
function backref(ref::UUID)
    return Dict(
        :type=>backref_sym,
        :version=>1, # ???
        :id=>string(ref),
    )
end

################################################################################
# Version info

function get_version_info()
    result = Dict{Symbol, Any}(
        :Oscar => ["https://github.com/oscar-system/Oscar.jl", VERSION_NUMBER]
    )
    return result
end
const versionInfo = get_version_info()


################################################################################
# (De|En)coding types
#
# Right now this is still quite primitive, however there is functionality for
# encoding Vectors without explicitly adding them to the typeMap. This will
# have to be adapted though to make it more generic for matrices, arrays, etc.
const typeMap = Dict{Type, String}([
    Int => "Base.Int",
    Graphs.Graph{Graphs.Undirected} => "Oscar.Graphs.Graph{Oscar.Graphs.Undirected}",
    Graphs.Graph{Graphs.Directed} => "Oscar.Graphs.Graph{Oscar.Graphs.Directed}",
    SimplicialComplex => "Oscar.SimplicialComplex",
    Vector => "Vector",
    Nemo.GaloisField => "Nemo.GaloisField",
    Nemo.GaloisFmpzField => "Nemo.GaloisFmpzField",
    fmpz => "fmpz",
    fmpq => "fmpq",
    gfp_elem => "Nemo.gfp_elem",
    gfp_fmpz_elem => "Nemo.gfp_fmpz_elem",
    fmpq_mpoly => "fmpq_mpoly",
    fmpz_mpoly => "fmpz_mpoly",
    Nemo.NmodRing => "Nemo.NmodRing",
    Hecke.NfRel{nf_elem} => "Hecke.NfRel{nf_elem}",
])

const reverseTypeMap = Dict{String, Type}(value => key for (key, value) in typeMap)


function encodeType(::Type{T}) where T
    if haskey(typeMap, T)
        return typeMap[T]
    else
        # As a default just save the type as a string.
        string(T)
    end
end



function decodeType(input::String)
    if haskey(reverseTypeMap, input)
        return reverseTypeMap[input]
    else
        # As a default, parse the type from the string.
        #
        # WARNING: Never deserialize data from an untrusted source, as this
        # parsing is insecure and potentially malicious code could be
        # entered here.
        eval(Meta.parse(input))
    end
end




################################################################################
# High level
function save_type_dispatch(s::SerializerState, obj::T) where T
    if !isprimitivetype(T) && !Base.issingletontype(T) && T !== Symbol
        # if obj is already serialzed, just output
        # a backref
        ref = get(s.objmap, obj, nothing)
        ref !== nothing && return backref(ref)
        # otherwise, 
        ref = s.objmap[obj] = uuid4()
    else
        ref = -1
    end
    # invoke the actual serializer
    encodedType = encodeType(T)
    s.depth += 1
    result = Dict{Symbol, Any}(
        :type => encodedType,
        :id => string(ref),
        :data => save_internal(s, obj),
    )
    s.depth -= 1
    if s.depth == 0
        result[:_ns] = versionInfo
    end
    return result
end

function load_type_dispatch(s::DeserializerState, ::Type{T}, dict::Dict) where T
    # TODO: deal with versions? enforce their presence?
    if dict[:type] == string(backref_sym)
        return s.objs[UUID(dict[:id])]
    end

    # TODO: compare T against decodedType ???
    #decodedType = decodeType(dict[:type])
    result = load_internal(s, T, dict[:data])
    id = dict[:id]
    
    if id != "-1"
        s.objs[UUID(id)] = result
    end
    
    return result
end

function load_type_dispatch(s::DeserializerState, dict::Dict; check_namespace=true)
    if check_namespace
        haskey(dict, :_ns) || throw(ArgumentError("Namespace is missing"))
        _ns = dict[:_ns]
        if haskey(_ns, :polymake)
            # If this is a polymake file
            return load_from_polymake(dict)
        end
        haskey(_ns, :Oscar) || throw(ArgumentError("Not an Oscar object"))
    end

    if dict[:type] == string(backref_sym)
        return s.objs[UUID(dict[:id])]
    end

    T = decodeType(dict[:type])
    Base.issingletontype(T) && return T()

    # TODO: offer a generic handler for primitive
    # types e.g. storing them as byte strings or so
    #Base.isprimitivetype(T) && return load_primitivetype(T, dict) #
    return load_type_dispatch(s, T, dict)
end


################################################################################
# Default generic save_internal, load_internal
function save_internal(s::SerializerState, obj::T) where T
    result = Dict{Symbol, Any}()
    for n in fieldnames(T)
        if n != :__attrs
            result[n] = save_type_dispatch(s, getfield(obj, n))
        end
    end
    return result
end

function load_internal(s::DeserializerState, ::Type{T}, dict::Dict) where T
    fields = []
    for (n,t) in zip(fieldnames(T), fieldtypes(T))
        if n!= :__attrs
            push!(fields, load_type_dispatch(s, t, dict[n]))
        end
    end
    return T(fields...)
end


################################################################################
# Interacting with files
@doc Markdown.doc"""
    save(obj::T, filename::String) where T

Save an object `T` to a file.
"""
function save(obj::T, filename::String) where T
    state = SerializerState()
    jsoncompatible = save_type_dispatch(state, obj)
    jsonstr = json(jsoncompatible)
    open(filename, "w") do file
        write(file, jsonstr)
    end
end

@doc Markdown.doc"""
    load(filename::String)

Load the object stored in a file.
"""
function load(filename::String)
    state = DeserializerState()
    # Check for type of file somewhere here?
    jsondict = JSON.parsefile(filename, dicttype=Dict{Symbol, Any})
    return load_type_dispatch(state, jsondict)
end




include("basic_types.jl")
include("containers.jl")
include("PolyhedralGeometry.jl")
include("Combinatorics.jl")
include("Fields.jl")
include("ToricGeometry.jl")
include("Rings.jl")
include("polymake.jl")


@deprecate save_cone(Obj::Cone, filename::String) save(Obj, filename)
@deprecate load_cone(filename::String) load(Obj, filename)
@deprecate save_simplicialcomplex(Obj::SimplicialComplex, filename::String) save(Obj, filename)
@deprecate load_simplicialcomplex(filename::String) load(Obj, filename)
@deprecate save_subdivisionofpoints(Obj::SubdivisionOfPoints, filename::String) save(Obj, filename)
@deprecate load_subdivisionofpoints(filename::String) load(Obj, filename)
@deprecate save_polyhedron(Obj::Polyhedron, filename::String) save(Obj, filename)
@deprecate load_polyhedron(filename::String) load(Obj, filename)
@deprecate save_polyhedralcomplex(Obj::PolyhedralComplex, filename::String) save(Obj, filename)
@deprecate load_polyhedralcomplex(filename::String) load(Obj, filename)
@deprecate save_polyhedralfan(Obj::PolyhedralFan, filename::String) save(Obj, filename)
@deprecate load_polyhedralfan(filename::String) load(Obj, filename)
@deprecate save_linearprogram(Obj::LinearProgram, filename::String) save(Obj, filename)
@deprecate load_linearprogram(filename::String) load(Obj, filename)
@deprecate save_graph(Obj::Graphs.Graph{Graphs.Directed}, filename::String) save(Obj, filename)
@deprecate save_graph(Obj::Graphs.Graph{Graphs.Undirected}, filename::String) save(Obj, filename)
@deprecate load_graph(filename::String) load(Obj, filename)
