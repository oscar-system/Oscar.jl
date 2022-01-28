using JSON

# struct which tracks state for (de)serialization
struct SerializerState
    # dict to track already serialized objects
    objmap::IdDict{Any,Int}

    # TODO: if we don't want to produce intermediate dictionaries (which is a lot of overhead), we would probably store an `IO` object here
    # io::IO
end

function SerializerState()
    return SerializerState(IdDict{Any, Int}())
end

struct DeserializerState
    objs::Vector{Any}  # or perhaps Dict{Int,Any} to be resilient against corrupts/malicious files using huge ids
end

function DeserializerState()
    return DeserializerState(Vector{Any}())
end


const backref_sym = Symbol("#backref")
function backref(ref::Int)
   return Dict(
     :type=>backref_sym,
     :version=>1, # ???
     :id=>ref,
     )
end

################################################################################
# Some setup for encoding types as strings and back and for getting the version
# info.
const typeMap = Dict{Type, String}([
    Int => "Base.Int",
    Cone => "Oscar.Cone",
    Polyhedron => "Oscar.Polyhedron",
    PolyhedralFan => "Oscar.PolyhedralFan",
    PolyhedralComplex => "Oscar.PolyhedralComplex",
    SubdivisionOfPoints => "Oscar.SubdivisionOfPoints",
    Graphs.Graph{Graphs.Undirected} => "Oscar.Graphs.Graph{Oscar.Graphs.Undirected}",
    Graphs.Graph{Graphs.Directed} => "Oscar.Graphs.Graph{Oscar.Graphs.Directed}",
    SimplicialComplex => "Oscar.SimplicialComplex",
    LinearProgram => "Oscar.LinearProgram"
])

const reverseTypeMap = Dict{String, Type}(value => key for (key, value) in typeMap)

function get_version_info()
    result = Dict{Symbol, Any}(
        :Oscar => ["https://github.com/oscar-system/Oscar.jl", VERSION_NUMBER],
        :Julia => ["https://julialang.org/", VERSION]
    )
    return result
end
const versionInfo = get_version_info()


################################################################################
# High level
function save(s::SerializerState, obj::T) where T
    if !isprimitivetype(T) && !Base.issingletontype(T)
        # if obj is already serialzed, just output
        # a backref
        ref = get(s.objmap, obj, nothing)
        ref !== nothing && return backref(ref)
        # otherwise, 
        ref = s.objmap[obj] = length(s.objmap)
    else
        ref = -1
    end
    # invoke the actual serializer
    typestr = typeMap[T]
    return Dict(
        :_ns => versionInfo,
        :type => typestr,
        :id => ref,
        :data => save_intern(s, obj)
    )
end


function load(s::DeserializerState, dict::Dict)
    haskey(dict, :_ns) || throw(ArgumentError("Namespace is missing"))
    _ns = dict[:_ns]
    if haskey(_ns, :polymake)
        # If this is a polymake file
        return load_from_polymake(dict)
    end
    haskey(_ns, :Oscar) || throw(ArgumentError("Not an Oscar object"))
    typestr = dict[:type]
    id = dict[:id]
    # TODO: deal with versions? enforce their presence?

    if typestr == backref_sym
        return s.objs[id]
    end

    T = reverseTypeMap[typestr]

    Base.issingletontype(T) && return T()

    # TODO: offer a generic handler for primitive
    # types e.g. storing them as byte strings or so
    #Base.isprimitivetype(T) && return load_primitivetype(T, dict) #

    return load_intern(s, T, dict[:data])
end

################################################################################
# Interacting with files
@doc Markdown.doc"""
    save(obj::T, filename::String) where T

Save an object `T` to a file.
"""
function save(obj::T, filename::String) where T
    state = SerializerState()
    jsoncompatible = save(state, obj)
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
    return load(state, jsondict)
end





include("PolyhedralGeometry.jl")
include("Combinatorics.jl")

