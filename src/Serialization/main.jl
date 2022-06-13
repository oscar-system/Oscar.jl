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
const typeMap = Dict{Type, String}()
const reverseTypeMap = Dict{String, Type}()

function registerSerializationType(@nospecialize(T::Type), str::String)
  isconcretetype(T) || error("Currently only concrete types can be registered, but $T is not")
  if haskey(typeMap, T) && typeMap[T] != str
    error("type $T already registered with a different encoding: $str versus $(typeMap[T])")
  end

  if haskey(reverseTypeMap, str) && reverseTypeMap[str] != T
    error("encoded type $str already registered for a different type: $T versus $(reverseTypeMap[str])")
  end
  typeMap[T] = str
  reverseTypeMap[str] = T
end

# registerSerializationType is a macro to ensure that the string we generate
# matches exactly the expression passed as first argument, and does not change
# in unexpected ways when import/export statements are adjusted.
macro registerSerializationType(ex::Any, str::Union{String,Nothing} = nothing)
  if str === nothing
    str = string(ex)
  end
  return :( registerSerializationType($ex, $str) )
end

for (T, str) in (
    AbstractAlgebra.Generic.Frac{fmpq_poly} => "AbstractAlgebra.Generic.Frac{fmpq_poly}",
    AbstractAlgebra.Generic.FracField{fmpq_poly} => "AbstractAlgebra.Generic.FracField{fmpq_poly}",
    AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_poly}} => "AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_poly}}",
    AbstractAlgebra.Generic.MPoly{Hecke.NfRelElem{nf_elem}} => "AbstractAlgebra.Generic.MPoly{Hecke.NfRelElem{nf_elem}}",
    AbstractAlgebra.Generic.MPoly{Hecke.NfRelNSElem{nf_elem}} => "AbstractAlgebra.Generic.MPoly{Hecke.NfRelNSElem{nf_elem}}",
    AbstractAlgebra.Generic.MPoly{nf_elem} => "AbstractAlgebra.Generic.MPoly{nf_elem}",
    AbstractAlgebra.Generic.MPoly{NfAbsNSElem} => "AbstractAlgebra.Generic.MPoly{NfAbsNSElem}",
    AbstractAlgebra.Generic.MPolyRing{AbstractAlgebra.Generic.Frac{fmpq_poly}} => "AbstractAlgebra.Generic.MPolyRing{AbstractAlgebra.Generic.Frac{fmpq_poly}}",
    AbstractAlgebra.Generic.MPolyRing{Hecke.NfRelElem{nf_elem}} => "AbstractAlgebra.Generic.MPolyRing{Hecke.NfRelElem{nf_elem}}",
    AbstractAlgebra.Generic.MPolyRing{Hecke.NfRelNSElem{nf_elem}} => "AbstractAlgebra.Generic.MPolyRing{Hecke.NfRelNSElem{nf_elem}}",
    AbstractAlgebra.Generic.MPolyRing{nf_elem} => "AbstractAlgebra.Generic.MPolyRing{nf_elem}",
    AbstractAlgebra.Generic.MPolyRing{NfAbsNSElem} => "AbstractAlgebra.Generic.MPolyRing{NfAbsNSElem}",
    AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.Frac{fmpq_poly}} => "AbstractAlgebra.Generic.Poly{AbstractAlgebra.Generic.Frac{fmpq_poly}}",
    AbstractAlgebra.Generic.Poly{Hecke.NfRelElem{nf_elem}} => "AbstractAlgebra.Generic.Poly{Hecke.NfRelElem{nf_elem}}",
    AbstractAlgebra.Generic.Poly{Hecke.NfRelNSElem{nf_elem}} => "AbstractAlgebra.Generic.Poly{Hecke.NfRelNSElem{nf_elem}}",
    AbstractAlgebra.Generic.Poly{nf_elem} => "AbstractAlgebra.Generic.Poly{nf_elem}",
    AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.Frac{fmpq_poly}} => "AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.Frac{fmpq_poly}}",
    AbstractAlgebra.Generic.PolyRing{Hecke.NfRelElem{nf_elem}} => "AbstractAlgebra.Generic.PolyRing{Hecke.NfRelElem{nf_elem}}",
    AbstractAlgebra.Generic.PolyRing{Hecke.NfRelNSElem{nf_elem}} => "AbstractAlgebra.Generic.PolyRing{Hecke.NfRelNSElem{nf_elem}}",
    AbstractAlgebra.Generic.PolyRing{nf_elem} => "AbstractAlgebra.Generic.PolyRing{nf_elem}",
    Hecke.NfRel{nf_elem} => "Hecke.NfRel{nf_elem}",
    Hecke.NfRelElem{nf_elem} => "Hecke.NfRelElem{nf_elem}",
    Hecke.NfRelNSElem{nf_elem} => "Hecke.NfRelNSElem{nf_elem}",
    MPolyIdeal{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_poly}}} => "MPolyIdeal{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_poly}}}",
    MPolyIdeal{AbstractAlgebra.Generic.MPoly{Hecke.NfRelElem{nf_elem}}} => "MPolyIdeal{AbstractAlgebra.Generic.MPoly{Hecke.NfRelElem{nf_elem}}}",
    MPolyIdeal{AbstractAlgebra.Generic.MPoly{Hecke.NfRelNSElem{nf_elem}}} => "MPolyIdeal{AbstractAlgebra.Generic.MPoly{Hecke.NfRelNSElem{nf_elem}}}",
    MPolyIdeal{AbstractAlgebra.Generic.MPoly{nf_elem}} => "MPolyIdeal{AbstractAlgebra.Generic.MPoly{nf_elem}}",
    MPolyIdeal{fmpq_mpoly} => "MPolyIdeal{fmpq_mpoly}",
    MPolyIdeal{fmpz_mpoly} => "MPolyIdeal{fmpz_mpoly}",
    MPolyIdeal{fq_nmod_mpoly} => "MPolyIdeal{fq_nmod_mpoly}",
    MPolyIdeal{nmod_mpoly} => "MPolyIdeal{nmod_mpoly}",
    NfRelNS{nf_elem} => "NfRelNS{nf_elem}",
    Polymake.BigObjectAllocated => "Polymake.BigObject",
    Vector{AbstractAlgebra.Generic.Poly{nf_elem}} => "Vector{AbstractAlgebra.Generic.Poly{nf_elem}}",
    Vector{AbstractAlgebra.Ring} => "Vector{AbstractAlgebra.Ring}",
    Vector{Any} => "Vector{Any}",
    Vector{Float16} => "Vector{Float16}",
    Vector{Float32} => "Vector{Float32}",
    Vector{Float64} => "Vector{Float64}",
    Vector{fmpq_poly} => "Vector{fmpq_poly}",
    Vector{gfp_elem} => "Vector{gfp_elem}",
    Vector{gfp_fmpz_elem} => "Vector{gfp_fmpz_elem}",
    Vector{Int8} => "Vector{Int8}",
    Vector{Int16} => "Vector{Int16}",
    Vector{Int32} => "Vector{Int32}",
    Vector{Int64} => "Vector{Int64}",
    Vector{Int128} => "Vector{Int128}",
    Vector{LinearProgram{fmpq}} => "Vector{LinearProgram{fmpq}}",
    Vector{ToricDivisor} => "Vector{ToricDivisor}",
    Vector{UInt8} => "Vector{UInt8}",
    Vector{UInt16} => "Vector{UInt16}",
    Vector{UInt32} => "Vector{UInt32}",
    Vector{UInt64} => "Vector{UInt64}",
    Vector{UInt128} => "Vector{UInt128}",
    Vector{Union{LinearProgram, Polyhedron}} => "Vector{Union{LinearProgram, Polyhedron}}",
    )

  registerSerializationType(T, str)
end


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
        # entered here. (also computationally expensive)
        # Standard Oscar tests should never pass this line
        @warn "Serialization: Generic Decoding of Type: $input"
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
        if ref !== nothing
            return Dict{Symbol, Any}(
              :type => backref_sym,
              :id => string(ref),
              :version => 1, # ???
              )
        end
        # otherwise,
        ref = s.objmap[obj] = uuid4()
    else
        ref = nothing
    end
    
    result = Dict{Symbol, Any}(:type => encodeType(T))
    if ref !== nothing
        result[:id] = string(ref)
    end
    if !Base.issingletontype(T)
        s.depth += 1
        # invoke the actual serializer
        result[:data] = save_internal(s, obj)
        s.depth -= 1
    end
    if s.depth == 0
        result[:_ns] = versionInfo
    end
    return result
end

function load_type_dispatch(s::DeserializerState, ::Type{T}, dict::Dict) where T
    # File version to be dealt with on first breaking change
    # A file without version number is treated as the "first" version
    
    if dict[:type] == string(backref_sym)
        backref = s.objs[UUID(dict[:id])]
        backref isa T || throw(ErrorException("Backref of incorrect type encountered: $backref !isa $T"))
        return backref
    end

    encodeType(T) == dict[:type] || throw(ErrorException("Type in file doesn't match target type: $(dict[:type]) != $T"))

    Base.issingletontype(T) && return T()

    result = load_internal(s, T, dict[:data])
    if haskey(dict, :id)
        s.objs[UUID(dict[:id])] = result
    end
    return result
end

function load_unknown_type(s::DeserializerState, dict::Dict; check_namespace=false)
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

    return load_type_dispatch(s, T, dict)
end


################################################################################
# Default generic save_internal, load_internal
function save_internal_generic(s::SerializerState, obj::T) where T
    result = Dict{Symbol, Any}()
    for n in fieldnames(T)
        if n != :__attrs
            result[n] = save_type_dispatch(s, getfield(obj, n))
        end
    end
    return result
end

function load_internal_generic(s::DeserializerState, ::Type{T}, dict::Dict) where T
    fields = []
    for (n,t) in zip(fieldnames(T), fieldtypes(T))
        if n!= :__attrs
            push!(fields, load_type_dispatch(s, t, dict[n]))
        end
    end
    return T(fields...)
end


################################################################################
# Interacting with IO streams and files

"""
    save(io::IO, obj::Any)
    save(filename::String, obj::Any)

Save an object `T` to the given io stream
respectively to the file `filename`.

See [`load`](@ref).

# Examples

```jldoctest
julia> save("fourtitwo.json", 42);

julia> load("fourtitwo.json")
42
```
"""
function save(io::IO, obj::Any)
    state = SerializerState()
    jsoncompatible = save_type_dispatch(state, obj)
    jsonstr = json(jsoncompatible)
    write(io, jsonstr)
end

function save(filename::String, obj::Any)
    open(filename, "w") do file
        save(file, obj)
    end
end

"""
    load(io::IO)
    load(filename::String)

Load the object stored in the given io stream
respectively in the file `filename`.

See [`save`](@ref).

# Examples

```jldoctest
julia> save("fourtitwo.json", 42);

julia> load("fourtitwo.json")
42
```
"""
function load(io::IO)
    state = DeserializerState()
    # Check for type of file somewhere here?
    jsondict = JSON.parse(io, dicttype=Dict{Symbol, Any})
    return load_unknown_type(state, jsondict; check_namespace=true)
end

function load(filename::String)
    open(filename) do file
        return load(file)
    end
end

"""
    load(io::IO, ::Type)
    load(filename::String, T::Type)

Load the object of the given type stored in the given io stream
respectively in the file `filename`.

This guarantees that the end result has the given type.

See [`save`](@ref).

# Examples

```jldoctest
julia> save("fourtitwo.json", 42);

julia> load("fourtitwo.json")
42

julia> load("fourtitwo.json", String)
ERROR: Type in file doesn't match target type: Base.Int != String
```
"""
function load(io::IO, T::Type)
    state = DeserializerState()
    # Check for type of file somewhere here?
    jsondict = JSON.parse(io, dicttype=Dict{Symbol, Any})
    return load_type_dispatch(state, T, jsondict)
end

function load(filename::String, T::Type)
    open(filename) do file
        return load(file, T)
    end
end

include("basic_types.jl")
include("containers.jl")
include("PolyhedralGeometry.jl")
include("Combinatorics.jl")
include("Fields.jl")
include("ToricGeometry.jl")
include("Rings.jl")
include("polymake.jl")
