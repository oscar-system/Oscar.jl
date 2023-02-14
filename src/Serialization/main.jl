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
const oscarSerializationVersion = get_version_info()


################################################################################
# (De|En)coding types
#
const reverseTypeMap = Dict{String, Type}()

function registerSerializationType(@nospecialize(T::Type), str::String)
  if haskey(reverseTypeMap, str) && reverseTypeMap[str] != T
    error("encoded type $str already registered for a different type: $T versus $(reverseTypeMap[str])")
  end
  reverseTypeMap[str] = T
end

# registerSerializationType is a macro to ensure that the string we generate
# matches exactly the expression passed as first argument, and does not change
# in unexpected ways when import/export statements are adjusted.
macro registerSerializationType(ex::Any, str::Union{String,Nothing} = nothing)
  if str === nothing
    str = string(ex)
  end
  return esc(quote
    registerSerializationType($ex, $str)
    encodeType(::Type{<:$ex}) = $str
  end)
end

function encodeType(::Type{T}) where T
    error("unsupported type '$T' for encoding")
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
        @warn "Serialization: Generic Decoding of type $input"
        eval(Meta.parse(input))

        error("unsupported type '$input' for decoding")
    end
end



################################################################################
# High level
function save_type_dispatch(s::SerializerState, obj::T) where T
    if is_basic_serialization_type(T) && s.depth != 0
        return save_internal(s, obj)
    end

    if !isprimitivetype(T) && !Base.issingletontype(T)
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
        result[:_ns] = oscarSerializationVersion
    end
    return result
end

function load_type_dispatch(s::DeserializerState,
                            ::Type{T}, str::String; parent=nothing) where T
    @assert is_basic_serialization_type(T)
    if parent !== nothing
        load_internal_with_parent(s, T, str, parent)
    end
    return load_internal(s, T, str)
end

function load_type_dispatch(s::DeserializerState, ::Type{T}, dict::Dict;
                            parent=nothing) where T
    # File version to be dealt with on first breaking change
    # A file without version number is treated as the "first" version
    if dict[:type] == string(backref_sym)
        backref = s.objs[UUID(dict[:id])]
        backref isa T || throw(ErrorException("Backref of incorrect type encountered: $backref !isa $T"))
        return backref
    end

    # Decode the stored type, and compare it to the type `T` supplied by the caller.
    # If they are identical, just proceed. If not, then we assume that either
    # `T` is concrete, in which case `T <: U` should hold; or else `U` is
    # concrete, and `U <: T` should hold.
    #
    # However, we actually do not currently check for the types being concrete,
    # to allow for things like decoding `Vector{Vector}` ... we can tighten or loosen
    # these checks later on, depending on what we actually need...
    U = decodeType(dict[:type])
    U <: T || U >: T || throw(ErrorException("Type in file doesn't match target type: $(dict[:type]) not a subtype of $T"))

    Base.issingletontype(T) && return T()

    if parent !== nothing
        result = load_internal_with_parent(s, T, dict[:data], parent)
    else
        result = load_internal(s, T, dict[:data])
    end

    if haskey(dict, :id)
        s.objs[UUID(dict[:id])] = result
    end
    return result
end

function load_unknown_type(s::DeserializerState,
                           dict::Dict;
                           parent=nothing)
    if dict[:type] == string(backref_sym)
        return s.objs[UUID(dict[:id])]
    end

    T = decodeType(dict[:type])
    Base.issingletontype(T) && return T()
    return load_type_dispatch(s, T, dict; parent=parent)
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
# Include serialization implementations for various types

include("basic_types.jl")
include("containers.jl")
include("PolyhedralGeometry.jl")
include("Combinatorics.jl")
include("Fields.jl")
include("ToricGeometry.jl")
include("Rings.jl")
include("polymake.jl")
include("TropicalGeometry.jl")
include("QuadForm.jl")

################################################################################
# Include upgrade scripts

include("upgrades/main.jl")

function get_file_version(dict::Dict{Symbol, Any})
    ns = dict[:_ns]
    version_dict = ns[:Oscar][2]
    return get_version_number(version_dict)
end

function get_version_number(dict::Dict{Symbol, Any})
    return VersionNumber(dict[:major], dict[:minor], dict[:patch])
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
julia> save("/tmp/fourtitwo.json", 42);

julia> load("/tmp/fourtitwo.json")
42
```
"""
function save(io::IO, obj::Any)
    state = SerializerState()
    jsoncompatible = save_type_dispatch(state, obj)
    jsonstr = json(jsoncompatible)
    write(io, jsonstr)
    return nothing
end

function save(filename::String, obj::Any)
    open(filename, "w") do file
        save(file, obj)
    end
end

"""
    load(io::IO; parent::Any = nothing, type::Any = nothing)
    load(filename::String; parent::Any = nothing, type::Any = nothing)

Load the object stored in the given io stream
respectively in the file `filename`.

If `parent` is specified, then the root object of the loaded data
either will have this as its parent, or it is a container such as
`Vector` and `Tuple`, then the parent of the entries will be set to `parent`.

If a type `T` is given then attempt to load the root object of the data
being loaded with this type; if this fails, an error is thrown.

See [`save`](@ref).

# Examples

```jldoctest
julia> save("/tmp/fourtitwo.json", 42);

julia> load("/tmp/fourtitwo.json")
42

julia> load("/tmp/fourtitwo.json"; type=Int64)
42

julia> R, x = QQ["x"]
(Univariate Polynomial Ring in x over Rational Field, x)

julia> p = x^2 - x + 1
x^2 - x + 1

julia> save("/tmp/p.json", p)

julia> p_loaded = load("/tmp/p.json", parent=R)
x^2 - x + 1

julia> parent(p_loaded) === R
true

julia> save("/tmp/p_v.json", [p, p])

julia> loaded_p_v = load("/tmp/p_v.json", parent=R)
2-element Vector{fmpq_poly}:
 x^2 - x + 1
 x^2 - x + 1

julia> parent(loaded_p_v[1]) === parent(loaded_p_v[2]) === R
true
```
"""
function load(io::IO; parent::Any = nothing, type::Any = nothing)
    state = DeserializerState()
    # Check for type of file somewhere here?
    jsondict = JSON.parse(io, dicttype=Dict{Symbol, Any})

    haskey(jsondict, :_ns) || throw(ArgumentError("Namespace is missing"))
    _ns = jsondict[:_ns]
    if haskey(_ns, :polymake)
        # If this is a polymake file
        return load_from_polymake(jsondict)
    end
    haskey(_ns, :Oscar) || throw(ArgumentError("Not an Oscar object"))

    file_version = get_file_version(jsondict)

    if file_version < VERSION_NUMBER
        jsondict = upgrade(jsondict, file_version)
    end

    if type !== nothing
        return load_type_dispatch(state, type, jsondict; parent=parent)
    end
    return load_unknown_type(state, jsondict; parent=parent)
end

function load(filename::String; parent::Any = nothing, type::Any = nothing)
    open(filename) do file
        return load(file; parent=parent, type=type)
    end
end
