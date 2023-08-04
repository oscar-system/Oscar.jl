using JSON
using UUIDs

include("serializers.jl")

const backref_sym = Symbol("#backref")
################################################################################
# Meta Data

@Base.kwdef struct MetaData
    author_orcid::Union{String, Nothing} = nothing
    name::Union{String, Nothing} = nothing
end

function metadata(;args...)
    return MetaData(;args...)
end

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
# It also sets the value of serialize_with_id, which determines
# whether or not the type can be back referenced.
# If omitted, the default is that no back references are allowed.
function registerSerializationType(ex::Any,
                                   uses_id::Bool,
                                   str::Union{String,Nothing} = nothing)
    if str === nothing
        str = string(ex)
    end
    return esc(
        quote
            registerSerializationType($ex, $str)
            encode_type(::Type{<:$ex}) = $str
            # There exist types where equality cannot be discerned from the serialization
            # these types require an id so that equalities can be forced upon load.
            # The ids are only necessary for parent types, checking for element type equality
            # can be done once the parents are known to be equal.
            # For example two serializations of QQ[x] require ids to check for equality.
            # Although they're isomorphic rings, they may want to be treated as seperate
            # This is done since other software might not use symbols in their serialization of QQ[x].
            # Which will then still allow for the distinction between QQ[x] and QQ[y], i.e.
            # whenever there is a possibility (amongst any software system) that the objects
            # cannot be distinguish on a syntatic level we use ids.
            # Types like ZZ, QQ, and ZZ/nZZ do not require ids since there is no syntatic
            # ambiguities in their encodings.

            serialize_with_id(obj::T) where T <: $ex = $uses_id 
            serialize_with_id(T::Type{<:$ex}) = $uses_id 
        end)
end

macro registerSerializationType(ex::Any, str::Union{String,Nothing} = nothing)
    return registerSerializationType(ex, false, str)
end

macro registerSerializationType(ex::Any, uses_id::Bool, str::Union{String,Nothing} = nothing)
    return registerSerializationType(ex, uses_id, str)
end


function encode_type(::Type{T}) where T
    error("unsupported type '$T' for encoding")
end

function decode_type(input::String)
    get(reverseTypeMap, input) do
        error("unsupported type '$input' for decoding")
    end
end

function decode_type(input::Dict{Symbol, Any})
    return decode_type(input[:name])
end

################################################################################
# Encoding helper functions

@doc raw"""
    is_basic_serialization_type(::Type)

During the serialization of types of the form `Vector{T}`, entries
of type `T` will either be serialized as strings if `is_basic_serialization_type`
returns `true`, or serialized as a dict provided the serialization for such a `T`
exists. If `Vector{T}` is serialized with `is_basic_serialization_type(T) = true`
then the `entry_type` keyword is used to store the type `T` as a property of
the vector.

# Examples

```jldoctest
julia> is_basic_serialization_type(ZZRingElem)
true
```
"""
is_basic_serialization_type(::Type) = false
is_basic_serialization_type(::Type{ZZRingElem}) = true
is_basic_serialization_type(::Type{QQFieldElem}) = true
is_basic_serialization_type(::Type{Bool}) = true
is_basic_serialization_type(::Type{String}) = true
is_basic_serialization_type(::Type{Symbol}) = true
# this deals with int32, int64 etc.:precision
is_basic_serialization_type(::Type{T}) where T <: Number = isconcretetype(T)

# ATTENTION
# We need to distinguish between data with a globally defined normal form and data where such a normal form depends on some parameters.
# In particular, this does NOT ONLY depend on the type; see, e.g., FqField.

# Objects having a basic encoding are fundamental objects. They are "atomic"
# with regards to the available types in Oscar (maybe in general), they are base cases when
# a parent tree is constructed via recursion, and they can be decoded from a string or array of
# strings provided their parent is known. All types with basic serialization have a basic encoding
# has_elem_basic_encoding also deals with basic parametrized types.

# a basic serialization type is a type that only requires a single property for serialization.
# examples: QQ, Symbol, Int

# a basic encoded type is a type that only requires a single property or a list of ordered properties
# of the same type and a known parent for serialization.
# example: elements of ZZ/nZZ. Say the parent is ZZ/7ZZ then we can (de)serialize "6", in a sense once the parent is known
# the serialization relies on a "basic" serialization or list of ordered "basic" serializations. We use the word encoding
# here to mean serialization knowing how elements are represented in the parent.

# has_elem_basic_encoding is used for parent types (and only makes sense for parent types)

function has_elem_basic_encoding(obj::T) where T <: Ring
    return is_basic_serialization_type(elem_type(obj))
end

has_basic_encoding(obj::T) where T = is_basic_serialization_type(T)

function has_basic_encoding(obj::T) where T <: RingElem
    return has_elem_basic_encoding(parent(obj))
end

has_elem_basic_encoding(obj::T) where T = false

# used for types that require parents when serialized
type_needs_params(::Type) where Type = false

################################################################################
# High level

function load_ref(s::DeserializerState, id::String)
    if haskey(s.objs, UUID(id))
        loaded_ref = s.objs[UUID(id)]
    else
        ref_dict = s.refs[Symbol(id)]
        ref_dict[:id] = id
        loaded_ref = load_unknown_type(s, ref_dict)
    end
    return loaded_ref
end

function save_as_ref(s::SerializerState, obj::T) where T
    # find ref or create one
    ref = get(s.objmap, obj, nothing)
    if ref !== nothing
        return string(ref)
    end
    
    ref = s.objmap[obj] = uuid4()
    s.refs[Symbol(ref)] = obj

    return string(ref)
end

function save_basic_encoded_type(s::SerializerState, obj::Any, key::Symbol = :data)
    @assert has_basic_encoding(obj)
    active_type = get_active_object_type(s)
    if active_type <: Dict{Symbol, Any}
        add_object(s, string(obj), key)
    else
        add_object(s, string(obj))
    end
end

function save_type_dispatch(s::SerializerState, obj::T, key::Symbol = :data) where T
    # this is used when serializing basic types like "3//4"
    # when s.depth == 0 file should know it belongs to QQ
    if s.depth != 0 && has_basic_encoding(obj)
        save_basic_encoded_type(s, obj, key)
    elseif s.depth == 0 && has_basic_encoding(obj)
        if is_type_serializing_parent(T)
            open_dict(s)
            s.depth += 1
            save_basic_encoded_type(s, obj)
            save_type_dispatch(s, parent(obj), :parent)
            close(s)
            s.depth -= 1
        else
            save_basic_encoded_type(s, obj)
        end
        add_object(s, encode_type(T), :type)
    elseif serialize_with_id(T)
        ref = save_as_ref(s, obj)
        add_object(s, ref, key)
    elseif T <: Vector
        open_dict(s)
        s.depth += 1
        save_internal(s, obj)
        s.depth -= 1
        add_object(s, encode_type(T), :type)
        close(s)
    elseif !Base.issingletontype(T)
        s.depth += 1
        open_dict(s)
        # invoke the actual serializer
        open_dict(s)

        if is_type_serializing_parent(T)
            save_parent_refs(s, parent(obj))
        end
        save_internal(s, obj)
        close(s)
        add_object(s, encode_type(T), :type)
        s.depth -= 1

        if s.depth != 0
            close(s, key)
        end
    else
        # added to catch singleton types
        if s.depth == 0
            # catch edge case
            add_object(s, encode_type(T), :type)
        else
            open_dict(s)
            add_object(s, encode_type(T), :type)
            close(s, key)
        end
    end

    if s.depth == 0
        add_object(s, oscarSerializationVersion, :_ns)

        if !isempty(s.refs)
            add_object(s, s.refs, :refs)
        end
    end
    
    #store_serialized(s, result, key)
end

function save_object(s::SerializerState, x::Vector)
    data_array(s) do
        for elem in x
            save_object(s, elem)
        end
    end
end

# this function is necessary so we can have arrays with entries of different types
function save_object(s::SerializerState, x::T) where T <: Tuple
    data_array(s) do
        for elem in x
            save_object(s, elem)
        end
    end
end

function save_object(s::SerializerState, x::Any, key::Symbol)
    s.key = key
    save_object(s, x)
end

function save_header(s::SerializerState, h::Dict{Symbol, Any}, key::Symbol)
    s.key = key
    data_dict(s) do
        for (k, v) in h
            save_object(s, v, k)
        end
    end
end

# calling this function directly should only happen for the root
function save_typed_object(s::SerializerState, x::T) where T
    if type_needs_params(T)
        save_type_params(s, x, :type)
        save_object(s, x, :data)
    elseif Base.issingletontype(T)
        save_object(s, encode_type(T), :type)
    elseif x isa AbstractVector
        save_object(s, "Vector", :type)
        nested_entry = get_nested_entry(x)
        nested_type = typeof(nested_entry)
        if type_needs_params(nested_type)
            save_type_params(s, nested_entry, :nested_type)
        else
            save_object(s, encode_type(nested_type), :nested_type)
        end
        save_object(s, x, :data)
    else
        save_object(s, encode_type(T), :type)
        save_object(s, x, :data)
    end
end

function save_typed_object(s::SerializerState, x::T, key::Symbol) where T
    s.key = key
    data_dict(s) do 
        save_typed_object(s, x)
    end
end

#function save_type_params(s::SerializerState, x::T, key::Symbol) where T
#    s.key = key
#    data_dict(s) do
#        if serialize_with_id(parent_x)
#            parent_refs = save_parents(s, parent_x)
#            save_object(s, parent_refs, :parents)
#            save_object(s, encode_type(T), :name)
#        else
#            save_typed_object(s, parent_x, :parent)
#            save_object(s, encode_type(T), :name)
#        end
#    end
#end

# ATTENTION
# The load mechanism needs to look at the serialized data first,
# in order to detect objects with a basic encoding.
function load_type_dispatch(s::DeserializerState,
                            ::Type{T}, str::String; parent=nothing) where T
    if parent !== nothing && has_elem_basic_encoding(parent)
        return load_internal_with_parent(s, T, str, parent)
    end
    @assert is_basic_serialization_type(T)
    return load_internal(s, T, str)
end

function load_type_dispatch(s::DeserializerState, ::Type{T}, dict::Dict;
                            parent=nothing) where T
    # File version to be dealt with on first breaking change
    # A file without version number is treated as the "first" version
    if dict[:type] == string(backref_sym)
        backref = s.objs[UUID(dict[:id])]
        backref isa T || error("Backref of incorrect type encountered: $backref !isa $T")
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
    U = decode_type(dict[:type])
    U <: T || U >: T || error("Type in file doesn't match target type: $(dict[:type]) not a subtype of $T")

    Base.issingletontype(T) && return T()

    if parent !== nothing
        result = load_internal_with_parent(s, T, dict[:data], parent)
    elseif type_needs_params(T)
        parent = load_parent_type(s, dict[:type])
        result = load_internal_with_parent(s, T, dict[:data], parent)
    else
        result = load_internal(s, T, dict[:data])
    end

    if haskey(dict, :id)
        s.objs[UUID(dict[:id])] = result
    end
    return result
end

function load_unknown_type(s::DeserializerState, dict::Dict; parent=nothing)
    T = decode_type(dict[:type])
    Base.issingletontype(T) && return T()
    return load_type_dispatch(s, T, dict; parent=parent)
end

# sole purpose of this function is to catch when obj is a ref
function load_unknown_type(s::DeserializerState, str::String)
    return load_ref(s, str)
end

function load_typed_object(s::DeserializerState, dict::Dict{Symbol, Any};
                           parent::Any = nothing)
    T = decode_type(dict[:type])
    # depending on the type, :params is either an object to be loaded or a
    # dict with keys and object values to be loaded
    if type_needs_params(T)
        params = load_type_params(s, T, dict[:type][:params])
        return load_object_with_params(s, T, dict[:data], params)
    elseif Base.issingletontype(T)
        return decode_type(dict[:type])()
    elseif T isa AbstractVector
        nested_type = decode_type(s, dict[:nested_type])
        # not yet sure how recursion will work here
        if type_needs_params(nested_type)
            nested_p_type = load_type_params(s, dict[:nested_type])
            #return load_object()
        else
            #return load_object()
        end
    else
        return load_object(s, T, dict[:data])
    end
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
# Utility functions for parent tree
function load_parent_type(s::DeserializerState, dict::Dict)
    if haskey(dict, :parent)
        return load_unknown_type(s, dict[:parent])
    elseif haskey(dict, :parents)
        return load_parents(s, dict[:parents])[end]
    end
    error("Can't load parent type for $dict[:name]")
end

# loads parent tree
function load_parents(s::DeserializerState, parent_ids::Vector)
    loaded_parents = []

    for id in parent_ids
        loaded_parent = load_ref(s, id)
        push!(loaded_parents, loaded_parent)
    end
    return loaded_parents
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

function get_version_number(v_number::String)
    return VersionNumber(v_number)
end

################################################################################
# Interacting with IO streams and files

"""
    save(io::IO, obj::Any; metadata::MetaData=nothing)
    save(filename::String, obj::Any, metadata::MetaData=nothing)

Save an object `T` to the given io stream
respectively to the file `filename`.

See [`load`](@ref).

# Examples

```jldoctest
julia> meta = metadata(author_orcid="0000-0000-0000-0042", name="the meaning of life the universe and everything")
Oscar.MetaData("0000-0000-0000-0042", "the meaning of life the universe and everything")

julia> save("/tmp/fourtitwo.json", 42; metadata=meta);

julia> load("/tmp/fourtitwo.json")
42
```
"""
function save(io::IO, obj::Any; metadata::Union{MetaData, Nothing}=nothing)
    state = serializer_open(io)
    data_dict(state) do
        save_typed_object(state, obj)

        # this should be handled by serializers in a later commit / PR
        state.key = :refs
        data_dict(state) do 
            for (id, ref_obj) in state.refs
                state.key = id
                data_dict(state) do
                    save_typed_object(state, ref_obj)
                end
            end
        end
        save_header(state, oscarSerializationVersion, :_ns)
    end
    serializer_close(state)
    return nothing
end

function save(filename::String, obj::Any; metadata::Union{MetaData, Nothing}=nothing)
    open(filename, "w") do file
        save(file, obj; metadata=metadata)
    end
    return nothing
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
(Univariate polynomial ring in x over QQ, x)

julia> p = x^2 - x + 1
x^2 - x + 1

julia> save("/tmp/p.json", p)

julia> p_loaded = load("/tmp/p.json", parent=R)
x^2 - x + 1

julia> parent(p_loaded) === R
true

julia> save("/tmp/p_v.json", [p, p])

julia> loaded_p_v = load("/tmp/p_v.json", parent=R)
2-element Vector{QQPolyRingElem}:
 x^2 - x + 1
 x^2 - x + 1

julia> parent(loaded_p_v[1]) === parent(loaded_p_v[2]) === R
true
```
"""
function load(io::IO; parent::Any = nothing, type::Any = nothing)
    state = deserializer_open(io)

    # this should be moved to the serializer at some point
    jsondict = JSON.parse(io, dicttype=Dict{Symbol, Any})

    # handle different namespaces
    @req haskey(jsondict, :_ns) "Namespace is missing"
    _ns = jsondict[:_ns]
    if haskey(_ns, :polymake)
        # If this is a polymake file
        return load_from_polymake(jsondict)
    end
    @req haskey(_ns, :Oscar) "Not an Oscar object"

    # deal with upgrades
    file_version = get_file_version(jsondict)

    if file_version < VERSION_NUMBER
        jsondict = upgrade(jsondict, file_version)
    end

    # add refs to state for referencing during recursion
    if haskey(jsondict, :refs)
        merge!(state.refs, jsondict[:refs])
    end

    # this is a nice to have for later
    # if type !== nothing
    #     return load_type_dispatch(state, type, jsondict; parent=parent)
    # end
    return load_typed_object(state, jsondict; parent=parent)
end

function load(filename::String; parent::Any = nothing, type::Any = nothing)
    open(filename) do file
        return load(file; parent=parent, type=type)
    end
end
