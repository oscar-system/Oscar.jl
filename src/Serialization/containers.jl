const MatVecType{T} = Union{Matrix{T}, Vector{T}, SRow{T}} where T
const ContainerTypes = Union{MatVecType, Set, Dict, Tuple, NamedTuple, Array}

# handle Vector and Matrix instances, i.e. Array instances with 1 or 2 dimensions
function type_params(obj::S) where {T, S <: MatVecType{T}}
  isempty(obj) && return TypeParams(S, TypeParams(T, nothing))
  
  params = type_params.(obj)
  @req allequal(params) "Not all params of the entries are the same, consider using a Tuple for serialization"
  return TypeParams(S, params[1])
end

# handle Array instances with >= 3 dimensions
function type_params(obj::S) where {N, T, S <: Array{T, N}}
  isempty(obj) && return TypeParams(S, :subtype_params => TypeParams(T, nothing), :dims => N)
  
  params = type_params.(obj)
  @req allequal(params) "Not all params of the entries are the same, consider using a Tuple for serialization"
  return TypeParams(S, :subtype_params => params[1], :dims => N)
end

has_empty_entries(obj) = false

has_empty_entries(obj::ContainerTypes) = isempty(obj) || any(has_empty_entries, obj)

# handle Vector and Matrix instances containing other containers
function type_params(obj::S) where {T <: ContainerTypes, S <: MatVecType{T}}
  isempty(obj) && return TypeParams(S, TypeParams(T, nothing))

  # empty entries can inherit params from the rest of the collection
  non_empty_entries = filter(!has_empty_entries, obj)
  params = type_params.(non_empty_entries)
  @req allequal(params) "Not all params of the entries are the same, consider using a Tuple for serialization"

  # need to check if any of the inner containers are non empty, and if there is at least one
  # then we can use that one to get the type for the entire nested container
  isempty(params) && return TypeParams(S, type_params(first(obj)))
  return TypeParams(S, params[1])
end

function save_type_params(s::SerializerState,
                          tp::TypeParams{T, <:Tuple{Vararg{Pair}}}) where {S, N,  T <: Array{S, N}}
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_dict(s, :params) do
      for (k, v) in params(tp)
        if k == :dims
          save_object(s, v, k)
        else
          save_type_params(s, v, k)
        end
      end
    end
  end
end

# this function is forced to deal with all array types
function load_type_params(s::DeserializerState, T::Type{<: Union{Array, MatVecType}})
  !haskey(s, :params) && return TypeParams(T, nothing)
  return load_node(s, :params) do _
    if haskey(s, :dims)
      U = load_node(s, :subtype_params) do _
        return decode_type(s)
      end

      dims = load_object(s, Int, :dims)
      sub_tp = load_type_params(s, U, :subtype_params)
      return TypeParams(T{type(sub_tp), dims}, Oscar.params(sub_tp))
    else
      sub_tp = load_type_params(s, decode_type(s))
      return TypeParams(T{type(sub_tp)}, Oscar.params(sub_tp))
    end
  end
end

################################################################################
# loads to handle

function load_object(s::DeserializerState, T::Type{Vector{S}},
                     key::Union{Symbol, Int}) where S
  load_node(s, key) do _
    load_object(s, T)
  end
end

################################################################################
# Saving and loading vectors
@register_serialization_type Vector

# see Array above for load_type_params for general Arrays which include Matrix
# and Vector

function save_object(s::SerializerState, x::AbstractVector{S}) where S
  save_data_array(s) do
    for elem in x
      if serialize_with_id(typeof(elem))
        ref = save_as_ref(s, elem)
        save_object(s, ref)
      else
        save_object(s, elem)
      end
    end
  end
end

function load_object(s::DeserializerState, T::Type{<: Vector{params}}) where params
  load_node(s) do _
    if serialize_with_id(params)
      loaded_v = load_array_node(s) do _
        load_ref(s)
      end::Vector{params}
    else
      isempty(s) && return params[]
      loaded_v = load_array_node(s) do _
        load_object(s, params)
      end
    end
    return loaded_v
  end
end


function load_object(s::DeserializerState, tp::TypeParams{Vector{T}}) where T
  p = Oscar.params(tp)
  isempty(s) && return T[]
  v = load_array_node(s) do _
    if serialize_with_id(T)
      load_ref(s)
    else
      load_object(s, TypeParams(T, p))
    end
  end
  return v
end

################################################################################
# Saving and loading matrices
@register_serialization_type Matrix

function save_object(s::SerializerState, mat::AbstractMatrix{S}) where {S}
  save_data_array(s) do
    for r in eachrow(mat)
      save_object(s, r)
    end
  end
end

# this needs to be Matrix to avoid ambiguities, this is also ok
# since types on load will be Matrix
function load_object(s::DeserializerState, T::Type{<:Matrix{S}}) where {S}
  load_node(s) do _
    if isempty(s)
      return T(undef, 0, 0)
    end
    len = length(s.obj)
    m = stack([
      load_object(s, Vector{S}, i) for i in 1:len
        ]; dims=1)
    return T(m)
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:Matrix{S}}) where S
  T = type(tp)
  p = Oscar.params(tp)
  load_node(s) do _
    if isempty(s)
      return T(undef, 0, 0)
    end
    len = length(s.obj)
    m = stack([
      load_object(s, TypeParams(Vector{S}, p), i) for i in 1:len
        ]; dims=1)
    return T(m)
  end
end

function load_object(s::DeserializerState, T::Type{<:Matrix{S}}, key::Union{Int, Symbol}) where {S}
  return load_node(s, key) do _
    load_object(s, T)
  end
end

################################################################################
# Multidimensional Arrays
@register_serialization_type Array "MultiDimArray"

function save_object(s::SerializerState, arr::AbstractArray{T, N}) where {T, N}
  save_data_array(s) do
    for slice in eachslice(arr; dims=1)
      save_object(s, slice)
    end
  end
end

function load_object(s::DeserializerState, T::Type{Array{S, N}}) where {S, N}
  load_node(s) do _
    if isempty(s)
      return T(undef, [0 for _ in 1:N]...)
    end
    len = length(s.obj)
    m = [load_object(s, Array{S, N - 1}, i) for i in 1:len]
    return stack(m; dims=1)
  end
end

function load_object(s::DeserializerState, tp::TypeParams{Array{S, N}}) where {S, N}
  T = type(tp)
  p = Oscar.params(tp)
  load_node(s) do _
    if isempty(s)
      return T(undef, [0 for _ in 1:N]...)
    end
    len = length(s.obj)
    m = [load_object(s, TypeParams(Array{S, N - 1}, p), i) for i in 1:len]
    return stack(m; dims=1)
  end
end

function load_object(s::DeserializerState, T::Type{Array{S, N}}, key::Union{Int, Symbol}) where {S, N}
  return load_node(s, key) do _
    load_object(s, T)
  end
end

################################################################################
# Saving and loading Tuple
@register_serialization_type Tuple

function type_params(obj::T) where T <: Tuple
  return TypeParams(T, type_params.(obj))
end

function save_type_params(s::SerializerState, tp::TypeParams{T}) where T <: Tuple
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_array(s, :params) do
      for (i, param_tp) in enumerate(params(tp))
        save_type_params(s, param_tp)
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type{Tuple})
  !haskey(s, :params) && return TypeParams(T{}, nothing)
  return load_node(s, :params) do _
    tuple_params = load_array_node(s) do _
      load_type_params(s, decode_type(s))
    end
    tuple_types = Tuple([type(x) for x in tuple_params])
    raw_params = Tuple(Oscar.params(x) for x in tuple_params)
    return TypeParams(T{tuple_types...}, raw_params)
  end
end

function save_object(s::SerializerState, obj::Tuple)
  save_data_array(s) do
    for entry in obj
      if serialize_with_id(typeof(entry))
        ref = save_as_ref(s, entry)
        save_object(s, ref)
      else
        save_object(s, entry)
      end
    end
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:Tuple, <:Tuple})
  T = type(tp)
  p = Oscar.params(tp)
  entries = load_array_node(s) do (i, entry)
    S = fieldtype(T, i)
    if serialize_with_id(S)
      return load_ref(s)
    else
      return load_object(s, TypeParams(S, p[i]))
    end
  end
  return Tuple(entries)
end

function load_object(s::DeserializerState, T::Type{<:Tuple})
  entries = load_array_node(s) do (i, entry)
    S = fieldtype(T, i)
    load_object(s, S)
  end
  return T(entries)
end

################################################################################
# Saving and loading NamedTuple
@register_serialization_type NamedTuple

function type_params(obj::T) where T <: NamedTuple
  return TypeParams(
    T, 
    NamedTuple(x.first => type_params(x.second) for x in pairs(obj))
  )
end

# Named Tuples need to preserve order so they are handled separate from Dict
function save_type_params(s::SerializerState, tp::TypeParams{<:NamedTuple})
  T = type(tp)
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_dict(s, :params) do
      save_data_array(s, :names) do
        for name in keys(params(tp))
          save_object(s, name)
        end
      end
      save_data_array(s, :tuple_params) do
        for (i, param_tp) in enumerate(values(params(tp)))
          save_type_params(s, param_tp)
        end
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type{NamedTuple})
  return load_node(s, :params) do _
    tuple_params = load_array_node(s, :tuple_params) do _
      load_type_params(s, decode_type(s))
    end
    tuple_types = Tuple([type(x) for x in tuple_params])
    raw_params = Tuple(Oscar.params(x) for x in tuple_params)
    names = load_object(s, Vector{Symbol}, :names)
    return TypeParams(T{tuple(names...), Tuple{tuple_types...}}, raw_params)
  end
end

function save_object(s::SerializerState, obj::NamedTuple)
  save_object(s, values(obj))
end

function load_object(s::DeserializerState, tp::TypeParams{<:NamedTuple, <:Tuple})
  T = type(tp)
  p = Oscar.params(tp)
  return T(load_object(s, TypeParams(Tuple{Base.fieldtypes(T)...}, p)))
end

################################################################################
# Saving and loading dicts
@register_serialization_type Dict

function type_params(obj::T) where {S <: Union{Symbol, Int, String},
                                    T <: Dict{S, Any}}
  return TypeParams(
    T,
    :key_params => TypeParams(S, nothing),
    (x.first => type_params(x.second) for x in pairs(obj))...
  )
end

function type_params(obj::T) where {U, S, T <: Dict{S, U}}
  if isempty(obj)
    return TypeParams(
      T,
      :key_params => TypeParams(S, nothing),
      :value_params => TypeParams(U, nothing)
    )
  end
  key_params = Oscar.type_params.(collect(keys(obj)))
  @req allequal(key_params) "Not all type parameters of keys in $obj are the same"

  value_params = type_params.(collect(values(obj)))
  @req allequal(value_params) "Not all type parameters of values in $obj are the same"

  return TypeParams(
    T, 
    :key_params => first(key_params),
    :value_params => first(value_params),
  )
end

function save_type_params(
  s::SerializerState,
  tp::TypeParams{Dict{S, T}, <:Tuple{Vararg{Pair}}}) where {T, S <: Union{Symbol, Int, String}}
  save_data_dict(s) do
    save_object(s, encode_type(Dict), :name)
    save_data_dict(s, :params) do
      isempty(params(tp)) && save_object(s, encode_type(T), :value_params)
      for (k, param_tp) in params(tp)
        save_type_params(s, param_tp, Symbol(k))
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type{Dict})
  return load_node(s, :params) do _
    if haskey(s, :value_params)
      key_tp = load_node(s, :key_params) do _
        if is_string(s)
          TypeParams(decode_type(s), nothing)
        else
          load_type_params(s, decode_type(s))
        end
      end
      S = type(key_tp)
      key_p = Oscar.params(key_tp)

      value_tp = load_node(s, :value_params) do _
        load_type_params(s, decode_type(s))
      end
      U = type(value_tp)
      value_p = Oscar.params(value_tp)

      return TypeParams(Dict{S, U}, Dict{Symbol, Any}(:key_params => key_p, :value_params => value_p))
    else
      S = load_node(s, :key_params) do _
        decode_type(s)
      end
      params_dict = Dict{S, Any}()
      value_types = Type[]
      for k in propertynames(s.obj)
        k == :key_params && continue
        key = S <: Integer ? parse(Int, string(k)) : S(k)
        params_dict[key] = load_node(s, Symbol(k)) do _
          return load_type_params(s, decode_type(s))
        end
        push!(value_types, type(params_dict[key]))
      end
      params_val = isempty(params_dict) ? nothing : params_dict
      return TypeParams(Dict{S, Union{value_types...}}, params_val)
    end
  end
end

function save_object(s::SerializerState, obj::Dict{S, T}) where {S <: Union{Symbol, String, Int}, T}
  save_data_dict(s) do
    for (k, v) in obj
      if !Base.issingletontype(typeof(v))
        save_object(s, v, Symbol(k))
      end
    end
  end
end

function save_object(s::SerializerState, obj::Dict{S, T}) where {S, T}
  save_data_array(s) do
    for (k, v) in obj
      save_object(s, (k, v))
    end
  end
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Dict{S, Any}, <:Dict}) where {S <: Union{Symbol, String, Int}}
  T = type(tp)
  p = Oscar.params(tp)
  dict = T()
  for k in keys(p)
    sk = S <: Symbol ? Symbol(k) : (S <: Int ? parse(Int, string(k)) : String(k))
    # has no data, hence no key was generated on the data side
    if Base.issingletontype(type(p[k]))
      dict[sk] = type(p[k])()
    else
      dict[sk] = load_object(s, p[k], Symbol(k))
    end
  end
  return dict
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Dict{S, U}, Nothing}) where {S <: Union{Symbol, String, Int}, U}
  T = type(tp)
  dict = Dict{S, U}()
  for k in propertynames(s.obj)
    key = S <: Integer ? parse(S, string(k)) : S(k)
    dict[key] = load_object(s, U, k)
  end
  return dict
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{S, U}}) where {S <: Union{Symbol, String, Int}, U}
  load_object(s, TypeParams(T, nothing))
end

function load_object(s::DeserializerState, tp::TypeParams{<:Dict{S, U}, Nothing}) where {S, U}
  T = type(tp)
  pairs = load_array_node(s) do _
    load_object(s, TypeParams(Tuple{S, U}, (nothing, nothing)))
  end::Vector{Tuple{S, U}}
  return T(k => v for (k, v) in pairs)
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Dict{S, U}, <:Dict}) where {S <: Union{Int, String, Symbol}, U}
  T = type(tp)
  p = Oscar.params(tp)
  if haskey(p, :value_params)
    pairs = Tuple{S, U}[]
    for k in propertynames(s.obj)
      key = S <: Integer ? parse(S, string(k)) : S(k)
      push!(pairs, (key, load_object(s, TypeParams(U, p[:value_params]), k)))
    end
    isempty(pairs) && return T()
    _, v = first(pairs)
    return Dict{S, typeof(v)}(k => v for (k, v) in pairs)
  else
    dict = Dict{S, Any}()
    value_types = Type[]
    for k in propertynames(s.obj)
      key = S <: Integer ? parse(S, string(k)) : S(k)
      param_key = haskey(p, key) ? key : string(k)
      v = load_object(s, p[param_key], k)
      dict[key] = v
      push!(value_types, typeof(v))
    end
    isempty(value_types) && return T()
    value_params_list = type_params.(collect(values(dict)))
    value_type = allequal(value_params_list) ? typejoin(unique(value_types)...) : Any
    return Dict{S, value_type}(dict)
  end
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Dict{S, U}, <:Dict}) where {S, U}
  T = type(tp)
  p = Oscar.params(tp)
  pairs = load_array_node(s) do _
    load_object(s, TypeParams(Tuple{S, U}, (p[:key_params], p[:value_params])))
  end
  isempty(pairs) && return T()
  k1, v1 = first(pairs)
  return Dict{typeof(k1), typeof(v1)}(k => v for (k, v) in pairs)
end

################################################################################
# Saving and loading sets
@register_serialization_type Set

function type_params(obj::T) where {S, T <: Set{S}}
  isempty(obj) && return TypeParams(T, TypeParams(S, nothing))
  return TypeParams(T, type_params(first(obj)))
end

function load_type_params(s::DeserializerState, T::Type{<: Set})
  !haskey(s, :params) && return TypeParams(T, nothing)
  return load_node(s, :params) do _
    sub_tp = load_type_params(s, decode_type(s))
    return TypeParams(T{type(sub_tp)}, Oscar.params(sub_tp))
  end
end

function save_object(s::SerializerState, x::Set)
  save_data_array(s) do
    for elem in x
      if serialize_with_id(typeof(elem))
        ref = save_as_ref(s, elem)
        save_object(s, ref)
      else
        save_object(s, elem)
      end
    end
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:Set{T}}) where T
  S = type(tp)
  p = Oscar.params(tp)
  elems = load_array_node(s) do _
    load_object(s, TypeParams(T, p))
  end
  return S(elems)
end

function load_object(s::DeserializerState, S::Type{<:Set{T}}) where T
  elems = load_array_node(s) do _
    load_object(s, T)
  end
  return S(elems)
end

################################################################################
# Sparse rows

@register_serialization_type SRow 

function save_object(s::SerializerState, obj::SRow)
  save_data_array(s) do
    for (i, v) in obj
      save_object(s, (i, v))
    end
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:SRow, <:NCRing})
  ring = Oscar.params(tp)
  pos = Int[]
  entry_type = elem_type(ring)
  values = entry_type[]
  load_array_node(s) do _
    push!(pos, load_object(s, Int, 1))
    push!(values, load_object(s, TypeParams(entry_type, ring), 2))
  end
  return sparse_row(ring, pos, values)
end
