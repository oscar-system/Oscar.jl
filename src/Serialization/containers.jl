const MatVecType{T} = Union{Matrix{T}, Vector{T}, SRow{T}}
const ContainerTypes = Union{MatVecType, Set, Dict, Tuple, NamedTuple}

function params_all_equal(params::MatVecType{<:TypeParams})
  all(map(x -> isequal(first(params), x), params))
end

function type_params(obj::S) where {T, S <:MatVecType{T}}
  if isempty(obj)
    return TypeParams(S, TypeParams(T, nothing))
  end
  
  params = type_params.(obj)
  @req params_all_equal(params) "Not all params of Vector or Matrix entries are the same, consider using a Tuple for serialization"
  return TypeParams(S, params[1])
end

function has_empty_entries(obj::T) where T
  return false
end

function has_empty_entries(obj::T) where T <: ContainerTypes
  isempty(obj) && return true
  any(has_empty_entries, obj)  && return true
  return false
end

function type_params(obj::S) where {T <: ContainerTypes, S <:MatVecType{T}}
  isempty(obj) && return TypeParams(S, TypeParams(T, nothing))

  # empty entries can inherit params from the rest of the collection
  params = type_params.(filter(!has_empty_entries, obj))
  @req params_all_equal(params) "Not all params of Vector or Matrix entries are the same, consider using a Tuple for serialization"
  return TypeParams(S, params[1])
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

function load_type_params(s::DeserializerState, T::Type{<: MatVecType})
  !haskey(s, :params) && return T, nothing
  subtype, params = load_node(s, :params) do _
    U = decode_type(s)
    subtype, params = load_type_params(s, U)
  end
  return T{subtype}, params
end

function save_object(s::SerializerState, x::Vector)
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
  load_node(s) do v
    if serialize_with_id(params)
      loaded_v::Vector{params} = load_array_node(s) do _
        load_ref(s)
      end
    else
      isempty(s.obj) && return params[]
      loaded_v = load_array_node(s) do obj
        load_object(s, params)
      end
    end
    return loaded_v
  end
end


function load_object(s::DeserializerState, ::Type{Vector{T}}, params::S) where {T, S}
  isempty(s.obj) && return T[]
  v = load_array_node(s) do _
    if serialize_with_id(T)
      load_ref(s)
    else
      load_object(s, T, params)
    end
  end
  return v
end

function load_object(s::DeserializerState, T::Type{Vector{U}}, ::Nothing) where U
  isempty(s.obj) && return U[]
  return load_array_node(s) do _
    load_object(s, U)
  end
end

################################################################################
# Saving and loading matrices
@register_serialization_type Matrix

function save_object(s::SerializerState, mat::Matrix)
  m, n = size(mat)
  save_data_array(s) do
    for i in 1:m
      save_object(s, [mat[i, j] for j in 1:n])
    end
  end
end

function load_object(s::DeserializerState, T::Type{<:Matrix{S}}) where S
  load_node(s) do entries
    if isempty(entries)
      return T(undef, 0, 0)
    end
    len = length(entries)
    m = reduce(vcat, [
      permutedims(load_object(s, Vector{S}, i)) for i in 1:len
        ])
    return T(m)
  end
end

function load_object(s::DeserializerState, T::Type{<:Matrix{S}}, params::NCRing) where S
  load_node(s) do entries
    if isempty(entries)
      return T(undef, 0, 0)
    end
    len = length(entries)
    m = reduce(vcat, [
      permutedims(load_object(s, Vector{S}, params, i)) for i in 1:len
        ])
    return T(m)
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
  subtype, params = load_node(s, :params) do _
    tuple_params = load_array_node(s) do _
      U = decode_type(s)
      load_type_params(s, U)
    end
    return Tuple([x[1] for x in tuple_params]), Tuple(x[2] for x in tuple_params)
  end
  return T{subtype...}, params
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

function load_object(s::DeserializerState, T::Type{<:Tuple}, params::Tuple)
  entries = load_array_node(s) do (i, entry)
    S = fieldtype(T, i)
    if serialize_with_id(S)
      return load_ref(s)
    else
      return load_object(s, S, params[i])
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
    NamedTuple(map(x -> x.first => type_params(x.second), collect(pairs(obj))))
  )
end

# Named Tuples need to preserve order so they are handled seperate from Dict
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
  subtype, params = load_node(s, :params) do obj
    tuple_params = load_array_node(s, :tuple_params) do _
      U = decode_type(s)
      load_type_params(s, U)
    end
    tuple_types, named_tuple_params = collect(zip(tuple_params...))
    names = load_object(s, Vector{Symbol}, :names)
    return (tuple(names...), Tuple{tuple_types...}), named_tuple_params
  end
  return T{subtype...}, params
end

function save_object(s::SerializerState, obj::NamedTuple)
  save_object(s, values(obj))
end

function load_object(s::DeserializerState, T::Type{<: NamedTuple}, params::Tuple)
  return T(load_object(s, Tuple{fieldtypes(T)...}, params))
end

################################################################################
# Saving and loading dicts
@register_serialization_type Dict

function type_params(obj::T) where {S <: Union{Symbol, Int, String},
                                    T <: Dict{S, Any}}
  return TypeParams(
    T,
    :key_params => TypeParams(S, nothing),
    map(x -> x.first => type_params(x.second), collect(pairs(obj)))...
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
  @req params_all_equal(key_params) "Not all params of keys in $obj are the same"

  value_params = type_params.(collect(values(obj)))
  @req params_all_equal(value_params) "Not all params of values in $obj are the same"

  return TypeParams(
    T, 
    :key_params => first(key_params),
    :value_params => first(value_params)
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
  subtype, params = load_node(s, :params) do obj
    if haskey(s, :value_params)
      S, key_params = load_node(s, :key_params) do params
        params isa String && return decode_type(s), nothing
        load_type_params(s, decode_type(s))
      end

      U, value_params = load_node(s, :value_params) do _
        load_type_params(s, decode_type(s))
      end

      isnothing(key_params) && return (S, U), value_params
      isnothing(key_params) && isnothing(value_params) && return (S, U), nothing
      return (S, U), Dict(:key_params => key_params, :value_params => value_params)
    else
      S, key_params = load_node(s, :key_params) do _
        decode_type(s), nothing
      end
      params_dict = Dict{S, Any}()
      value_types = Type[]
      for (k, _) in obj
        k == :key_params && continue
        key = S <: Integer ? parse(Int, string(k)) : S(k)
        params_dict[key] = load_node(s, k) do _
          value_type = decode_type(s)
          return load_type_params(s, value_type)
        end
        push!(value_types, params_dict[key][1])
      end
      params_dict = isempty(params_dict) ? nothing : params_dict
      return (S, Union{value_types...}), params_dict
    end
  end
  return Dict{subtype...}, params
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
                     T::Type{<:Dict{S, Any}},
                     params::Dict{S, Any}) where {S <: Union{Symbol, String, Int}}
  dict = T()
  for k in keys(params)
    # has no data, hence no key was generated on the data side
    if Base.issingletontype(params[k][1])
      dict[k] = params[k][1]()
    else
      dict[k] = load_object(s, params[k]..., Symbol(k))
    end
  end
  return dict
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{S, U}}) where {S <: Union{Symbol, String}, U}
  dict = T()
  for k in keys(s.obj)
    dict[S(k)] = load_object(s, U, Symbol(k))
  end
  return dict
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{S, U}}, ::Nothing) where {S <: Union{Symbol, String}, U}
  dict = T()
  for k in keys(s.obj)
    dict[S(k)] = load_object(s, U, Symbol(k))
  end
  return dict
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{Int, Int}})
  dict = T()
  for k in keys(s.obj)
    dict[parse(Int, string(k))] = load_object(s, Int, k)
  end
  return dict
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{Int, S}}, params::Nothing) where S
  dict = T()
  for k in keys(s.obj)
    dict[parse(Int, string(k))] = load_object(s, S, k)
  end
  return dict
end

# here to handle ambiguities
function load_object(s::DeserializerState, T::Type{Dict{Int, Int}}, key::Union{Symbol, Int})
  load_node(s, key) do _
    load_object(s, T)
  end
end

# here to handle ambiguities
function load_object(s::DeserializerState, T::Type{Dict{String, Int}}, key::Union{Symbol, Int})
  load_node(s, key) do _
    load_object(s, T)
  end
end

# here to handle ambiguities
function load_object(s::DeserializerState, T::Type{Dict{Symbol, Int}}, key::Union{Symbol, Int})
  load_node(s, key) do _
    load_object(s, T)
  end
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{S, U}},
                     params::Any) where {S, U}
  if params isa Dict
    if haskey(params, :value_params)
      pairs = load_array_node(s) do _
        load_object(s, Tuple{S, U}, (params[:key_params], params[:value_params]))
      end
      return T(k => v for (k, v) in pairs)
    else
      dict = Dict{S, Any}()
      value_types = Type[]
      for k in keys(s.obj)
        key = S <: Integer ? parse(S, string(k)) : S(k)
        value_type, param = params[key]
        v = load_object(s, value_type, param, k)
        dict[key] = load_object(s, value_type, param, k)
        push!(value_types, typeof(v))
      end
      isempty(value_types) && return T()
      value_params = type_params.(collect(values(dict)))
      value_type = params_all_equal(value_params) ? typejoin(unique(value_types)...) : Any
      return Dict{S, value_type}(dict)
    end
  else
    dict = Dict{S, Any}()
    value_types = Type[]
    for k in keys(s.obj)
      v = load_object(s, U, params, k)
      key = S <: Integer ? parse(S, string(k)) : S(k)
      dict[key] = v
      push!(value_types, typeof(v))
    end
    isempty(value_types) && return T()
    value_params = type_params.(collect(values(dict)))
    value_type = params_all_equal(value_params) ? typejoin(unique(value_types)...) : Any
    return Dict{S, value_type}(dict)
  end
end

################################################################################
# Saving and loading sets
@register_serialization_type Set

function type_params(obj::T) where {S, T <: Set{S}}
  if isempty(obj)
    return TypeParams(T, TypeParams(S, nothing))
  end
  return TypeParams(T, type_params(first(obj)))
end

function load_type_params(s::DeserializerState, T::Type{<: Set})
  !haskey(s, :params) && return T, nothing
  subtype, params = load_node(s, :params) do _
    U = decode_type(s)
    subtype, params = load_type_params(s, U)
  end
  return T{subtype}, params
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

function load_object(s::DeserializerState, S::Type{<:Set{T}}, params::Any) where T
  elems = load_array_node(s) do _
    load_object(s, T, params)
  end
  return S(elems)
end

function load_object(s::DeserializerState, S::Type{<:Set{T}}, ::Nothing) where T
  elems = load_array_node(s) do _
    load_object(s, T, nothing)
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

function load_object(s::DeserializerState, ::Type{<:SRow}, params::NCRing)
  pos = Int[]
  entry_type = elem_type(params)
  values = entry_type[]
  load_array_node(s) do _
    push!(pos, load_object(s, Int, 1))
    push!(values, load_object(s, entry_type, params, 2))
  end
  return sparse_row(params, pos, values)
end
