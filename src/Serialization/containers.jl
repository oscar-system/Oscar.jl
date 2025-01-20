################################################################################
# Saving and loading vectors

@register_serialization_type Vector uses_params

const MatVecType{T} = Union{Matrix{T}, Vector{T}, SRow{T}}

function type_params(obj::S) where {T, S <:MatVecType{T}}
  if isempty(obj)
    return S, nothing
  end

  params = type_params.(obj)
  params_all_equal = all(map(x -> isequal(first(params), x), params))
  @req params_all_equal "Not all params of Vector or Matrix entries are the same, consider using a Tuple for serialization"
  return S, params[1]
end

function save_type_params(s::SerializerState, T::Type{Vector{S}}, ::Nothing) where S
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_object(s, encode_type(S), :params)
  end
end

function save_type_params(s::SerializerState, T::Type, params::Vector)
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_array(s, :params) do
      for param in params
        save_type_params(s, param...)
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type{<: Union{MatVecType, Set}})
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

function load_object(s::DeserializerState, ::Type{<: Vector{params}}) where params
  load_node(s) do v
    if serialize_with_id(params)
      loaded_v::Vector{params} = load_array_node(s) do _
        load_ref(s)
      end
    else
      loaded_v = params[]
      for (i, entry) in enumerate(v)
        push!(loaded_v, load_object(s, params, i))
      end
    end
    return loaded_v
  end
end


function load_object(s::DeserializerState, ::Type{Vector{T}}, params::S) where {T, S}
  load_array_node(s) do _
    if serialize_with_id(T)
      load_ref(s)
    else
      load_object(s, T, params)
    end
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

function load_object(s::DeserializerState, T::Type{<:Matrix{S}}, params::Ring) where S
  load_node(s) do entries
    if isempty(entries)
      return T(undef, 0, 0)
    end

    len = length(entries)
    m = reduce(vcat, [
      permutedims(load_object(s, Vector{S}, params, i)) for i in 1:len
        ])
    return Matrix{elem_type(params)}(m)
  end
end

################################################################################
# Saving and loading Tuple
@register_serialization_type Tuple

function type_params(obj::T) where T <: Tuple
  n = fieldcount(T)
  return T, [
    isnothing(type_params(obj[i])) ? fieldtype(T, i) : type_params(obj[i])
    for i in 1:n]
end

function load_type_params(s::DeserializerState, T::Type{Tuple}) 
  subtype, params = load_node(s, :params) do _
    tuple_params = load_array_node(s) do _
      U = decode_type(s)
      load_type_params(s, U)
    end
    return collect(zip(tuple_params...))
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
  return T, NamedTuple(map(x -> x.first => type_params(x.second), collect(pairs(obj))))
end


# Named Tuples need to preserve order so they are handled seperate from Dict
function save_type_params(s::SerializerState, T::Type{<:NamedTuple}, params::Vector{<:Pair{Symbol, S}}) where S
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_dict(s, :params) do
      save_data_array(s, :names) do
        for param in params
          save_object(s, param.first)
        end
      end
      save_data_array(s, :tuple_params) do
        for param in params
          save_type_params(s, param.second...)
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
function type_params(obj::T) where T <: Dict
  return T, Dict(map(x -> x.first => type_params(x.second), collect(pairs(obj))))
end

function save_type_params(s::SerializerState, T::Type{<:Dict{S, U}}, params::Vector{<:Pair{S, V}}) where {U, V, S <: Union{Int, String}}
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_dict(s, :params) do
      save_object(s, encode_type(S), :key_type)
      isempty(params) && save_object(s, encode_type(U), :value_type)
      for param in params
        save_type_params(s, param.second..., Symbol(param.first))
      end
    end
  end
end

function save_type_params(s::SerializerState, T::Type{<:Dict{Symbol, U}},
                          params::Vector{<:Pair{Symbol, V}}) where {U, V}
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    save_data_dict(s, :params) do
      save_object(s, encode_type(Symbol), :key_type)
      isempty(params) && save_object(s, encode_type(U), :value_type)
      for param in params
        save_type_params(s, param.second..., Symbol(param.first))
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type{Dict}) 
  subtype, params = load_node(s, :params) do obj
    S = load_node(s, :key_type) do _
      decode_type(s)
    end
    params_dict = Dict{S, Any}()
    if S <: Union{String, Symbol, Int}
      value_types = Type[]
      for (k, _) in obj
        k == :key_type && continue
        if k == :value_type
          load_node(s, k) do _
            push!(value_types, decode_type(s))
          end
          continue
        end
        key = S == Int ? parse(Int, string(k)) : S(k)
        params_dict[key] = load_node(s, k) do _
          value_type = decode_type(s)
          return load_type_params(s, value_type)
        end
        push!(value_types, params_dict[key][1])
      end
      return (S, Union{value_types...}), params_dict
    else
      error{"not implemented yet"}
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

function load_object(s::DeserializerState,
                     T::Type{<:Dict{S, U}},
                     params::Dict{S, V}) where {S <: Union{Symbol, String, Int}, U, V}
  dict = T()
  for k in keys(params)
    dict[k] = load_object(s, params[k]..., Symbol(k))
  end
  return dict
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{S, U}}) where {S <: Union{Int, Symbol, String}, U <: Union{Symbol, String, Int}}
  dict = T()
  for k in keys(s.obj)
    dict[S(k)] = load_object(s, U, Symbol(k))
  end
  return dict
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{Int, S}}) where {S <: Union{Symbol, String, Int}}
  dict = T()
  for k in keys(s.obj)
    dict[parse(Int, string(k))] = load_object(s, S, k)
  end
  return dict
end

################################################################################
# Saving and loading sets
@register_serialization_type Set

function type_params(obj::T) where T <: Set
  if isempty(obj)
    return T, nothing
  end

  return T, type_params(first(obj))
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

function load_object(s::DeserializerState, S::Type{<:Set{T}}, params::Ring) where T
  elems = load_array_node(s) do _
    load_object(s, T, params)
  end
  return Set(elems)
end

function load_object(s::DeserializerState, S::Type{<:Set{T}}, ::Nothing) where T
  elems = load_array_node(s) do _
    load_object(s, T, nothing)
  end
  return Set(elems)
end

function load_object(s::DeserializerState, S::Type{<:Set{T}}) where T <: Union{String, Symbol, Int}
  elems = load_array_node(s) do _
    load_object(s, T, nothing)
  end
  return Set(elems)
end

################################################################################
# Sparse rows

@register_serialization_type SRow uses_params

function save_object(s::SerializerState, obj::SRow)
  save_data_array(s) do
    for (i, v) in obj
      save_object(s, (i, v))
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:SRow}, params::Ring)
  pos = Int[]
  entry_type = elem_type(params)
  values = entry_type[]
  load_array_node(s) do _
    push!(pos, load_object(s, Int, 1))
    if serialize_with_params(entry_type)
      push!(values, load_object(s, entry_type, params, 2))
    else
      push!(values, load_object(s, entry_type, 2))
    end
  end
  return sparse_row(params, pos, values)
end
