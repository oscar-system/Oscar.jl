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

# this should eventually become deprecated
function load_object(s::DeserializerState, ::Type{<: Vector}, params::Type)
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

function load_object(s::DeserializerState, ::Type{<: Vector{T}}, R::Ring) where T
  load_array_node(s) do _
    load_object(s, T, R)
  end
end

################################################################################
# Saving and loading matrices
@register_serialization_type Matrix uses_params

function save_object(s::SerializerState, mat::Matrix)
  m, n = size(mat)
  save_data_array(s) do
    for i in 1:m
      save_object(s, [mat[i, j] for j in 1:n])
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:Matrix}, params::Type)
  load_node(s) do entries
    if isempty(entries)
      return Matrix{params}(undef, 0, 0)
    end
    len = length(entries)
    m = reduce(vcat, [
      permutedims(load_object(s, Vector, params, i)) for i in 1:len
        ])
    return Matrix{params}(m)
  end
end

function load_object(s::DeserializerState, ::Type{<:Matrix}, params::Tuple)
  load_node(s) do entries
    if isempty(entries)
      return Matrix{params[1]}(undef, 0, 0)
    end

    len = length(entries)
    m = reduce(vcat, [
      permutedims(load_object(s, Vector, params, i)) for i in 1:len
        ])
    return Matrix{params[1]}(m)
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

################################################################################
# Saving and loading NamedTuple
@register_serialization_type NamedTuple

function type_params(obj::T) where T <: NamedTuple
  return T, NamedTuple(map(x -> x.first => type_params(x.second), collect(pairs(obj))))
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

function save_object(s::SerializerState, obj::Dict{S, T}) where {S <: Union{Symbol, String, Int}, T}
  save_data_dict(s) do
    for (k, v) in obj
      if !Base.issingletontype(typeof(v))
        save_object(s, v, Symbol(k))
      end
    end
  end
end

function save_object(s::SerializerState, obj::Dict) 
  save_data_array(s) do
    for (k, v) in obj
      save_object(s, (k, v))
    end
  end
end

function load_object(s::DeserializerState, ::Type{Dict{String, Int}})
  return Dict{String, Int}(string(k) => parse(Int, v) for (k,v) in s.obj)
end

function load_object(s::DeserializerState, ::Type{Dict{Int, Int}})
  return Dict{Int, Int}(parse(Int,string(k)) => parse(Int, v) for (k,v) in s.obj)
end

function load_object(s::DeserializerState,
                     T::Type{<:Dict{S, U}},
                     params::Dict{S, V}) where {S <: BasicKeyTypes, U, V}
  dict = T()
  for k in keys(params)
    dict[k] = load_object(s, params[k]..., Symbol(k))
  end
  return dict
end

################################################################################
# Saving and loading sets
@register_serialization_type Set uses_params

function save_type_params(s::SerializerState, obj::Set{T}) where T
  save_data_dict(s) do
    save_object(s, encode_type(Set), :name)
    if serialize_with_params(T) && !isempty(obj)
        save_type_params(s, first(obj), :params)
      else
        save_object(s, encode_type(T), :params)
      end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:Set})
  T = decode_type(s)
  if serialize_with_params(T) && haskey(s, :params)
    params = load_params_node(s)
    return (T, params)
  end
  return T
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

function load_object(s::DeserializerState, ::Type{<: Set}, params::Any)
  load_node(s) do v
    if serialize_with_id(params)
      loaded_v = params[load_ref(s, x) for x in v]
    else
      loaded_v = params[]
      for (i, entry) in enumerate(v)
        push!(loaded_v, load_object(s, params, i))
      end
    end
    return Set(loaded_v)
  end
end

# handles nested
function load_object(s::DeserializerState, ::Type{<: Set}, params::Tuple)
  T = params[1]
  load_node(s) do v
    if isempty(v)
      return Set{T}()
    else
      loaded_v = Set{T}()
      len = length(v)
      for i in 1:len
        load_node(s, i) do _
          push!(loaded_v, load_object(s, T, params[2]))
        end
      end
      return Set{typeof(first(loaded_v))}(loaded_v)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: Set}, params::Ring)
  T = elem_type(params)
  loaded_entries = load_array_node(s) do _
    if serialize_with_params(T)
      return load_object(s, T, params)
    else
      return load_object(s, T)
    end
  end
  return Set{T}(loaded_entries)
end

################################################################################
# Sparse rows

@register_serialization_type SRow uses_params

function save_object(s::SerializerState, obj::SRow)
  save_data_array(s) do
    for (i, v) in collect(obj)
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
