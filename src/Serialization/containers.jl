################################################################################
# Saving and loading vectors

@register_serialization_type Vector uses_params

const MatVecType{T} = Union{Matrix{T}, Vector{T}}

function save_type_params(s::SerializerState, obj::S) where {T, S <:MatVecType{T}}
  save_data_dict(s) do
    save_object(s, encode_type(S), :name)
    if serialize_with_params(T) && !isempty(obj)
      if hasmethod(parent, (T,))
        parents = map(parent, obj)
        parents_all_equal = all(map(x -> isequal(first(parents), x), parents))
        @req parents_all_equal "Not all parents of Vector or Matrix entries are the same, consider using a Tuple"
      end
      save_type_params(s, obj[1], :params)
    else
      save_object(s, encode_type(T), :params)
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:MatVecType})
  T = decode_type(s)
  if serialize_with_params(T) && haskey(s, :params)
    params = load_params_node(s)
    return (T, params)
  end
  return T
end

function load_type_params(s::DeserializerState, ::Type{<:MatVecType}, override_params::Any)
  return (elem_type(override_params), override_params)
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
      loaded_v = params[load_ref(s, x) for x in v]
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
      loaded_v = params[load_ref(s, x) for x in v]
    else
      loaded_v = params[]
      for (i, entry) in enumerate(v)
        push!(loaded_v, load_object(s, params, i))
      end
    end
    return loaded_v
  end
end

# handles nested Vectors
function load_object(s::DeserializerState, ::Type{<: Vector}, params::Tuple)
  T = params[1]
  load_node(s) do v
    if isempty(v)
      return T[]
    else
      loaded_v = []
      for i in 1:length(v)
        load_node(s, i) do _
          push!(loaded_v, load_object(s, T, params[2]))
        end
      end
      return Vector{typeof(loaded_v[1])}(loaded_v)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: Vector}, params::Ring)
  T = elem_type(params)
  loaded_entries = load_array_node(s) do _
    if serialize_with_params(T)
      return load_object(s, T, params)
    else
      return  load_object(s, T)
    end
  end
  return Vector{T}(loaded_entries)
end

################################################################################
# Saving and loading Tuple
@register_serialization_type Tuple uses_params

function save_type_params(s::SerializerState, tup::T) where T <: Tuple
  save_data_dict(s) do
    save_object(s, encode_type(Tuple), :name)
    n = fieldcount(T)
    save_data_array(s, :params) do 
      for i in 1:n
        U = fieldtype(T, i)
        if serialize_with_params(U)
          save_type_params(s, tup[i])
        else
          save_object(s, encode_type(U))
        end
      end
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{Tuple})
  loaded_params = Any[]
  load_array_node(s) do (_, param)
    T = decode_type(s)
    if serialize_with_params(T)
      push!(loaded_params, (T, load_params_node(s)))
    else
      push!(loaded_params, T)
    end
  end
  return loaded_params
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

function load_object(s::DeserializerState, ::Type{<:Tuple}, params::Vector)
  entries = load_array_node(s) do (i, entry)
    if params[i] isa Type
      if serialize_with_id(params[i])
        return load_ref(s)
      else
        return load_object(s, params[i])
      end
    else
      return load_object(s, params[i][1], params[i][2])
    end
  end
  return Tuple(entries)
end

################################################################################
# Saving and loading NamedTuple
@register_serialization_type NamedTuple uses_params

function save_type_params(s::SerializerState, obj::T) where T <: NamedTuple
  save_data_dict(s) do
    save_object(s, encode_type(NamedTuple), :name)
    save_data_dict(s, :params) do
      save_data_array(s, :tuple_params) do 
        for (i, value) in enumerate(values(obj))
          U = fieldtype(T, i)
          if serialize_with_params(U)
            save_type_params(s, value)
          else
            save_object(s, encode_type(U))
          end
        end
      end
      save_object(s, keys(obj), :names)
    end
  end
end

function save_object(s::SerializerState, obj::NamedTuple)
  save_object(s, values(obj))
end

function load_type_params(s::DeserializerState, ::Type{<:NamedTuple})
  loaded_params = Any[]
  load_array_node(s, :tuple_params) do (_, param)
    if param isa String
      push!(loaded_params, decode_type(s))
    else
      T = decode_type(s)
      params = load_params_node(s)
      push!(loaded_params, (T, params))
    end
  end
  load_node(s, :names) do names
    return (names, loaded_params)
  end
end

function load_object(s::DeserializerState, ::Type{<: NamedTuple}, params::Tuple)
  keys, tuple_params = params
  tuple = load_object(s, Tuple, tuple_params)
  keys = Symbol.(keys)
  return NamedTuple{Tuple(keys), typeof(tuple)}(tuple)
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
      return zero_matrix(parent_type(params)(), 0, 0)
    end
    m = reduce(vcat, [
      permutedims(load_object(s, Vector, params, i)) for i in 1:length(entries)
        ])
    return Matrix{params}(m)
  end
end

function load_object(s::DeserializerState, ::Type{<:Matrix}, params::Tuple)
  load_node(s) do entries
    m = reduce(vcat, [
      permutedims(load_object(s, Vector, params, i)) for i in 1:length(entries)
        ])
    return Matrix{params[1]}(m)
  end
end
