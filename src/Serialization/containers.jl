################################################################################
# Saving and loading vectors

@register_serialization_type Vector uses_params

const MatVecType{T} = Union{Matrix{T}, Vector{T}, SRow{T}}

function save_type_params(s::SerializerState, obj::S) where {T, S <:MatVecType{T}}
  save_data_dict(s) do
    save_object(s, encode_type(S), :name)
    if serialize_with_params(T) && !isempty(obj)
      if !(T <: MatVecType) && hasmethod(parent, (T,))
        parents = parent.(obj)
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

# handles nested Vectors
function load_object(s::DeserializerState, ::Type{<: Vector}, params::Tuple)
  T = params[1]
  load_node(s) do v
    if isempty(v)
      return T[]
    else
      loaded_v = []
      len = length(v)
      for i in 1:len
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
# Saving and loading dicts
@register_serialization_type Dict uses_params

function save_type_params(s::SerializerState, obj::Dict{S, Any}) where S <: Union{Symbol, String, Int}
  save_data_dict(s) do
    save_object(s, encode_type(Dict), :name)
    save_data_dict(s, :params) do
      save_object(s, encode_type(S), :key_type)
      for (k, v) in obj
        U = typeof(v)
        if serialize_with_params(U)
          save_type_params(s, v, Symbol(k))
        else
          save_object(s, encode_type(U), Symbol(k))
        end
      end
    end
  end
end

function save_type_params(s::SerializerState, obj::Dict{S, T}) where {S <: Union{Symbol, String, Int}, T}
  save_data_dict(s) do
    save_object(s, encode_type(Dict), :name)
    save_data_dict(s, :params) do
      save_object(s, encode_type(S), :key_type)
      
      if serialize_with_params(T)
        if isempty(obj)
          save_object(s, encode_type(T), :value_type)
        else
          v = first(values(obj))
          save_object(s, encode_type(T), :value_type)
          save_type_params(s, v, :value_params)
        end
      else
        save_object(s, encode_type(T), :value_type)
      end
    end
  end
end

function save_type_params(s::SerializerState, obj::Dict{T, S}) where {T, S}
  save_data_dict(s) do
    save_object(s, encode_type(Dict), :name)
    save_data_dict(s, :params) do
      if serialize_with_params(S)
        if isempty(obj)
          save_object(s, encode_type(S), :value_type)
        else
          v = first(values(obj))
          save_object(s, encode_type(S), :value_type)
          save_type_params(s, v, :value_params)
        end
      else
        save_object(s, encode_type(S), :value_type)
      end
      
      if serialize_with_params(T)
        if isempty(obj)
          save_object(s, encode_type(T), :key_type)
        else
          k = first(keys(obj))
          save_object(s, encode_type(T), :key_type)
          save_type_params(s, k, :key_params)
        end
      else
        save_object(s, encode_type(T), :key_type)
      end
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:Dict})
  if haskey(s, :value_type)
    key_type = load_node(_ -> decode_type(s), s, :key_type)
    value_type = load_node(_ -> decode_type(s), s, :value_type)
    d = Dict{Symbol, Any}(:key_type => key_type, :value_type => value_type)

    if serialize_with_params(value_type)
      d[:value_params] = load_node(s, :value_params) do _
        load_params_node(s)
      end
    end

    if serialize_with_params(key_type)
      d[:key_params] = load_node(s, :key_params) do _
        load_params_node(s)
      end
    end

    return d
  end

  params_dict = Dict{Symbol, Any}()
  for (k, _) in s.obj
    load_node(s, k) do _
      value_type = decode_type(s)

      if serialize_with_params(value_type)
        params_dict[k] = Dict{Symbol, Any}(
          type_key => value_type,
          :params => load_params_node(s)
        )
      else
        params_dict[k] = value_type
      end
    end
  end
  return params_dict
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

function load_object(s::DeserializerState, ::Type{<:Dict}, params::Dict{Symbol, Any})
  key_type = params[:key_type]
  value_type = haskey(params, :value_type) ? params[:value_type] : Any
  dict = Dict{key_type, value_type}()
  
  for (i, (k, _)) in enumerate(s.obj)
    if k == :key_type
      continue
    end
    
    if key_type == Int
      key = parse(Int, string(k))
    elseif haskey(params, :key_params) # type is not Int, String or Symbol
      load_node(s, i) do _
        # 1 is for first entry of tuple which is the key in this case
        key = load_object(s, key_type, params[:key_params], 1) 
      end
    else
      key = key_type(k)
    end

    if value_type != Any
      if serialize_with_params(value_type)
        if haskey(params, :key_params) # key type is not Int, String or Symbol
          load_node(s, i) do _
            # 2 is for second entry of tuple which is the value in this case
            dict[key] = load_object(s, value_type, params[:value_params], 2) 
          end
        else
          dict[key] = load_object(s, value_type, params[:value_params], k)
        end
      else
        if haskey(params, :key_params) # key type is not Int, String or Symbol
          load_node(s, i) do _
            # 2 is for second entry of tuple which is the value in this case
            dict[key] = load_object(s, value_type, 2) 
          end
        else
          dict[key] = load_object(s, value_type, k)
        end
      end
    elseif params[k] isa Type
      dict[key] = load_object(s, params[k], k)
    else
      dict[key] = load_object(s, params[k][type_key], params[k][:params], k)
    end
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
