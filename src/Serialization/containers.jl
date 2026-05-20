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
      for (k, v) in parameters(tp)
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
  return load_node(s, :params) do
    if haskey(s, :dims)
      dims = load_object(s, Int, :dims)
      sub_tp = load_node(s, :subtype_params) do
        load_type_params(s, decode_type(s))
      end
      return TypeParams(T{type(sub_tp), dims}, parameters(sub_tp))
    else
      sub_tp = load_type_params(s, decode_type(s))
      return TypeParams(T{type(sub_tp)}, parameters(sub_tp))
    end
  end
end

################################################################################
# loads to handle

function load_object(s::DeserializerState, T::Type{Vector{S}},
                     key::Union{Symbol, Int}) where S
  load_node(s, key) do
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
  load_node(s) do
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
  p = parameters(tp)
  elem_tp = p isa TypeParams ? p : TypeParams(T, p)
  isempty(s) && return T[]
  v = load_array_node(s) do _
    if serialize_with_id(T)
      load_ref(s)
    else
      load_object(s, elem_tp)
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
  load_node(s) do
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
  p = parameters(tp)
  elem_tp = p isa TypeParams ? p : TypeParams(S, p)
  load_node(s) do
    if isempty(s)
      return T(undef, 0, 0)
    end
    len = length(s.obj)
    m = stack([
      load_object(s, TypeParams(Vector{S}, elem_tp), i) for i in 1:len
        ]; dims=1)
    return T(m)
  end
end

function load_object(s::DeserializerState, T::Type{<:Matrix{S}}, key::Union{Int, Symbol}) where {S}
  return load_node(s, key) do
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
  load_node(s) do
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
  p = parameters(tp)
  sub_tp = if p isa Tuple{Vararg{Pair}}
    tp[:subtype_params]
  elseif p isa TypeParams
    p
  else
    TypeParams(S, p)
  end
  load_node(s) do
    if isempty(s)
      return T(undef, [0 for _ in 1:N]...)
    end
    len = length(s.obj)
    m = [load_object(s, TypeParams(Array{S, N - 1}, sub_tp), i) for i in 1:len]
    return stack(m; dims=1)
  end
end

function load_object(s::DeserializerState, T::Type{Array{S, N}}, key::Union{Int, Symbol}) where {S, N}
  return load_node(s, key) do
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
      for (i, param_tp) in enumerate(parameters(tp))
        save_type_params(s, param_tp)
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type{Tuple})
  !haskey(s, :params) && return TypeParams(T{}, nothing)
  return load_node(s, :params) do
    tuple_params = load_array_node(s) do _
      load_type_params(s, decode_type(s))
    end
    tuple_types = Tuple([type(x) for x in tuple_params])
    raw_params = Tuple(parameters(x) for x in tuple_params)
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
  p = parameters(tp)
  entries = load_array_node(s) do i
    S = fieldtype(T, i)
    if serialize_with_id(S)
      return load_ref(s)
    else
      elem_tp = p[i] isa TypeParams ? p[i] : TypeParams(S, p[i])
      return load_object(s, elem_tp)
    end
  end
  return Tuple(entries)
end

function load_object(s::DeserializerState, T::Type{<:Tuple})
  entries = load_array_node(s) do i
    S = fieldtype(T, i)
    load_object(s, S)
  end
  return T(entries)
end

################################################################################
# Saving and loading NamedTuple
@register_serialization_type NamedTuple

function type_params(obj::T) where T <: NamedTuple
  return TypeParams(T, (x.first => type_params(x.second) for x in pairs(obj))...)
end

function load_type_params(s::DeserializerState, T::Type{NamedTuple})
  return load_node(s, :params) do
    pairs_list = Pair{Symbol, Any}[]
    for k in propertynames(s.obj)
      tp = load_node(s, k) do
        is_string(s) ? TypeParams(decode_type(s), nothing) : load_type_params(s, decode_type(s))
      end
      push!(pairs_list, k => tp)
    end
    names = Tuple(Symbol(k) for (k, _) in pairs_list)
    types = Tuple([type(v) for (_, v) in pairs_list])
    return TypeParams(T{names, Tuple{types...}}, pairs_list...)
  end
end

function save_object(s::SerializerState, obj::NamedTuple)
  save_object(s, values(obj))
end

function load_object(s::DeserializerState, tp::TypeParams{<:NamedTuple, <:Tuple{Vararg{Pair}}})
  T = type(tp)
  tp_values = Tuple(v for (_, v) in parameters(tp))
  return T(load_object(s, TypeParams(Tuple{Base.fieldtypes(T)...}, tp_values)))
end

################################################################################
# Saving and loading dicts
@register_serialization_type Dict

function type_params(obj::T) where {S <: Union{Symbol, Int, String},
                                    T <: Dict{S, Any}}
  return TypeParams(
    T,
    :key_params => TypeParams(S, nothing),
    :value_params => Tuple((x.first => type_params(x.second) for x in pairs(obj)))
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
  tp::TypeParams{Dict{S, T}, <:Tuple{Pair, Pair}}) where {T, S <: Union{Symbol, Int, String}}
  save_data_dict(s) do
    save_object(s, encode_type(Dict), :name)
    save_data_dict(s, :params) do
      isempty(parameters(tp)) && save_object(s, encode_type(T), :value_params)
      save_type_params(s, tp[:key_params], :key_params)

      if tp[:value_params] isa Tuple
        save_data_dict(s, :value_params) do
          for (k, param_tp) in tp[:value_params]
            save_type_params(s, param_tp, Symbol(k))
          end
        end
      else
        save_type_params(s, tp[:value_params], :value_params)
      end
    end
  end
end

function load_type_params(s::DeserializerState, T::Type{Dict})
  return load_node(s, :params) do
    key_tp = load_node(s, :key_params) do
      is_string(s) ? TypeParams(decode_type(s), nothing) : load_type_params(s, decode_type(s))
    end
    S = type(key_tp)

    load_node(s, :value_params) do
      if is_object(s) && !haskey(s, :name) && !haskey(s, type_key)
        # Heterogeneous Dict{S, Any} — per-key type params stored as dict
        pairs_list = Pair[]
        for k in propertynames(s.obj)
          key = S <: Integer ? parse(Int, string(k)) : S(k)
          tp = load_node(s, Symbol(k)) do
            load_type_params(s, decode_type(s))
          end
          push!(pairs_list, key => tp)
        end
        TypeParams(Dict{S, Any}, :key_params => key_tp, :value_params => Tuple(pairs_list))
      else
        value_tp = load_type_params(s, decode_type(s))
        TypeParams(Dict{S, type(value_tp)}, :key_params => key_tp, :value_params => value_tp)
      end
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
                     T::Type{<:Dict{S, U}}) where {S <: Union{Symbol, String, Int}, U}
  load_object(s, TypeParams(T, nothing))
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Dict{S, Any}, <:Dict}) where {S <: Union{Symbol, String, Int}}
  T = type(tp)
  p = parameters(tp)
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

function load_object(s::DeserializerState, tp::TypeParams{<:Dict{S, U}, Nothing}) where {S, U}
  T = type(tp)
  pairs = load_array_node(s) do _
    load_object(s, TypeParams(Tuple{S, U}, (nothing, nothing)))
  end::Vector{Tuple{S, U}}
  return T(k => v for (k, v) in pairs)
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Dict{S, U}, <:Tuple{Pair, Pair}}) where {S <: Union{Int, String, Symbol}, U}
  T = type(tp)
  if tp[:value_params] isa Tuple
    # Heterogeneous Dict{S, Any} — per-key type params
    per_key_params = Dict(k => v for (k, v) in tp[:value_params])
    dict = Dict{S, Any}()
    for k in propertynames(s.obj)
      key = S <: Integer ? parse(S, string(k)) : S(k)
      dict[key] = load_object(s, per_key_params[key], k)
    end
    return Dict{S, Any}(dict)
  else
    # Homogeneous Dict{S, U} — single value type
    pairs = Tuple{S, U}[]
    for k in propertynames(s.obj)
      key = S <: Integer ? parse(S, string(k)) : S(k)
      push!(pairs, (key, load_object(s, tp[:value_params], k)))
    end
    isempty(pairs) && return T()
    _, v = first(pairs)
    return Dict{S, typeof(v)}(k => v for (k, v) in pairs)
  end
end

function load_object(s::DeserializerState,
                     tp::TypeParams{<:Dict{S, U}, <:Tuple{Vararg{Pair}}}) where {S, U}
  T = type(tp)
  pairs = load_array_node(s) do _
    load_object(s, TypeParams(Tuple{S, U}, (tp[:key_params], tp[:value_params])))
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
  return load_node(s, :params) do
    sub_tp = load_type_params(s, decode_type(s))
    return TypeParams(T{type(sub_tp)}, parameters(sub_tp))
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
  p = parameters(tp)
  elem_tp = p isa TypeParams ? p : TypeParams(T, p)
  elems = load_array_node(s) do _
    load_object(s, elem_tp)
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
  ring = parameters(tp)
  pos = Int[]
  entry_type = elem_type(ring)
  values = entry_type[]
  load_array_node(s) do _
    push!(pos, load_object(s, Int, 1))
    push!(values, load_object(s, TypeParams(entry_type, ring), 2))
  end
  return sparse_row(ring, pos, values)
end
