###############################################################################
#
#   Root systems
#
###############################################################################

@register_serialization_type RootSystem uses_id

function save_object(s::SerializerState, R::RootSystem)
  save_data_dict(s) do
    save_object(s, cartan_matrix(R), :cartan_matrix)
    if has_root_system_type(R)
      type, type_ordering = root_system_type_with_ordering(R)
      save_object(s, type, :type)
      if !issorted(type_ordering) # don't save if it's the default
        save_object(s, type_ordering, :type_ordering)
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{RootSystem})
  cm = load_object(s, Matrix, Int, :cartan_matrix)
  R = root_system(cm; check=false, detect_type=false)
  if haskey(s, :type)
    type = Vector{Tuple{Symbol,Int}}(load_object(s, Vector, (Tuple, [Symbol, Int]), :type)) # coercion needed due to https://github.com/oscar-system/Oscar.jl/issues/3983
    if haskey(s, :type_ordering)
      type_ordering = load_object(s, Vector, Int, :type_ordering)
      set_root_system_type!(R, type, type_ordering)
    else
      set_root_system_type!(R, type)
    end
  end
  return R
end

@register_serialization_type RootSpaceElem uses_params
@register_serialization_type DualRootSpaceElem uses_params

function save_object(s::SerializerState, r::Union{RootSpaceElem,DualRootSpaceElem})
  save_object(s, _vec(coefficients(r)))
end

function load_object(
  s::DeserializerState, T::Type{<:Union{RootSpaceElem,DualRootSpaceElem}}, R::RootSystem
)
  return T(R, load_object(s, Vector, QQ))
end

function save_type_params(s::SerializerState, r::Union{RootSpaceElem,DualRootSpaceElem})
  save_data_dict(s) do
    save_object(s, encode_type(typeof(r)), :name)
    rs_x = root_system(r)
    rs_ref = save_as_ref(s, rs_x)
    save_object(s, rs_ref, :params)
  end
end

function load_type_params(
  s::DeserializerState, ::Type{<:Union{RootSpaceElem,DualRootSpaceElem}}
)
  return load_typed_object(s)
end
