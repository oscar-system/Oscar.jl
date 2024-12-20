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

###############################################################################
#
#   Weight lattices
#
###############################################################################

@register_serialization_type WeightLattice uses_id

function save_object(s::SerializerState, P::WeightLattice)
  save_data_dict(s) do
    save_typed_object(s, root_system(P), :root_system)
  end
end

function load_object(s::DeserializerState, ::Type{WeightLattice})
  R = load_typed_object(s, :root_system)
  return weight_lattice(R)
end

@register_serialization_type WeightLatticeElem uses_params

function save_object(s::SerializerState, w::WeightLatticeElem)
  save_object(s, _vec(coefficients(w)))
end

function load_object(s::DeserializerState, ::Type{WeightLatticeElem}, P::WeightLattice)
  return WeightLatticeElem(P, load_object(s, Vector, ZZ))
end

function save_type_params(s::SerializerState, w::WeightLatticeElem)
  save_data_dict(s) do
    save_object(s, encode_type(typeof(w)), :name)
    parent_w = parent(w)
    parent_ref = save_as_ref(s, parent_w)
    save_object(s, parent_ref, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{WeightLatticeElem})
  return load_typed_object(s)
end

###############################################################################
#
#   Weyl groups
#
###############################################################################

@register_serialization_type WeylGroup uses_id

function save_object(s::SerializerState, W::WeylGroup)
  save_data_dict(s) do
    save_typed_object(s, root_system(W), :root_system)
  end
end

function load_object(s::DeserializerState, ::Type{WeylGroup})
  R = load_typed_object(s, :root_system)
  return weyl_group(R)
end

@register_serialization_type WeylGroupElem uses_params

function save_object(s::SerializerState, x::WeylGroupElem)
  save_object(s, word(x))
end

function load_object(s::DeserializerState, ::Type{WeylGroupElem}, W::WeylGroup)
  return W(load_object(s, Vector, UInt8); normalize=false)
end

function save_type_params(s::SerializerState, x::WeylGroupElem)
  save_data_dict(s) do
    save_object(s, encode_type(typeof(x)), :name)
    parent_x = parent(x)
    parent_ref = save_as_ref(s, parent_x)
    save_object(s, parent_ref, :params)
  end
end

function load_type_params(s::DeserializerState, ::Type{WeylGroupElem})
  return load_typed_object(s)
end
