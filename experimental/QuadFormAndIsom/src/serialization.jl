## TODO
# it might be beneficial to eventually implement a save_type_params / load_type_params
# especially if one wants to serialize vectors of these objects
# also ZZLatWithIsom currently always serialized with ambient space, this may not always be
# desired

@register_serialization_type QuadSpaceWithIsom

function save_object(s::SerializerState, QS::QuadSpaceWithIsom)
  save_data_dict(s) do
    save_typed_object(s, space(QS), :quad_space)
    save_typed_object(s, isometry(QS), :isom)

    save_object(s, order_of_isometry(QS), :order)
  end
end

function load_object(s::DeserializerState, ::Type{QuadSpaceWithIsom})
  quad_space = load_typed_object(s, :quad_space)
  isom = load_typed_object(s, :isom)

  # not quite sure how to deal with IntExt/PosInf yet..
  # we could add it to the basic type section
  # then we could use
  # load_object(s, IntExt, dict[:order]) 
  n = load_object(s, Int, :order)
  
  return QuadSpaceWithIsom(quad_space, isom, n)
end

@register_serialization_type ZZLatWithIsom

function save_object(s::SerializerState, L::ZZLatWithIsom)
  save_data_dict(s) do
    save_typed_object(s, ambient_space(L), :ambient_space)
    save_typed_object(s, basis_matrix(L), :basis)
  end
end

function load_object(s::DeserializerState, ::Type{ZZLatWithIsom})
  quad_space = load_typed_object(s, :ambient_space)
  B = load_typed_object(s, :basis)
  return lattice(quad_space, B; check=false)
end
