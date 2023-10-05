
##############################################################################
#
# General ideas of the (de)serialization of groups and group elements:
#
# - Each `GAPGroupElem` object gets serialized together with its `parent`.
#
# - We request `uses_id` for the (de)serialization of group objects
#   because parent references to the same group must point to the same object.
#
# - We do not request `uses_id` for the (de)serialization of group elements,
#   deserializing the same group element twice may yield two nonidentical
#   objects.
#
# - There are cases where the comparison of Oscar group elements relies on
#   the comparison of underlying GAP elements,
#   such that the object identity of the `GAPWrap.FamilyObj` values of these
#   GAP elements is required.
#   `FPGroupElem` and `PcGroupElem` are examples.
#   In these cases, we install a (de)serialization of the GAP group object,
#   with requirement `uses_id`.
#   This way, we can (de)serialize Oscar group elements and Oscar groups
#   as well as GAP groups.
#
#   (Note that from GAP's viewpoint, it would look more logical
#   to (de)serialize the `GAPWrap.FamilyObj` object,
#   because the group element in GAP stores a reference to this object
#   but not to a group object.
#   However, technically the creation of the family objects in GAP
#   is encapsulated inside functions such as `GAPWrap.FreeGroup`
#   that create groups.
#   Thus it would require extra code on the Oscar side if we would want to
#   first deserialize some elements of a free group by creating the family
#   object on the GAP side, and then later create a free group with prescribed
#   elements family.
#   In particular, we do not support the (de)serialization of GAP's
#   family objects.)


##############################################################################
# `GAPGroupElem` objects get serialized together with their parents.

function save_type_params(s::SerializerState, p::T) where T <: GAPGroupElem
  # this has just been more or less copied from the Rings section
  # and might be removed from this file during a future refactor
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)
    parent_p = parent(p)
    if serialize_with_id(parent_p)
      parent_ref = save_as_ref(s, parent_p)
      save_object(s, parent_ref, :params)
    else
      save_typed_object(s, parent_p, :params)
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:GAPGroupElem},
                          dict::Dict{Symbol, Any})
  return load_typed_object(s, dict)
end


##############################################################################
# PermGroup

@register_serialization_type PermGroup uses_id

function save_object(s::SerializerState, G::PermGroup)
  n = degree(G)

  save_data_dict(s) do
    save_object(s, n, :degree)

    save_data_array(s, :gens) do
      for x in gens(G)
        # We don't need any type information inside the array.
        save_data_array(s) do
          for i in Vector{Int}(GAP.Globals.ListPerm(x.X))
            save_object(s, i)
          end
        end
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{PermGroup}, dict::Dict)
  n = parse(Int, dict[:degree])
  generators = load_object(s, Vector, dict[:gens], (Vector{Int}, Int))

  return permutation_group(n, [perm(x) for x in generators])
end


################################################################################
# PermGroupElem
# We need just the array of image points up to the largest moved point.
# (Or would it be easier for other systems if the length of the array
# is equal to the degree of the parent group?)

@register_serialization_type PermGroupElem uses_params

function save_object(s::SerializerState, p::PermGroupElem)
  vector_int = Vector{Int}(GAP.Globals.ListPerm(p.X))
  # again, here we dont need to store as a vector
  save_data_array(s) do
    for i in vector_int
      save_object(s, i)
    end
  end
end

function load_object(s::DeserializerState, T::Type{PermGroupElem},
                     imgs_data::Vector,
                     parent_group::PermGroup)
  imgs = load_object(s, Vector, imgs_data, Int)
  return perm(parent_group, imgs)
end


################################################################################
# FPGroup
# We do the same for full free groups, subgroups of free groups,
# full f.p. groups, and subgroups of f.p. groups.

@register_serialization_type FPGroup uses_id

function save_object(s::SerializerState, G::FPGroup)
  save_data_dict(s) do
    save_object(s, G.X, :X)
  end
end

function load_object(s::DeserializerState, ::Type{FPGroup}, dict::Dict)
  return FPGroup(load_object(s, GapObj, dict[:X]))
end


################################################################################
# FPGroupElem
# We need the parent and a description of the word that defines the element.

@register_serialization_type FPGroupElem uses_params

function save_object(s::SerializerState, g::FPGroupElem)
  vector_int = Vector{Int}(vcat([[x[1], x[2]] for x in syllables(g)]...))
  # again, here we don't need to store as a vector.
  save_data_array(s) do
    for i in vector_int
      save_object(s, i)
    end
  end
end

function load_object(s::DeserializerState, ::Type{FPGroupElem},
                     word_data::Vector,
                     parent_group::FPGroup)
  lo = load_object(s, Vector, word_data, Int)
  fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(parent_group.X))
  if GAP.Globals.IsElementOfFpGroupFamily(fam)
    # go via the underlying free group
    free = GAP.getbangproperty(fam, :freeGroup)
    freefam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(free))
    freeelm = GAPWrap.ObjByExtRep(freefam, GapObj(lo, true))
    gapelm = GAPWrap.ElementOfFpGroup(fam, freeelm)
  else
    # element in a free group
    gapelm = GAPWrap.ObjByExtRep(fam, GapObj(lo, true))
  end
  return Oscar.group_element(parent_group, gapelm)
end
