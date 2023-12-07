
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
# - Not all subobjects of GAP objects get serialized, only the ones that are
#   needed to define the object.
#   For example, we do not attempt to serialize the known properties and
#   attributes.
#   (We may decide that some of them shall better get serialized because
#   they are expensive to compute, such as finiteness and group order of
#   finitely presented groups.)
#   Moreover, some subobjects are intended just for "internal purposes".
#   For example, a pc group in GAP stores (via its family object)
#   a rewriting system in terms of elements of a free group.
#   We do not serialize this information, it gets created anew when the
#   pc group gets deserialized; this way, the original pc group and the
#   deserialized one store nonidentical free group objects.
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
#
# - Currently it is not our aim to provide a (de)serialization of as many
#   GAP objects as possible.
#   We provide methods for those GAP objects that are needed for the
#   (de)serialization of Oscar objects.
#   For example, we need a serialization of free groups in GAP, because of
#   the object identity requirements, but (de)serializing elements of
#   free groups in Oscar is done by serializing the underlying word and
#   the parent object, without serializing the underlying element in GAP.
#
# - Remark:
#   In those cases where the object identity of `GAP.Globals.FamilyObj`
#   objects must be preserved in the deserialization, we cannot simply
#   leave it to the (de)serialization of Oscar objects to deal with "their"
#   underlying GAP objects.
#
#   For example, serialize two subgroups `H`, `K` of a finitely presented
#   group `G` in Oscar.
#   In order to achieve that the deserialized `H` and `K` in a new Julia
#   session are compatible (that is, elements of `H` and elements of `K`
#   can be multiplied with each other), we have to make sure that
#   `GAP.Globals.FamilyObj(H.X) === GAP.Globals.FamilyObj(K.X)` are
#   identical.
#   Since `H` and `K` know about these objects only via `H.X` and `K.X`,
#   the mechanism that automatically takes care of object identity in the
#   (de)serialization can be used only if (de)serialization methods for `H.X`
#   and `K.X` are provided.
#
#   (The situation would be different if `H` and `K` would store references
#   to the "full group" `G`.
#   In this case, we could force that serializing `H` and `K` involves also
#   a serialization of `G`, then deserializing `H` and `K` would recreate the
#   common `G`, and the object identity of the GAP family object could be
#   forced via `G`.)


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

function load_type_params(s::DeserializerState, ::Type{<:GAPGroupElem})
  return load_typed_object(s)
end


##############################################################################
# PermGroup

@register_serialization_type PermGroup uses_id

function save_object(s::SerializerState, G::PermGroup)
  n = degree(G)

  save_data_dict(s) do
    save_object(s, n, :degree)
    save_object(s, [Vector{Int}(GAP.Globals.ListPerm(x.X)) for x in gens(G)], :gens)
  end
end

function load_object(s::DeserializerState, ::Type{PermGroup})
  n = load_object(s, Int, :degree)
  generators = load_object(s, Vector, (Vector{Int}, Int), :gens)

  return permutation_group(n, [perm(x) for x in generators])
end


##############################################################################
# PermGroupElem
# We need just the array of image points up to the largest moved point.
# (Or would it be easier for other systems if the length of the array
# is equal to the degree of the parent group?)

@register_serialization_type PermGroupElem uses_params

function save_object(s::SerializerState, p::PermGroupElem)
  save_object(s, Vector{Int}(GAP.Globals.ListPerm(p.X)))
end

function load_object(s::DeserializerState, T::Type{PermGroupElem}, parent_group::PermGroup)
  imgs = load_object(s, Vector, Int)
  return perm(parent_group, imgs)
end


##############################################################################
# FPGroup
# We do the same for full free groups, subgroups of free groups,
# full f.p. groups, and subgroups of f.p. groups.

@register_serialization_type FPGroup uses_id

function save_object(s::SerializerState, G::FPGroup)
  save_data_dict(s) do
    save_object(s, G.X, :X)
  end
end

function load_object(s::DeserializerState, ::Type{FPGroup})
  return FPGroup(load_object(s, GapObj, :X))
end


##############################################################################
# FPGroupElem
# We need the parent and a description of the word that defines the element.

@register_serialization_type FPGroupElem uses_params

function save_object(s::SerializerState, g::FPGroupElem)
  save_object(s, Vector{Int}(vcat([[x[1], x[2]] for x in syllables(g)]...)))
end

function load_object(s::DeserializerState, ::Type{FPGroupElem}, parent_group::FPGroup)
  lo = load_object(s, Vector, Int)
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


##############################################################################
# PcGroup

@register_serialization_type PcGroup uses_id

function save_object(s::SerializerState, G::PcGroup)
  save_data_dict(s) do
    save_object(s, G.X, :X)
  end
end

function load_object(s::DeserializerState, ::Type{PcGroup})
  return PcGroup(load_object(s, GapObj, :X))
end


##############################################################################
# PcGroupElem
# We need the parent and a description of the exponent vector
# that defines the element.

@register_serialization_type PcGroupElem uses_params

function save_object(s::SerializerState, g::PcGroupElem)
  elfam = GAPWrap.FamilyObj(g.X)
  fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)
  save_object(s, Vector{Int}(GAP.Globals.ExponentsOfPcElement(fullpcgs, g.X)))
end

function load_object(s::DeserializerState, ::Type{PcGroupElem}, parent_group::PcGroup)
  lo = load_object(s, Vector, Int)
  elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(parent_group.X))
  fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)::GapObj
  gapelm = GAP.Globals.PcElementByExponentsNC(fullpcgs, GapObj(lo, true))::GapObj
  return Oscar.group_element(parent_group, gapelm)
end
