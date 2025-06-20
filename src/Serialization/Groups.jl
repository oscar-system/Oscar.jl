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
#   In those cases where the object identity of `GAPWrap.FamilyObj`
#   objects must be preserved in the deserialization, we cannot simply
#   leave it to the (de)serialization of Oscar objects to deal with "their"
#   underlying GAP objects.
#
#   For example, serialize two subgroups `H`, `K` of a finitely presented
#   group `G` in Oscar.
#   In order to achieve that the deserialized `H` and `K` in a new Julia
#   session are compatible (that is, elements of `H` and elements of `K`
#   can be multiplied with each other), we have to make sure that
#   `GAPWrap.FamilyObj(GapObj(H)) === GAPWrap.FamilyObj(GapObj(K))`
#   are identical.
#   Since `H` and `K` know about these objects only via `GapObj(H)` and
#   `GapObj(K)`,
#   the mechanism that automatically takes care of object identity in the
#   (de)serialization can be used only if (de)serialization methods for
#   `GapObj(H)` and `GapObj(K)` are provided.
#
#   (The situation would be different if `H` and `K` would store references
#   to the "full group" `G`.
#   In this case, we could force that serializing `H` and `K` involves also
#   a serialization of `G`, then deserializing `H` and `K` would recreate the
#   common `G`, and the object identity of the GAP family object could be
#   forced via `G`.)

import Oscar: GAPGroup, _coeff

##############################################################################
# `GAPGroupElem` objects get serialized together with their parents.
const GrpElemUnionType = Union{GAPGroupElem, FinGenAbGroupElem}

type_params(p::T) where T <: GrpElemUnionType = TypeParams(T, parent(p))

#############################################################################
# attributes handling
const GAPGroup_attributes = [
  :order, :is_abelian, :is_nilpotent, :is_perfect, :is_simple, :is_solvable
]

function save_attrs(s::SerializerState, G::T) where T <: GAPGroup
  save_data_dict(s, :attrs) do 
    for attr in attrs_list(T)
      func = Symbol(string("has_", attr))
      if @eval $func($G)
        attr_value = @eval $attr($G)
        save_typed_object(s, attr_value, attr)
      end
    end
  end
end

function load_attrs(s::DeserializerState, G::T) where T <: GAPGroup
  !with_attrs(s) && return
  haskey(s, :attrs) && load_node(s, :attrs) do d
    for attr in attrs_list(T)
      if haskey(d, attr)
        func = Symbol(string("set_", attr))
        attr_value = load_typed_object(s, attr)
        @eval $func($G, $attr_value)
      end
    end
  end
end

##############################################################################
# PermGroup

@register_serialization_type PermGroup uses_id [GAPGroup_attributes;]

function save_object(s::SerializerState, G::PermGroup)
  n = degree(G)
  save_data_dict(s) do
    save_object(s, n, :degree)
    save_object(s, [Vector{Int}(GAPWrap.ListPerm(GapObj(x))) for x in gens(G)], :gens)
  end
end

function load_object(s::DeserializerState, ::Type{PermGroup})
  n = load_object(s, Int, :degree)
  generators = load_object(s, Vector{Vector{Int}}, :gens)
  G = permutation_group(n, [perm(x) for x in generators])
  return G
end


##############################################################################
# PermGroupElem
# We need just the array of image points up to the largest moved point.
# (Or would it be easier for other systems if the length of the array
# is equal to the degree of the parent group?)

@register_serialization_type PermGroupElem

type_params(x::T) where T <: GroupElem = TypeParams(T, parent(x))

function save_object(s::SerializerState, p::PermGroupElem)
  save_object(s, Vector{Int}(GAPWrap.ListPerm(GapObj(p))))
end

function load_object(s::DeserializerState, T::Type{PermGroupElem}, parent_group::PermGroup)
  imgs = load_object(s, Vector{Int})
  return perm(parent_group, imgs)
end


##############################################################################
# FPGroup, SubFPGroup
# PcGroup, SubPcGroup
# We do the same for full free groups, subgroups of free groups,
# full f.p. groups, and subgroups of f.p. groups.

@register_serialization_type FPGroup uses_id
@register_serialization_type SubFPGroup uses_id
@register_serialization_type PcGroup uses_id
@register_serialization_type SubPcGroup uses_id

function type_params(G::T) where T <: Union{FPGroup, SubFPGroup, PcGroup, SubPcGroup}
  TypeParams(T, GapObj(G))
end

function save_object(s::SerializerState,
                     G::T) where T <: Union{FPGroup, SubFPGroup, PcGroup, SubPcGroup}
  #needs place holder
  save_data_array(() -> (),  s)
end

function load_object(s::DeserializerState, ::Type{T},
                     G::GapObj) where T <: Union{FPGroup, SubFPGroup, PcGroup, SubPcGroup}
  return T(G)
end

##############################################################################
# FPGroupElem, SubFPGroupElem
# We need the parent and a description of the word that defines the element.

@register_serialization_type FPGroupElem
@register_serialization_type SubFPGroupElem

function save_object(s::SerializerState, g::Union{FPGroupElem, SubFPGroupElem})
  save_object(s, Vector{Int}(vcat([[x[1], x[2]] for x in syllables(g)]...)))
end

typecombinations = (
    (:FPGroupElem, :FPGroup),
    (:SubFPGroupElem, :SubFPGroup),
)

for (eltype, type) in typecombinations
  @eval function load_object(s::DeserializerState, ::Type{$eltype}, parent_group::$type)
    lo = load_object(s, Vector{Int})
    fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(GapObj(parent_group)))
    if GAPWrap.IsElementOfFpGroupFamily(fam)
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
end

##############################################################################
# PcGroupElem
# We need the parent and a description of the exponent vector
# that defines the element.

@register_serialization_type PcGroupElem
@register_serialization_type SubPcGroupElem

function save_object(s::SerializerState, g::Union{PcGroupElem, SubPcGroupElem})
  elfam = GAPWrap.FamilyObj(GapObj(g))
  fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)
  save_object(s, Vector{Int}(GAPWrap.ExponentsOfPcElement(fullpcgs, GapObj(g))))
end

typecombinations = (
    (:PcGroupElem, :PcGroup),
    (:SubPcGroupElem, :SubPcGroup),
)

for (eltype, type) in typecombinations
  @eval function load_object(s::DeserializerState, ::Type{$eltype}, parent_group::$type)
    lo = load_object(s, Vector{Int})
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(GapObj(parent_group)))
    fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)::GapObj
    gapelm = GAPWrap.PcElementByExponentsNC(fullpcgs, GapObj(lo, true))::GapObj
    return Oscar.group_element(parent_group, gapelm)
  end
end


##############################################################################
# FinGenAbGroup

@register_serialization_type FinGenAbGroup uses_id

function save_object(s::SerializerState, G::FinGenAbGroup)
  save_object(s, rels(G))
end

function load_object(s::DeserializerState, ::Type{FinGenAbGroup})
  return abelian_group(load_object(s, Matrix{ZZRingElem}))
end

# elems
@register_serialization_type FinGenAbGroupElem

function save_object(s::SerializerState, g::FinGenAbGroupElem)
  save_object(s, _coeff(g))
end

function load_object(s::DeserializerState, ::Type{FinGenAbGroupElem}, G::FinGenAbGroup)
  return G(vec(load_object(s, Matrix{ZZRingElem})))
end
