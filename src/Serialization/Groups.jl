##############################################################################
#
# General ideas of the (de)serialization of groups and group elements:
#
# - Each `GAPGroupElem` object gets serialized together with its `parent`,
#   via the `type_params(x::T)` method for `T <: SetElem`.
#
# - We request `uses_id` for the (de)serialization of group objects
#   because parent references to the same group must point to the same object.
#
# - We do not request `uses_id` for the (de)serialization of group elements,
#   deserializing the same group element twice may yield two nonidentical
#   objects.
#
# - Not all subobjects of groups get serialized, only the ones that are
#   needed to define the object.
#   In particular, we do not attempt to serialize all known properties and
#   attributes.
#   (We may decide that some of them shall better get serialized because
#   they are expensive to compute, such as finiteness and group order of
#   finitely presented groups.)
#
# - For Oscar objects that rely on underlying GAP objects,
#   we do not implement the serialization recursively
#   via a serialization of these GAP objects.
#   For example, consider the situation that two groups in Oscar point to
#   the same GAP object.
#   Serializing and then deserializing these two groups may yield Oscar groups
#   that point to different GAP objects.
#
#   There are cases where the equality check for Oscar group elements
#   relies on the equality check for the underlying GAP elements,
#   such that the object identity of the `GAPWrap.FamilyObj` values of these
#   GAP elements is required for the equality of the elements.
#   In these cases, we assume that the analogous object identity is satisfied
#   also for the Oscar group objects in question.
#   Thus the (de)serialization code for groups and group elements does not
#   rely on a (de)serialization for the underlying GAP objects.
#
#   The type `FPGroupElem` and `PcGroupElem` are examples.
#   Serializing a tuple containing several subgroups and quotients of a
#   finitely presented group stores the information about the full group
#   and the underlying full free group,
#   such that the deserialization yields a tuple of groups which have the
#   same relations to each other as the groups in the original tuple.
#

using Oscar: GAPGroup, _coeff

##############################################################################
# `GAPGroupElem` objects get serialized together with their parents.
const GrpElemUnionType = Union{GAPGroupElem, FinGenAbGroupElem}

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

function save_object(s::SerializerState, p::PermGroupElem)
  save_object(s, Vector{Int}(GAPWrap.ListPerm(GapObj(p))))
end

function load_object(s::DeserializerState, T::Type{PermGroupElem}, parent_group::PermGroup)
  imgs = load_object(s, Vector{Int})
  return perm(parent_group, imgs)
end


##############################################################################
# FPGroup

@register_serialization_type FPGroup uses_id

function type_params(G::FPGroup)
  F = free_group(G)
  if F === G
    # A full *free* group does not need type parameters.
    return TypeParams(FPGroup, nothing)
  else
    # A group with relators needs the underlying free group
    # as a type parameter,
    # since the same free group object can be used for creating
    # several f.p. groups.
    return TypeParams(FPGroup, F)
  end
end

# There are several internal representations of elements of free groups in GAP.
# Let the serialized data not depend on the names of the GAP filters.
_word_filters = Dict(
  :IsSyllableWordsFamily => :syllables,
  :IsLetterWordsFamily => :letters,
  :IsWLetterWordsFamily => :wletters,
  :IsBLetterWordsFamily => :bletters,
  )
_word_filters_inv = Dict()
for x in keys(_word_filters)
  _word_filters_inv[_word_filters[x]] = x
end

function save_object(s::SerializerState, G::FPGroup)
  F = free_group(G)
  if F === G
    # A full *free* group is represented by the names of the generators,
    # and the information about the GAP filter that defines the internal
    # representation of the elements; there are four such filters.
    # (Currently we do not support free groups on infinitely many generators.)
    X = GapObj(G)
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
    @assert GAP.Globals.HasIsWholeFamily(X) && GAPWrap.IsWholeFamily(X)
    save_data_dict(s) do
      # (rank and) names of generators
      Xnames = GAP.getbangproperty(elfam, :names)::GapObj
      names = Vector{String}(Xnames)
      save_object(s, names, :names)
      # the internal representation of elements
      if GAP.Globals.IsSyllableWordsFamily(elfam)::Bool
        wfilt = :IsSyllableWordsFamily
      elseif GAP.Globals.IsLetterWordsFamily(elfam)::Bool
        wfilt = :IsLetterWordsFamily
      elseif GAP.Globals.IsWLetterWordsFamily(elfam)::Bool
        wfilt = :IsWLetterWordsFamily
      elseif GAP.Globals.IsBLetterWordsFamily(elfam)::Bool
        wfilt = :IsBLetterWordsFamily
      else
        error("not supported internal representation")
      end
      save_object(s, _word_filters[wfilt], :rep)
    end
  else
    # A group with relators is represented by the underlying free group
    # (via the type parameter) and a description of the relators.
    X = GapObj(G)
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
    relators = GAP.getbangproperty(elfam, :relators)::GapObj
    save_data_dict(s) do
      save_object(s, [Vector{Int}(GAPWrap.ExtRepOfObj(x)) for x in relators], :relators)
    end
  end
end

function load_object(s::DeserializerState, ::Type{FPGroup}, ::Nothing)
  # Without type parameters, the object is a free group.
  load_node(s) do d
    # Create a new full free group.
    wfilt = getproperty(GAP.Globals, _word_filters_inv[load_object(s, Symbol, :rep)])::GapObj

    init = load_node(s, :names) do names
      GapObj(names; recursive = true)
    end
    G = GAPWrap.FreeGroup(wfilt, init)
    res = FPGroup(G)
    res.free_group = res
    return res
  end
end

function load_object(s::DeserializerState, ::Type{FPGroup}, F::FPGroup)
  # The object is a f.p. group, `F` is the underlying free group.
  # Create a new full f.p. group.
  FX = GapObj(F)
  relators = load_object(s, Vector{Vector{Int}}, :relators)
  elfreefam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(FX))
  rels = [GAPWrap.ObjByExtRep(elfreefam, GapObj(x, true)) for x in relators]
  G = FX / GapObj(rels)::GapObj
  return fp_group(G, F)
end


##############################################################################
# PcGroup
# The full group is described by its collector, no type parameters are needed.

@register_serialization_type PcGroup uses_id

function save_object(s::SerializerState, G::PcGroup)
  X = GapObj(G)
  elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
  fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)::GapObj
  @assert fullpcgs === GAPWrap.Pcgs(X)
  save_data_dict(s) do
    # relative orders
    relord = [GAPWrap.RelativeOrderOfPcElement(fullpcgs, x) for x in fullpcgs]
    save_object(s, relord, :relord)
    # power relators
    rels = Tuple{Int, Vector{Int}}[]
    for i in 1:length(relord)
      ne = fullpcgs[i]^relord[i]
      if ! GAPWrap.IsOne(ne)
        push!(rels, (i, _free_group_extrep_from_exponents(
          Vector{Int}(GAPWrap.ExponentsOfPcElement(fullpcgs, ne)))))
      end
    end
    save_object(s, rels, :power_rels)
    # commutator relators
    rels = Tuple{Int, Int, Vector{Int}}[]
    for i in 1:(length(relord)-1)
      for j in (i+1):length(relord)
        ne = GAP.Globals.Comm(fullpcgs[j], fullpcgs[i])::GapObj
        if ! GAPWrap.IsOne(ne)
          push!(rels, (j, i, _free_group_extrep_from_exponents(
            Vector{Int}(GAPWrap.ExponentsOfPcElement(fullpcgs, ne)))))
        end
      end
    end
    save_object(s, rels, :comm_rels)
  end
end

function load_object(s::DeserializerState, ::Type{PcGroup}, ::Nothing)
  relord = load_object(s, Vector{Int}, :relord)
  F = GAPWrap.FreeGroup(GAP.Globals.IsSyllableWordsFamily, length(relord))
  fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(F))
  rws = GAP.Globals.SingleCollector(F, GapObj(relord))::GapObj
  for (i, elm) in load_object(s, Vector{Tuple{Int, Vector{Int}}}, :power_rels)
    GAP.Globals.SetPower(rws, i, GAPWrap.ObjByExtRep(fam, GapObj(elm)))
  end
  for (j, i, elm) in load_object(s, Vector{Tuple{Int, Int, Vector{Int}}}, :comm_rels)
    GAP.Globals.SetCommutator(rws, j, i, GAPWrap.ObjByExtRep(fam, GapObj(elm)))
  end
  G = GAP.Globals.GroupByRwsNC(rws)::GapObj
  return PcGroup(G)
end


##############################################################################
# SubFPGroup, SubPcGroup
# The group is described by full group (via a reference) and generators.

@register_serialization_type SubFPGroup uses_id
@register_serialization_type SubPcGroup uses_id

function type_and_params(G::T) where T <: Union{SubFPGroup, SubPcGroup}
  # The subgroup needs a reference to its full group.
  return TypeParams(T, G.full_group)
end

function save_object(s::SerializerState, G::SubFPGroup)
  save_data_dict(s) do
    save_object(s, [Vector{Int}(GAPWrap.ExtRepOfObj(GapObj(x)::GapObj)::GapObj) for x in gens(G)], :gens)
  end
end

function save_object(s::SerializerState, G::SubPcGroup)
  X = GapObj(G)
  elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
  fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)::GapObj

  save_data_dict(s) do
    save_object(s, [Vector{Int}(GAPWrap.ExponentsOfPcElement(fullpcgs, x))
                    for x in GAP.Globals.InducedPcgsWrtHomePcgs(X)::GapObj], :gens)
  end
end

function load_object(s::DeserializerState, ::Type{SubFPGroup}, F::FPGroup)
  @assert haskey(s, :gens)
  Ffam = GAPWrap.FamilyObj(F)
  elfam = GAPWrap.ElementsFamily(Ffam)
  freegroup = GAP.getbangproperty(elfam, :freeGroup)::GapObj
  freefam = GAPWrap.FamilyObj(freegroup)
  elfreefam = GAPWrap.ElementsFamily(freefam)
  load_node(s) do d
    generators = load_object(s, Vector{Vector{Int}}, :gens)
    gens = [GAPWrap.ObjByExtRep(elfreefam, GapObj(x, true)) for x in generators]
    Ggens = [Oscar.group_element(F, GAPWrap.ElementOfFpGroup(elfam, x)) for x in gens]
    return sub(F, Ggens)[1]
  end
end

function load_object(s::DeserializerState, ::Type{SubPcGroup}, parent_group::PcGroup)
  X = GapObj(parent_group)
  elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
  fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)::GapObj
  load_node(s) do d
    generators = load_object(s, Vector{Vector{Int}}, :gens)
    Ggens = [Oscar.group_element(parent_group, GAP.Globals.PcElementByExponentsNC(fullpcgs, GapObj(x, true))::GapObj)
             for x in generators]
    return sub(parent_group, Ggens)[1]
  end
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

# homomorphisms
@register_serialization_type FinGenAbGroupHom uses_id

type_and_params(X::FinGenAbGroupHom) = TypeAndParams(
  FinGenAbGroupHom,
  :domain => domain(X),
  :codomain => codomain(X)
)

function save_object(s::SerializerState, h::FinGenAbGroupHom)
  save_object(s, matrix(h))
end

function load_object(s::DeserializerState, ::Type{FinGenAbGroupHom}, params::Dict)
  map_matrix = load_object(s, Matrix{ZZRingElem})
  return hom(params[:domain], params[:codomain], matrix(ZZ, map_matrix))
end

##############################################################################
# MatGroup

@register_serialization_type MatGroup uses_id

type_and_params(G::MatGroup) = TypeAndParams(MatGroup,
                                             :base_ring => base_ring(G),
                                             :degree => degree(G))

function save_object(s::SerializerState, G::MatGroup)
  save_data_dict(s) do
    save_object(s, matrix.(gens(G)), :gens)

    if isdefined(G, :descr)
      save_object(s, G.descr, :descr)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:MatGroup}, params::Dict)
  R = params[:base_ring]
  d = params[:degree]
  generators = load_object(s, Vector{dense_matrix_type(R)}, matrix_space(R, d, d), :gens)
  G = matrix_group(R, d, generators; check=false)

  if haskey(s, :descr)
    G.descr = load_object(s, Symbol, :descr)
  end
  return G
end

@register_serialization_type MatGroupElem

save_object(s::SerializerState, g::MatGroupElem) = save_object(s, matrix(g))

function load_object(s::DeserializerState, ::Type{<:MatGroupElem}, G::MatGroup)
  R = base_ring(G)
  d = degree(G)
  return G(matrix(R, load_object(s, dense_matrix_type(R), matrix_space(R, d, d))); check = false)
end
