
# Cache families of GAP objects that have been serialized,
# in order to be able to deserialize the objects in the *same* GAP session,
# and to get equal objects.
# If deserialization takes place in a *different* GAP session then
# we have to create the families in question anew.

const _GAP_Families = Dict{String, Dict{String, GapObj}}()
const _GAP_session_id = Ref{String}("")

function __init_group_serialization()
    _GAP_session_id[] = string(objectid(GAP))
    _GAP_Families[_GAP_session_id[]] = Dict{UInt64, GapObj}()
end


################################################################################
# PermGroup

@registerSerializationType(PermGroup)

function save_internal(s::SerializerState, G::PermGroup)
    n = degree(G)
    generators = [Vector{Int}(GAP.Globals.ListPerm(x.X)) for x in gens(G)]

    return Dict(
        :degree => save_type_dispatch(s, n),
        :gens => save_type_dispatch(s, generators),
    )
end

function load_internal(s::DeserializerState, ::Type{PermGroup}, dict::Dict)
    n = load_type_dispatch(s, Int, dict[:degree])
    generators = load_type_dispatch(s, Vector{Vector{Int}}, dict[:gens])

    return permutation_group(n, [perm(x) for x in generators])
end


################################################################################
# PermGroupElem

@registerSerializationType(PermGroupElem)

function save_internal(s::SerializerState, p::PermGroupElem)
    return Dict(
        :parent => save_type_dispatch(s, parent(p)),
        :imgs => save_type_dispatch(s, Vector{Int}(GAP.Globals.ListPerm(p.X))),
    )
end

function load_internal(s::DeserializerState, T::Type{PermGroupElem}, dict::Dict)
    parent_group = load_unknown_type(s, dict[:parent])
    return load_internal_with_parent(s, T, dict, parent_group)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{PermGroupElem},
                                   dict::Dict,
                                   parent_group::PermGroup)
    imgs = load_type_dispatch(s, Vector{Int}, dict[:imgs])

    return perm(parent_group, imgs)
end


################################################################################
# FPGroup

@registerSerializationType(FPGroup)

function save_internal(s::SerializerState, G::FPGroup)
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(G.X))
    id = string(objectid(elfam))
    D = _GAP_Families[_GAP_session_id[]]
    if !haskey(D, id)
      D[id] = elfam
    end

    free = GAP.getbangproperty(elfam, :freeGroup)::GapObj
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(free))
    names = GAP.getbangproperty(elfam, :names)::GapObj

    res = Dict(
            :session => _GAP_session_id[],
            :family => string(id),
            :symbols => save_type_dispatch(s, Vector{Symbol}(names)),
          )

    if !is_full_fp_group(G)
      # We have no Oscar object corresponding to the full f.p. group.
      whole = GAP.getbangproperty(GAPWrap.FamilyObj(G.X), :wholeGroup)::GapObj
      generators = Vector{Vector{Int}}(GAP.Globals.List(GAPWrap.GeneratorsOfGroup(G.X), GAPWrap.ExtRepOfObj))
      res[:generators] = save_type_dispatch(s, generators)
      if !GAP.Globals.IsFreeGroup(G.X)
        rels = Vector{Vector{Int}}(GAP.Globals.List(GAPWrap.RelatorsOfFpGroup(whole), GAPWrap.ExtRepOfObj))
        res[:relators] = save_type_dispatch(s, rels)
      end
    elseif !GAP.Globals.IsFreeGroup(G.X)
      rels = map(syllables, Oscar.relators(G))
      rels = [vcat([[x[1], x[2]] for x in l]...) for l in rels]
      res[:relators] = save_type_dispatch(s, rels)
    end

    return res
end

function load_internal(s::DeserializerState, ::Type{FPGroup}, dict::Dict)
    if !haskey(_GAP_Families, dict[:session])
      # This is the first time we deserialize an object from that session id.
      _GAP_Families[dict[:session]] = Dict{UInt64, GapObj}()
    end
    D = _GAP_Families[dict[:session]]

    if !haskey(dict, :relators)
      # (subgroup of) a free group,
      if haskey(D, dict[:family])
        # Use the stored elements family.
        elfam = D[dict[:family]]
        G = FPGroup(GAP.getbangproperty(elfam, :freeGroup))
      else
        # Create the family anew, and store it (under the *old* key).
        names = load_type_dispatch(s, Vector{Symbol}, dict[:symbols])
        G = free_group(names)
        elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(G.X))
        D[dict[:family]] = elfam
      end
    else
      # (subgroup of) quotient of a free group
      if haskey(D, dict[:family])
        # Use the stored elements family.
        elfam = D[dict[:family]]
        G = FPGroup(GAP.getbangproperty(GAPWrap.CollectionsFamily(elfam), :wholeGroup))
      else
        # Create the family anew, and store it (under the *old* key).
        names = load_type_dispatch(s, Vector{Symbol}, dict[:symbols])
        F = free_group(names)
        rels = load_type_dispatch(s, Vector{Vector{Int}}, dict[:relators])
        rels_pairs = Vector{Pair{Int, Int}}[]
        for l in rels
          rel = Pair{Int, Int}[]
          for i in 1:2:(length(l)-1)
            push!(rel, Pair{Int, Int}(l[i], l[i+1]))
          end
          push!(rels_pairs, rel)
        end
        G = quo(F, [F(l) for l in rels_pairs])[1]
        elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(G.X))
        D[dict[:family]] = elfam
      end
    end

    if haskey(dict, :generators)
      generators = load_type_dispatch(s, Vector{Vector{Int}}, dict[:generators])
      gens_pairs = Vector{Pair{Int, Int}}[]
      for l in generators
        gen = Pair{Int, Int}[]
        for i in 1:2:(length(l)-1)
          push!(gen, Pair{Int, Int}(l[i], l[i+1]))
        end
        push!(gens_pairs, gen)
      end

      G = sub(G, [G(l) for l in gens_pairs])[1]
    end
    return G
end

################################################################################
# FPGroupElem

@registerSerializationType(FPGroupElem)

function save_internal(s::SerializerState, g::FPGroupElem)
    return Dict(
        :parent => save_type_dispatch(s, parent(g)),
        :extrep => save_type_dispatch(s, vcat([[x[1], x[2]] for x in syllables(g)]...))
    )
end

function load_internal(s::DeserializerState, T::Type{FPGroupElem}, dict::Dict)
    parent_group = load_unknown_type(s, dict[:parent])
    return load_internal_with_parent(s, T, dict, parent_group)
end

function load_internal_with_parent(s::DeserializerState,
                                   ::Type{FPGroupElem},
                                   dict::Dict,
                                   parent_group::FPGroup)
    l = load_type_dispatch(s, Vector{Int}, dict[:extrep])
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(parent_group.X))
    free = GAP.getbangproperty(elfam, :freeGroup)::GapObj
    freefam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(free))
    v = GapObj(l)
    w = GAPWrap.ObjByExtRep(freefam, v)
    return group_element(parent_group, GAP.Globals.ElementOfFpGroup(elfam, w))
end

