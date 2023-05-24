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
    free = GAP.getbangproperty(elfam, :freeGroup)::GapObj
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(free))
    names = GAP.getbangproperty(elfam, :names)::GapObj
    if is_full_fp_group(G)
      rels = map(syllables, Oscar.relators(G))
      rels = [vcat([[x[1], x[2]] for x in l]...) for l in rels]
      return Dict(
          :symbols => save_type_dispatch(s, Vector{Symbol}(names)),
          :relators => save_type_dispatch(s, rels),
      )
    else
      # We have no Oscar object corresponding to the full f.p. group.
      whole = GAP.getbangproperty(GAPWrap.FamilyObj(G.X), :wholeGroup)::GapObj
      rels = Vector{Vector{Int}}(GAP.Globals.List(GAPWrap.RelatorsOfFpGroup(whole), GAPWrap.ExtRepOfObj))
      generators = Vector{Vector{Int}}(GAP.Globals.List(GAPWrap.GeneratorsOfGroup(G.X), GAPWrap.ExtRepOfObj))
      return Dict(
          :symbols => save_type_dispatch(s, Vector{Symbol}(names)),
          :relators => save_type_dispatch(s, rels),
          :generators => save_type_dispatch(s, generators),
      )
    end
end

function load_internal(s::DeserializerState, ::Type{FPGroup}, dict::Dict)
    names = load_type_dispatch(s, Vector{Symbol}, dict[:symbols])
    F = free_group(names)
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(F.X))
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

