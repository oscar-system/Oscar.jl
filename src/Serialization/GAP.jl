##############################################################################
#
# General ideas of the (de)serialization of GAP objects:
#
# - We define (de)serialization methods for GAP objects depending on
#   GAP filters, and store them via calls to `install_GAP_serialization`
#   and `install_GAP_deserialization`, at compile time.
#
# - We install these methods at runtime (when the OscarInterface package
#   of GAP gets loaded) as methods for the GAP operations `SerializeInOscar`
#   and `DeserializeInOscar`.

#############################################################################
#
# utilities for installing methods depending on GAP filters
#
const _GAP_serializations = Pair{Symbol, Function}[]
const _GAP_deserializations = Pair{Symbol, Function}[]

function install_GAP_serialization(filtsymbol::Symbol, meth::Function)
  push!(_GAP_serializations, filtsymbol => meth)
  return
end

function install_GAP_deserialization(filtsymbol::Symbol, meth::Function)
  push!(_GAP_deserializations, filtsymbol => meth)
  return
end

#############################################################################
#
# the Oscar (de)serialization methods that delegate to GAP's method selection
#
@register_serialization_type GapObj uses_id

function save_object(s::SerializerState, X::GapObj)
  GAP.Globals.SerializeInOscar(X, s)
end

function load_object(s::DeserializerState, T::Type{GapObj})
  load_node(s) do d
    @req haskey(s, :GapType) "cannot deserialize GapObj without key :GapType"
    GAP_T = load_node(s, :GapType) do gap_type_data
      return GapObj(gap_type_data)
    end
    return GAP.Globals.DeserializeInOscar(GAPWrap.ValueGlobal(GAP_T), s, T)
  end
end

#############################################################################
#
# the individual methods

# - `IsObject`:
#   generic method, throw an exception
#
install_GAP_serialization(:IsObject,
  function(X::GapObj, s::SerializerState)
    error("serialization of GAP object $X is not yet supported")
  end)

# - `IsFamily`:
#   These objects are hidden inside GAP,
#   and typically their construction is encapsulated in functions which
#   create some domain.
#   In this situation, one cannot create the family object separately
#   and then later call a function that takes this object and creates
#   a domain based on it.
#   Thus (de)serializing family objects does not help us.
#
install_GAP_serialization(:IsFamily,
  function(X::GapObj, s::SerializerState)
    error("serialization of GAP family object $X is deliberately not supported")
  end)

# - `IsFreeGroup`:
#   full free group or subgroup of it,
#   distinguished by presence of `:freeGroup` and `:gens` in case of a subgroup
install_GAP_serialization(:IsFreeGroup,
  function(X::GapObj, s::SerializerState)
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
    if GAP.Globals.HasIsWholeFamily(X) && GAPWrap.IsWholeFamily(X)
      # full free group: Save the defining data.
      save_data_dict(s) do
        save_object(s, "IsFreeGroup", :GapType)
        # the internal representation of elements
        if GAP.Globals.IsSyllableWordsFamily(elfam)::Bool
          wfilt = "IsSyllableWordsFamily"
        elseif GAP.Globals.IsLetterWordsFamily(elfam)::Bool
          wfilt = "IsLetterWordsFamily"
        elseif GAP.Globals.IsWLetterWordsFamily(elfam)::Bool
          wfilt = "IsWLetterWordsFamily"
        elseif GAP.Globals.IsBLetterWordsFamily(elfam)::Bool
          wfilt = "IsBLetterWordsFamily"
        else
          error("not supported internal representation")
        end
        save_object(s, wfilt, :wfilt)
        # rank and names of generators
        Xnames = GAP.getbangproperty(elfam, :names)::GapObj
        if !GAPWrap.IsFinite(Xnames)
          # store the initial names and the prefix
          prefix = GAP.getbangindex(Xnames, 1)::GapObj
          save_object(s, string(prefix), :nameprefix)
          names = Vector{String}(GAP.getbangindex(Xnames, 2)::GapObj)
        else
          # store the names
          names = Vector{String}(Xnames)
        end
        save_object(s, names, :names)
      end
    else
      # subgroup of a full free group: save the full group and generators
      save_data_dict(s) do
        F = GAP.getbangproperty(elfam, :freeGroup)::GapObj
        save_object(s, "IsFreeGroup", :GapType)
        save_typed_object(s, F, :freeGroup)
        # store generators
        save_object(s, [Vector{Int}(GAPWrap.ExtRepOfObj(x)) for x in GAPWrap.GeneratorsOfGroup(X)], :gens)
      end
    end
  end)

install_GAP_deserialization(
  :IsFreeGroup,
  function(filt::GapObj, s::DeserializerState, T)
    load_node(s) do d
      if haskey(s, :freeGroup) && haskey(s, :gens)
        # Deserialize the full free group.
        F = load_typed_object(s, :freeGroup)
        # Deserialize the generators.
        generators = load_object(s, Vector, (Vector{Int}, Int), :gens)
        fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(F))
        Ggens = [GAPWrap.ObjByExtRep(fam, GapObj(x, true)) for x in generators]
        # Create the subgroup.
        G = GAP.Globals.SubgroupNC(F, GapObj(Ggens))::GapObj
      else
        # Create a new full free group.
        wfilt = getproperty(GAP.Globals, load_object(s, Symbol, :wfilt))::GapObj
        if haskey(s, :nameprefix)
          # infinite rank
          prefix = load_node(s, :nameprefix) do nameprefix
            GapObj(nameprefix)
          end
          init = load_node(s, :names) do names
            GapObj(names; recursive = true)
          end
          G = GAP.Globals.FreeGroup(wfilt, GAP.Globals.infinity, prefix, init)::GapObj
        else
          init = load_node(s, :names) do names
            GapObj(names; recursive = true)
          end
          G = GAP.Globals.FreeGroup(wfilt, init)::GapObj
        end
      end
      return G
    end
  end)

# - `IsSubgroupFpGroup`:
#   full f.p. group or subgroup of it,
#   distinguished by presence of `:wholeGroup` and `:gens` in case of a subgroup
install_GAP_serialization(:IsSubgroupFpGroup,
  function(X::GapObj, s::SerializerState)
    Xfam = GAPWrap.FamilyObj(X)
    elfam = GAPWrap.ElementsFamily(Xfam)
    if GAP.Globals.HasIsWholeFamily(X) && GAPWrap.IsWholeFamily(X)
      # full f.p. group: Save the defining data.
      save_data_dict(s) do
        save_object(s, "IsSubgroupFpGroup", :GapType)
        # underlying free group
        freegroup = GAP.getbangproperty(elfam, :freeGroup)::GapObj
        save_typed_object(s, freegroup, :freeGroup)
        # relators
        relators = GAP.getbangproperty(elfam, :relators)::GapObj
        save_object(s, [Vector{Int}(GAPWrap.ExtRepOfObj(x)) for x in relators], :relators)
      end
    else
      # subgroup of a full f.p. group: save the full group and generators
      save_data_dict(s) do
        F = GAP.getbangproperty(Xfam, :wholeGroup)::GapObj
        save_object(s, "IsSubgroupFpGroup", :GapType)
        save_typed_object(s, F, :wholeGroup)
        # store generators
        save_object(s, [Vector{Int}(GAPWrap.ExtRepOfObj(x)) for x in GAPWrap.GeneratorsOfGroup(X)], :gens)
      end
    end
  end)

install_GAP_deserialization(
  :IsSubgroupFpGroup,
  function(filt::GapObj, s::DeserializerState, T)
    load_node(s) do d 
      if haskey(s, :wholeGroup) && haskey(s, :gens)
        # Deserialize the full f.p. group.
        F = load_typed_object(s, :wholeGroup)
        Ffam = GAPWrap.FamilyObj(F)
        elfam = GAPWrap.ElementsFamily(Ffam)
        freegroup = GAP.getbangproperty(elfam, :freeGroup)::GapObj
        freefam = GAPWrap.FamilyObj(freegroup)
        elfreefam = GAPWrap.ElementsFamily(freefam)
        # Deserialize the generators.
        generators = load_object(s, Vector, (Vector{Int}, Int), :gens)
        gens = [GAPWrap.ObjByExtRep(elfreefam, GapObj(x, true)) for x in generators]
        Ggens = [GAPWrap.ElementOfFpGroup(elfam, x) for x in gens]
        # Create the subgroup.
        G = GAP.Globals.SubgroupNC(F, GapObj(Ggens))::GapObj
      else
        # Create a new full f.p. group.
        F = load_typed_object(s, :freeGroup)
        relators = load_object(s, Vector,  (Vector{Int}, Int), :relators)
        elfreefam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(F))
        rels = [GAPWrap.ObjByExtRep(elfreefam, GapObj(x, true)) for x in relators]
        G = F/GapObj(rels)
      end
      return G
    end
  end)


# - `IsPcGroup`:
#   full pc group or subgroup of it
#   we do not support (de)serialization of the stored rws,
#   thus we need not (de)serialize its underlying free group etc.

# utility, turn an exponent vector `[a_1, a_2, ..., a_n]`
# into `[1, a_1, 2, a_2, ..., n, a_n]`
function _free_group_extrep_from_exponents(exps::Vector{Int})
  res = Int[]
  for i in 1:length(exps)
    if exps[i] != 0
      push!(res, i)
      push!(res, exps[i])
    end
  end
  return res
end

install_GAP_serialization(:IsPcGroup,
  function(X::GapObj, s::SerializerState)
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
    fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)::GapObj
    if fullpcgs === GAPWrap.Pcgs(X)
      # full pc group: Save the defining data.
      save_data_dict(s) do
        save_object(s, "IsPcGroup", :GapType)
        # relative orders
        relord = [GAP.Globals.RelativeOrderOfPcElement(fullpcgs, x)::GAP.Obj
                  for x in fullpcgs]
        save_object(s, relord, :relord)
        # power relators
        rels = Tuple{Int, Vector{Int}}[]
        for i in 1:length(relord)
          ne = fullpcgs[i]^relord[i]
          if ! GAPWrap.IsOne(ne)
            push!(rels, (i, _free_group_extrep_from_exponents(
              Vector{Int}(GAP.Globals.ExponentsOfPcElement(fullpcgs, ne)::GapObj))))
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
                Vector{Int}(GAP.Globals.ExponentsOfPcElement(fullpcgs, ne)::GapObj))))
            end
          end
        end
        save_object(s, rels, :comm_rels)
      end
    else
      # save full group and generators
      save_data_dict(s) do
        save_object(s, "IsPcGroup", :GapType)
        G = GAP.getbangproperty(fullpcgs, :GroupOfPcgs)::GapObj
        save_typed_object(s, G, :fullGroup)
        save_object(s, [Vector{Int}(GAP.Globals.ExponentsOfPcElement(fullpcgs, x)::GapObj)
                        for x in GAP.Globals.InducedPcgsWrtHomePcgs(X)::GapObj], :gens)
      end
    end
  end)

install_GAP_deserialization(
  :IsPcGroup,
  function(filt::GapObj, s::DeserializerState, T)
    if haskey(s, :relord)
      # full pc group
      relord = load_object(s, Vector, Int, :relord)
      F = GAP.Globals.FreeGroup(GAP.Globals.IsSyllableWordsFamily,
                                length(relord))::GapObj
      fam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(F))
      rws = GAP.Globals.SingleCollector(F, GapObj(relord))::GapObj
      for (i, elm) in load_object(s, Vector, (Tuple, [Int, [Vector, Int]]), :power_rels)
        GAP.Globals.SetPower(rws, i, GAPWrap.ObjByExtRep(fam, GapObj(elm)))
      end
      for (j, i, elm) in load_object(s, Vector, (Tuple, [Int, Int, [Vector, Int]]), :comm_rels)
        GAP.Globals.SetCommutator(rws, j, i, GAPWrap.ObjByExtRep(fam, GapObj(elm)))
      end
      G = GAP.Globals.GroupByRwsNC(rws)::GapObj
    else
      # Deserialize the full pc group.
      F = load_typed_object(s, :fullGroup)
      elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(F))
      fullpcgs = GAP.getbangproperty(elfam, :DefiningPcgs)::GapObj
      # Deserialize the generators.
      generators = load_object(s, Vector, (Vector{Int}, Int), :gens)
      Ggens = [GAP.Globals.PcElementByExponentsNC(fullpcgs, GapObj(x, true))::GapObj
               for x in generators]
      # Create the subgroup.
      G = GAP.Globals.SubgroupNC(F, GapObj(Ggens))::GapObj
    end
    return G
  end)
