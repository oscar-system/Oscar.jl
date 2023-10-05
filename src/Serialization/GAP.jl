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

function _describe_GAP_flags(flags::GapObj)
  return [string(GAPWrap.NameFunction(x))
          for x in GAPWrap.ELMS_LIST(GAP.Globals.FILTERS,
                     GAP.Globals.TRUES_FLAGS(flags))]
#TODO: omit implied filters
end

function _GAP_flags(filtnames::Vector)
  # We need `GAP.evalstr` because of function names involving brackets.
  l = [GAP.Globals.FLAGS_FILTER(GAP.evalstr(x))::GapObj
       for x in filtnames]
  fl = l[1]
  for f in l[2:end]
    fl = GAP.Globals.AND_FLAGS(fl, f)
  end
  return fl
end

#############################################################################
#
# the Oscar (de)serialization methods that delegate to GAP's method selection
#
@register_serialization_type GAP.GapObj uses_id

function save_object(s::SerializerState, X::GapObj)
  GAP.Globals.SerializeInOscar(X, s)
end

function load_object(s::DeserializerState, T::Type{GapObj}, d::Dict{Symbol, Any})
  @req haskey(d, :GapType) "cannot deserialize GapObj without key :GapType"
  GAP_T = GapObj(d[:GapType])
  return GAP.Globals.DeserializeInOscar(GAPWrap.ValueGlobal(GAP_T), s, T, d)
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
#   full free group or subgroup of it
install_GAP_serialization(:IsFreeGroup, 
  function(X::GapObj, s::SerializerState)
    # save the defining data
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
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
      if GAP.Globals.Length(Xnames)::GAP.Obj == GAP.Globals.infinity
        # store the initial names and the prefix
        prefix = GAP.getbangposition(Xnames, 1)::GapObj
        save_object(s, string(prefix), :nameprefix)
        names = Vector{String}(GAP.getbangposition(Xnames, 2)::GapObj)
      else
        # store the names
        names = Vector{String}(Xnames)
      end
      save_data_array(s, :names) do
        for i in names
          save_object(s, i)
        end
      end
      if !(GAP.Globals.HasIsWholeFamily(X) &&
           GAPWrap.IsWholeFamily(X))
        # store generators
        save_data_array(s, :gens) do
          for x in GAPWrap.GeneratorsOfGroup(X)
            v = Vector{Int}(GAPWrap.ExtRepOfObj(x))
            save_data_array(s) do
              for i in v
                save_object(s, i)
              end
            end
          end
        end
      end
    end
  end)

install_GAP_deserialization(:IsFreeGroup,
  function(filt::GapObj, s::DeserializerState, T, d::Dict)
    # We have to create a new free group.
    wfilt = GAP.evalstr(d[:wfilt])
    if haskey(d, :nameprefix)
      # infinite rank
#TODO: support these groups in Oscar?
      prefix = string(d[:nameprefix])
      init = GapObj(d[:names], true)
      G = GAP.Globals.FreeGroup(wfilt, GAP.Globals.infinity, prefix, init)
    else
      names = GapObj(d[:names], true)
      G = GAP.Globals.FreeGroup(wfilt, names)
    end

    if haskey(d, :gens)
      generators = load_object(s, Vector, d[:gens], (Vector{Int}, Int))
      F = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(G))
      Ggens = [GAPWrap.ObjByExtRep(F, GapObj(x, true)) for x in generators]
      G = GAP.Globals.SubgroupNC(G, GapObj(Ggens))
    end

    return G
  end)

# # - IsElementOfFreeGroup
# install_GAP_serialization(:IsElementOfFreeGroup, 
#   function(X::GapObj, s::SerializerState)
#     cfam = GAPWrap.FamilyObj(X)::GapObj
#     efam = GAPWrap.ElementsFamily(cfam)::GapObj
#     G = GAP.getbangproperty(efam, :freeGroup)
#     save_data_dict(s) do
#       save_object(s, G, :freeGroup)
#     end
#   end)
# 
# install_GAP_deserialization(:IsElementOfFreeGroup,
#   function(filt::GapObj, s::DeserializerState, T, d::Dict)
# 
#   end)


# - `IsSubgroupFpGroup`:
#   full f.p. group or subgroup of it
install_GAP_serialization(:IsSubgroupFpGroup, 
  function(X::GapObj, s::SerializerState)
    # save the defining data
    elfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(X))
    save_data_dict(s) do
      save_object(s, "IsSubgroupFpGroup", :GapType)

      # underlying free group
      freegroup = GAP.getbangproperty(elfam, :freeGroup)::GapObj
      save_object(s, freegroup, :freeGroup)

      # relators
      relators = GAP.getbangproperty(elfam, :relators)::GapObj
      save_data_array(s, :relators) do
        for x in relators
          v = Vector{Int}(GAPWrap.ExtRepOfObj(x))
          save_data_array(s) do
            for i in v
              save_object(s, i)
            end
          end
        end
      end

      if !(GAP.Globals.HasIsWholeFamily(X) &&
           GAPWrap.IsWholeFamily(X))
        # generators
        save_data_array(s, :gens) do
          for x in GAPWrap.GeneratorsOfGroup(X)
            v = Vector{Int}(GAPWrap.ExtRepOfObj(x))
            save_data_array(s) do
              for i in v
                save_object(s, i)
              end
            end
          end
        end
      end
    end
  end)

install_GAP_deserialization(:IsSubgroupFpGroup,
  function(filt::GapObj, s::DeserializerState, T, d::Dict)
    # Create the free group.
    F = load_object(s, GapObj, d[:freeGroup])
    freefam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(F))
    relators = load_object(s, Vector, d[:relators], (Vector{Int}, Int))
    rels = [GAPWrap.ObjByExtRep(freefam, GapObj(x, true)) for x in relators]

    # Create the f.p. group.
    G = F/GapObj(rels)
    if haskey(d, :gens)
      # Create the subgroup.
      generators = load_object(s, Vector, d[:gens], (Vector{Int}, Int))
      gens = [GAPWrap.ObjByExtRep(freefam, GapObj(x, true)) for x in generators]
      fpfam = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(G))
      Ggens = [GAPWrap.ElementOfFpGroup(fpfam, x) for x in gens]
      G = GAP.Globals.SubgroupNC(G, GapObj(Ggens))
    end

    return G
  end)
