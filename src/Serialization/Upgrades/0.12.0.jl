################################################################################
# Upgrade Summary
# This upgrade restricts the use of ids to types that can not be algorithmic
# determined to be equal. Hence any if two instances of a type be decided to
# be equal given their encodings no longer have ids, hence they can no longer
# be back referenced.

push!(upgrade_scripts_set, UpgradeScript(
  v"0.12.0",
  function upgrade_0_12_0(s::UpgradeState, dict::Dict)
    # moves down tree to point where type exists in dict
    # since we are only doing updates based on certain types
    # no :type key implies the dict is data
    if !haskey(dict, :type)
      return upgrade_data(upgrade_0_12_0, s, dict)
    end

    if dict[:type] == string(backref_sym)
      backrefed_dict = s.id_to_dict[Symbol(dict[:id])]
      T = reverse_type_map[backrefed_dict[:type]]
      
      if serialize_with_id(T)
        #backref are still ok
        return dict
      else
        # backref cannot be used, since this type no longer
        # has an id
        return backrefed_dict
      end
    end
    
    upgraded_dict = dict

    if contains(upgraded_dict[:type], "MPolyRingElem")
      upgraded_dict[:type] = "MPolyRingElem"
    elseif contains(upgraded_dict[:type], "MPolyRing")
      upgraded_dict[:type] = "MPolyRing"
    elseif contains(upgraded_dict[:type], "PolyRingElem")
      upgraded_dict[:type] = "PolyRingElem"
    elseif contains(upgraded_dict[:type], "PolyRing")
      upgraded_dict[:type] = "PolyRing"
    end


    if upgraded_dict[:type] == "Nemo.fpFieldElem"
      upgraded_dict[:type] = "fpFieldElem"
    end

    # handle data upgrade (recurse on sub tree)
    if haskey(dict, :data) && dict[:data] isa Dict
      upgraded_dict[:data] = upgrade_0_12_0(s, dict[:data])
    end

    T = reverse_type_map[upgraded_dict[:type]]
    
    # remove ids according to update
    if serialize_with_id(T)
      s.id_to_dict[Symbol(dict[:id])] = upgraded_dict
    elseif haskey(dict, :id)
      id = pop!(upgraded_dict, :id)
      s.id_to_dict[Symbol(id)] = upgraded_dict
    end

    return upgraded_dict
  end
))
