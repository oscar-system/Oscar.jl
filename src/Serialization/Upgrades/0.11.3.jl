################################################################################
# Upgrade Summary
# When `is_basic_serialization_type(T) == true` the type `Vector{T}` the entries
# are no longer serialized as dicts, The type is encoded on the Vector at the same level as
# the entry data. The entries no longer being dictionaries prevents the use of backrefs.


push!(upgrade_scripts_set, UpgradeScript(
  v"0.11.3", # version this script upgrades to
   function upgrade_0_11_3(s::UpgradeState, dict::Dict)
    # moves down tree to point where type exists in dict
    # since we are only doing updates based on certain types
    # no :type key implies the dict is data
    if !haskey(dict, :type)
      return upgrade_data(upgrade_0_11_3, s, dict)
    end

    if dict[:type] == string(backref_sym)
      backrefed_object = s.id_to_dict[Symbol(dict[:id])]

      # if the backref points to a string, just use that string
      # instead of the backref
      backrefed_object isa String && return backrefed_object

      return dict
    end

    # this handles types like QQField that have no associated data
    if !haskey(dict, :data)
      # adds object to instance in case it is backrefed
      if haskey(dict, :id)
        s.id_to_dict[Symbol(dict[:id])] = dict
      end

      return dict
    end

    # apply upgrade to vector entries
    # we only do this here since upgrade_0_11_13 is only defined for dicts
    # we could make this its own function outside of this script
    if dict[:type] == "Vector"
      upgraded_vector = []
      entry_type = nothing
      for entry in dict[:data][:vector]
        if entry isa String
          result = entry
        else
          result = upgrade_0_11_3(s, entry)

          # store values in state that are vector entries that aren't backrefs
          if entry[:type] != string(backref_sym)
            if haskey(entry, :id)
              s.id_to_dict[Symbol(entry[:id])] = result
            end
            @assert entry_type === nothing || entry_type == entry[:type]
            entry_type = entry[:type]
          end
        end
        push!(upgraded_vector, result)
      end

      # if the objects in the vector are not basic, they are represented
      # here by dicts, and we do nothing else; if the objects are basic,
      # they are represented by strings, and we'll add the `entry_type`
      # to the vector data below...
      if upgraded_vector[1] isa Dict
        upgraded_dict = Dict(
          :type => "Vector",
          :id => dict[:id],
          :data => Dict(
            :vector => upgraded_vector
          )
        )
        # add to state
        s.id_to_dict[Symbol(dict[:id])] = upgraded_dict
        return upgraded_dict
      end

      upgraded_dict = Dict(
        :type => "Vector",
        :id => dict[:id],
        :data => Dict(
          :vector => upgraded_vector,
          :entry_type => entry_type
        )
      )
      s.id_to_dict[Symbol(dict[:id])] = upgraded_dict
      return upgraded_dict
    end

    dict_type = dict[:type]
    if contains(dict_type, "MPolyRingElem")
      dict_type = "MPolyRingElem"
    elseif contains(dict_type, "MPolyRing")
      dict_type = "MPolyRing"
    elseif contains(dict_type, "PolyRingElem")
      dict_type = "PolyRingElem"
    elseif contains(dict_type, "PolyRing")
      dict_type = "PolyRing"
    end

    # Upgrades basic types
     if dict_type in ["ZZRingElem", "String", "Symbol", "Bool", "QQFieldElem", "zzModRing"] || contains(dict_type, "Int")
dict_type
      if dict_type == "QQFieldElem"
        num = dict[:data][:num][:data]
        den = dict[:data][:den][:data]
        updated_fmpq = "$num//$den"
        #add to state
        s.id_to_dict[Symbol(dict[:id])] = updated_fmpq
        
        return updated_fmpq
      else
        updated_basic_value = dict[:data]
        #add to state
        if haskey(dict, :id)
          s.id_to_dict[Symbol(dict[:id])] = updated_basic_value
        end
        
        return updated_basic_value
      end
    end

    upgraded_data = upgrade_0_11_3(s, dict[:data])
    upgraded_dict = Dict(
      :type => dict_type,
      :data => upgraded_data,
      :id => dict[:id]
    )

    # adds updated object in case it is referenced somewhere
    if haskey(dict, :id)
      s.id_to_dict[Symbol(dict[:id])] = upgraded_dict
    end

    return upgraded_dict
  end
))
