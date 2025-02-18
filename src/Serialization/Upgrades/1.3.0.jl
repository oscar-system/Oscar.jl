function upgrade_NamedTuple(d::Dict)
  for (i, entry) in enumerate(d[:data])
    if entry isa Dict
      func = Symbol(string("upgrade_", entry[:_type][:name]))
      upgraded_entry = @eval $func($dict)
    end
  end
end

function upgrade_PolyRingElem(d::Dict)
  return d
end

function upgrade_PolyRing(d::Dict)
  return d
end

function upgrade_fqPolyRepField(d::Dict)
  upgraded_dict = d
  upgraded_dict[:_type] = Dict(
    :name => d[:_type],
    :params => d[:data][:def_pol][:_type][:params]
  )
  upgraded_dict[:data] = d[:data][:def_pol][:data]
  return upgraded_dict
end

push!(upgrade_scripts_set, UpgradeScript(
  v"1.3.0",
  function upgrade_1_3_0(s::UpgradeState, dict::Dict)
    if haskey(dict, :_type)
      if dict[:_type] isa Dict && haskey(dict[:_type], :name)
        type_name = dict[:_type][:name]
      else
        type_name = dict[:_type]
      end
      func = Symbol(string("upgrade_", type_name))
      upgraded_dict = @eval $func($dict)
        
    elseif haskey(dict, :data) && dict[:data] isa Dict
      upgraded_dict[:data] = upgrade_1_3_0(s, dict[:data])
    end

    if haskey(dict, :_refs)
      upgraded_refs = Dict()
      for (k, v) in dict[:_refs]
        upgraded_refs[k] = upgrade_1_3_0(s, v)
      end
      upgraded_dict[:_refs] = upgraded_refs
    end

    return upgraded_dict
  end
))
