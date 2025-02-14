function upgrade_NamedTuple(d::Dict)
  println(json(d, 2))
  for (i, entry) in enumerate(d[:data])
    println(entry)
    if entry isa Dict
      func = Symbol(string("upgrade_", entry[:_type][:name]))
      upgraded_entry = @eval $func($dict)
    end
  end
end

push!(upgrade_scripts_set, UpgradeScript(
  v"1.3.0",
  function upgrade_1_3_0(s::UpgradeState, dict::Dict)
    if haskey(dict, :_type) 
      func = Symbol(string("upgrade_", dict[:_type][:name]))
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
