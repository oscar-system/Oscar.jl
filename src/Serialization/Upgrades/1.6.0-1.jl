push!(upgrade_scripts_set, UpgradeScript(
  v"1.6.0-1",
  function upgrade_1_6_0_1(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    # recurse upgrade on containers
    upgrade_recursive(upgrade_1_6_0_1, s, dict)

    if dict[:_type] isa AbstractDict && haskey(dict[:_type], :name)
      type_name = dict[:_type][:name]
    else
      type_name = dict[:_type]
    end

    # Upgrades
    if type_name == "RationalFunctionField"
      if dict[:data][:symbols] isa String
        dict[:data][:symbol] = dict[:data][:symbols]
        delete!(dict[:data], :symbols)
      end
    end
    return dict
  end
))
