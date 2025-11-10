push!(upgrade_scripts_set, UpgradeScript(
  v"1.6.0-1",
  function upgrade_1_6_0_1(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    # recurse upgrade on containers
    upgrade_containers(upgrade_1_6_0_1, s, dict)

    # Upgrades
    if dict[:_type] isa AbstractDict && get(dict[:_type], :name, nothing) == "RationalFunctionField"
      if dict[:data][:symbols] isa String
        dict[:data][:symbol] = dict[:data][:symbols]
        delete!(dict[:data], :symbols)
      end
    end
    return dict
  end
))
