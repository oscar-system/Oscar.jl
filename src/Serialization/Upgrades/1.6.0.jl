push!(upgrade_scripts_set, UpgradeScript(
  v"1.6.0",
  function upgrade_1_6_0(s::UpgradeState, dict::Dict)
    upgrade_containers(upgrade_1_6_0, s, dict)
    
    return dict
  end
))
