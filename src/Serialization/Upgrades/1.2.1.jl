push!(upgrade_scripts_set, UpgradeScript(
  v"1.2.1",
  function upgrade_1_2_1(s::UpgradeState, dict::Dict)
    upgraded_dict = dict
    if haskey(dict, :_type) && dict[:_type] == "FqField"
      if dict[:data] isa Dict
        if !(haskey(dict[:data], :def_pol))
          upgraded_dict[:data][:def_pol] = copy(dict[:data])
        end
      end
    end
    return upgraded_dict
  end
))
