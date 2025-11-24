push!(upgrade_scripts_set, UpgradeScript(
  v"1.3.0",
  function upgrade_1_3_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    upgraded_dict = dict
    if haskey(dict, :_type) && dict[:_type] == "FqField"
      if dict[:data] isa AbstractDict
        if !(haskey(dict[:data], :def_pol))
          upgraded_dict[:data][:def_pol] = copy(dict[:data])
        end
      end
    elseif haskey(dict, :data) && dict[:data] isa AbstractDict
      upgraded_dict[:data] = upgrade_1_3_0(s, dict[:data])
    end
    return upgraded_dict
  end
))
