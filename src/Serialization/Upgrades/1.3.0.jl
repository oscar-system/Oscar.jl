push!(upgrade_scripts_set, UpgradeScript(
  v"1.3.0",
  function upgrade_1_3_0(s::UpgradeState, dict::Dict)
    upgraded_dict = dict
    if haskey(dict, :_type) && dict[:_type] == "FqField"
      if dict[:data] isa Dict
        if !(haskey(dict[:data], :def_pol)) 
          upgraded_dict[:data][:def_pol] = dict[:data]
        end
      end
    end
  end
))
