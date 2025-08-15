push!(upgrade_scripts_set, UpgradeScript(
  v"1.5.0",
  function upgrade_1_5_0(s::UpgradeState, dict::Dict)
    upgraded_dict = dict

    if haskey(dict, :_type) && dict[:_type] == "Matroid"
      # Now we store basis, before not. But basis is stored in :lattice
      #if !haskey(dict[:data], :basis)
      #  upgraded_dict[:data][:basis] = dict[:data][:lattice][:data][:basis]
      #end
      error("TODO")
    end

    return upgraded_dict
  end
))
