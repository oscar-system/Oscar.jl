################################################################################
# Upgrade Summary
# 
# - `ZZLatWithIsom`: we now only serialized the ambient `QuadSpaceWithIsom` and
# a basis matrix of the underlying lattice; the upgrade recovers the basis matrix
# in the old tree.

push!(upgrade_scripts_set, UpgradeScript(
  v"1.1.0",
  function upgrade_1_1_0(s::UpgradeState, dict::Dict)

    upgraded_dict = dict

    if haskey(dict, :_type) && dict[:_type] == "ZZLatWithIsom"
      # Now we store basis, before not. But basis is stored in :lattice
      if !haskey(dict[:data], :basis)
        upgraded_dict[:data][:basis] = dict[:data][:lattice][:data][:basis]
      end
    end

    return upgraded_dict
  end
))
