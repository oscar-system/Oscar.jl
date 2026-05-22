################################################################################
# Upgrade Summary
# - `ZZLatWithIsom`: we now only serialized the ambient `QuadSpaceWithIsom` and
# a basis matrix of the underlying lattice; the upgrade recovers the basis matrix
# in the old tree.

push!(upgrade_scripts_set, UpgradeScript(
  v"1.1.0",
  function upgrade_1_1_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})

    upgraded_dict = dict

    if haskey(dict, :_type) && dict[:_type] == "ZZLatWithIsom"
      # Now we store basis, before not. But basis is stored in :lattice
      if !haskey(dict[:data], :basis)
        upgraded_dict[:data][:basis] = dict[:data][:lattice][:data][:basis]
      end

      # In the past, we saved only finite order isometries, but the order
      # was not typed. Now we need to reference the type otherwise we cannot
      # read old files.
      if !(dict[:data][:ambient_space][:data][:order] isa AbstractDict)
        n = Base.parse(Int, dict[:data][:ambient_space][:data][:order])
        upgraded_dict[:data][:ambient_space][:data][:order] =
        Dict{Symbol, Any}(:_type => "Base.Int", :data => "$n")
      end
    end

    return upgraded_dict
  end
))
