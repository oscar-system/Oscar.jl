import Oscar.Serialization: upgrade_scripts_set, UpgradeScript, UpgradeState, upgrade_recursive

push!(upgrade_scripts_set, UpgradeScript(
  v"1.8.0+2",
  function upgrade_1_8_0_2_phylo(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    upgrade_recursive(upgrade_1_8_0_2_phylo, s, dict)

    if dict[:_type] isa AbstractDict && haskey(dict[:_type], :name)
      type_name = dict[:_type][:name]
    else
      type_name = dict[:_type]
    end

    if type_name == "PhylogeneticModel"
      p = dict[:_type][:params]
      if p isa AbstractDict && haskey(p, :graph_type)
        p[:graph]             = p[:graph_params]
        p[:transition_matrix] = p[:transition_matrix_params]
        p[:root_distribution] = p[:root_distribution_params]
        p[:model_parameter]   = p[:model_parameter_name_type]
        for k in (:graph_type, :graph_params,
                  :transition_matrix_entry_type, :transition_matrix_params,
                  :root_distribution_entry_type, :root_distribution_params,
                  :model_parameter_name_type)
          delete!(p, k)
        end
      elseif p isa AbstractDict && haskey(p, :model_parameter_name)
        p[:model_parameter] = p[:model_parameter_name]
        delete!(p, :model_parameter_name)
      end
      if haskey(dict[:data], :model_parameter_name)
        dict[:data][:model_parameter] = dict[:data][:model_parameter_name]
        delete!(dict[:data], :model_parameter_name)
      end
    elseif type_name == "GroupBasedPhylogeneticModel"
      p = dict[:_type][:params]
      if p isa AbstractDict && haskey(p, :model_parameter_name_type)
        p[:model_parameter] = p[:model_parameter_name_type]
        delete!(p, :model_parameter_name_type)
      elseif p isa AbstractDict && haskey(p, :model_parameter_name)
        p[:model_parameter] = p[:model_parameter_name]
        delete!(p, :model_parameter_name)
      end
      if haskey(dict[:data], :varnames_group_based)
        dict[:data][:model_parameter] = dict[:data][:varnames_group_based]
        delete!(dict[:data], :varnames_group_based)
      end
    end
    return dict
  end
))
