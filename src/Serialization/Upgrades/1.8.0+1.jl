# Upgrades are to allow the serialization of more subtypes of phylogenetic models
# since the only models that were serializable were the ones with
# the following set of parameters we hardcode them here in the upgrade
push!(upgrade_scripts_set, UpgradeScript(
  v"1.8.0+1",
  function upgrade_1_8_0_1(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    # recurse upgrade on containers
    upgrade_recursive(upgrade_1_8_0_1, s, dict)

    if dict[:_type] isa AbstractDict && haskey(dict[:_type], :name)
      type_name = dict[:_type][:name]
    else
      type_name = dict[:_type]
    end

    # Upgrades
    if type_name == "PhylogeneticModel"
      dict[:_type][:params][:model_parameter_name_type] = "Symbol"
      dict[:_type][:params][:transition_matrix_entry_type] = "Symbol"
      dict[:_type][:params][:transition_matrix_params] = Dict(:name => "Matrix", :params => "Symbol")
      dict[:_type][:params][:root_distribution_entry_type] = "QQFieldElem"
      dict[:_type][:params][:root_distribution_params] = Dict(:name => "Vector",
                                                              :params => Dict(:name => "QQFieldElem",
                                                                              :params => Dict(:_type => "QQField")))
    elseif type_name == "GroupBasedPhylogeneticModel"
      dict[:_type][:params][:model_parameter_name_type] = "Symbol"
    end
    return dict
  end
))
    
