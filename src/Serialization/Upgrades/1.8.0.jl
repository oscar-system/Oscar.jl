push!(upgrade_scripts_set, UpgradeScript(
  v"1.8.0",
  function upgrade_1_8_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    # recurse upgrade on containers
    upgrade_recursive(upgrade_1_8_0, s, dict)
    if dict[:_type] isa AbstractDict && haskey(dict[:_type], :name)
      type_name = dict[:_type][:name]
    else
      type_name = dict[:_type]
    end
    # Upgrades
    if type_name == "Bool"
      dict[:data] = parse(Bool, dict[:data])
    elseif type_name == "MatroidRealizationSpace"
      dict[:data][:one_realization] = parse(Bool, dict[:data][:one_realization])
    elseif type_name == "SelfProjectingMatroidRealizations"
      dict[:data][:equality_of_realizationspaces] = parse(Bool, dict[:data][:equality_of_realizationspaces])

      upgraded_mrs = upgrade_1_8_0(s, Dict{Symbol, Any}(
        :_type => Dict(:name => "MatroidRealizationSpace",
                       :params => Dict(
                         :matrix_space => dict[:_type][:params][:matrix_space_mrs],
                         :ideal_ring => dict[:_type][:params][:ideal_ring_rs],
                         :ground_ring => dict[:_type][:params][:ground_ring]
                        )),
        :data => dict[:data][:realization_space]
      ))
      dict[:data][:realization_space] = upgraded_mrs[:data]
    elseif type_name == "IdealGens"
      dict[:data][:is_gb] = parse(Bool, dict[:data][:is_gb])
      dict[:data][:is_reduced] = parse(Bool, dict[:data][:is_reduced])
      dict[:data][:keep_ordering] = parse(Bool, dict[:data][:keep_ordering])
    elseif type_name == "MonomialOrdering"
      if haskey(dict[:data], :is_total)
        dict[:data][:is_total] = parse(Bool, dict[:data][:is_total])
      end
    end
    return dict
  end
))
