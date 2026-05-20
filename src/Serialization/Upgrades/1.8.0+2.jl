# Upgrades for NamedTuple and IdealGens encoding changes:
# - NamedTuple params changed from {names:[...], tuple_params:[...]} to flat field-name dict
# - IdealGens params changed from {base_ring:<typed_obj>, ordering_type:"MonomialOrdering"} to <typed_obj> directly
push!(upgrade_scripts_set, UpgradeScript(
  v"1.8.0+2",
  function upgrade_1_8_0_2(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    upgrade_recursive(upgrade_1_8_0_2, s, dict)

    if dict[:_type] isa AbstractDict && haskey(dict[:_type], :name)
      type_name = dict[:_type][:name]
    else
      type_name = dict[:_type]
    end

    if type_name == "NamedTuple"
      old_params = dict[:_type][:params]
      if old_params isa AbstractDict && haskey(old_params, :names)
        new_params = Dict{Symbol, Any}()
        for (name, tp) in zip(old_params[:names], old_params[:tuple_params])
          new_params[Symbol(name)] = tp
        end
        dict[:_type][:params] = new_params
      end
    elseif type_name == "IdealGens"
      old_params = dict[:_type][:params]
      if old_params isa AbstractDict && haskey(old_params, :base_ring)
        dict[:_type][:params] = old_params[:base_ring]
      end
    end
    return dict
  end
))
