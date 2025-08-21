push!(upgrade_scripts_set, UpgradeScript(
  v"1.5.0",
  function upgrade_1_5_0(s::UpgradeState, dict::Dict)
    if haskey(dict, :_refs)
      s.id_to_dict = dict[:_refs]
    end
    if haskey(dict, :_ns)
      if haskey(dict[:_ns], :polymake)
        return dict
      end
    end
    if haskey(dict, :_type)

      if dict[:_type] isa Dict && haskey(dict[:_type], :name)
        type_name = dict[:_type][:name]
      else
        type_name = dict[:_type]
      end
      upgraded_dict = dict

      if type_name in ["GlobalTateModel", "HypersurfaceModel", "WeierstrassModel"]

        upgraded_attr_dict = upgrade_1_5_0(s, dict[:data][:__attrs])
        upgraded_dict[:attrs] = Dict()
        for k in keys(upgraded_attr_dict[:_type][:params])
          k == :key_params && continue
          upgraded_dict[:attrs][k] = Dict(:_type => upgraded_attr_dict[:_type][:params][k], :data => upgraded_attr_dict[:data][k])
        end

        if type_name == "WeierstrassModel" && haskey(dict[:data], :weierstrass_polynomial)
          upgraded_poly = upgrade_1_4_0(s, dict[:data][:weierstrass_polynomial])
          upgraded_dict[:_type][:params][:hypersurface_equation_ring] = upgraded_poly[:_type]
          upgraded_dict[:data][:hypersurface_equation] = upgraded_poly[:data]
        end

        if type_name == "GlobalTateModel" && haskey(dict[:data], :tate_polynomial)
          upgraded_poly = upgrade_1_4_0(s, dict[:data][:tate_polynomial])
          upgraded_dict[:_type][:params][:hypersurface_equation_ring] = upgraded_poly[:_type]
          upgraded_dict[:data][:hypersurface_equation] = upgraded_poly[:data]
        end

      end

    elseif haskey(dict, :data) && dict[:data] isa Dict
      upgraded_dict[:data] = upgrade_1_5_0(s, dict[:data])
    end

    if haskey(dict, :_refs)
      upgraded_refs = Dict()
      for (k, v) in dict[:_refs]
        upgraded_refs[k] = upgrade_1_5_0(s, v)
      end
      upgraded_dict[:_refs] = upgraded_refs
    end
    return upgraded_dict
  end
))
