# Upgrades for NamedTuple, IdealGens, FracField, MPolyQuoRing, polyhedral object, and Dict encoding changes:
# - NamedTuple: {names:[...], tuple_params:[...]} with array data changed to
#   {names:[...], entry_params:{name: type}} with data keyed by name
# - IdealGens params changed from {base_ring:<typed_obj>, ordering_type:"MonomialOrdering"} to <typed_obj> directly
# - FracField data no longer stores base_ring (redundant with type params)
# - MPolyQuoRing params simplified from {base_ring:..., ordering:...} to base_ring directly
# - Polyhedral objects with complex fields: pm_params wrapper flattened to top-level params
# - Heterogeneous Dict: per-key value type params moved from top-level into entry_params
push!(upgrade_scripts_set, UpgradeScript(
  v"1.9.0+1",
  function upgrade_1_9_0_1(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    upgrade_recursive(upgrade_1_9_0_1, s, dict)

    if dict[:_type] isa AbstractDict && haskey(dict[:_type], :name)
      type_name = dict[:_type][:name]
    else
      type_name = dict[:_type]
    end

    if type_name == "NamedTuple"
      old_params = dict[:_type][:params]
      if old_params isa AbstractDict && haskey(old_params, :tuple_params)
        names = Symbol.(old_params[:names])
        data = get(dict, :data, [])
        dict[:_type][:params] = Dict{Symbol, Any}(
          :names => old_params[:names],
          :entry_params => Dict{Symbol, Any}(zip(names, old_params[:tuple_params]))
        )
        dict[:data] = Dict{Symbol, Any}(zip(names, data))
      end
    elseif type_name == "IdealGens"
      old_params = dict[:_type][:params]
      if old_params isa AbstractDict && haskey(old_params, :base_ring)
        dict[:_type][:params] = old_params[:base_ring]
      end
    elseif type_name == "FracField"
      if dict[:data] isa AbstractDict && haskey(dict[:data], :base_ring)
        delete!(dict[:data], :base_ring)
      end
    elseif type_name == "MPolyQuoRing"
      old_params = dict[:_type][:params]
      if old_params isa AbstractDict && haskey(old_params, :base_ring)
        dict[:_type][:params] = old_params[:base_ring]
      end
    elseif type_name == "Dict"
      old_params = dict[:_type][:params]
      if old_params isa AbstractDict && haskey(old_params, :key_params) &&
         !haskey(old_params, :value_params) && !haskey(old_params, :entry_params)
        dict[:_type][:params] = Dict{Symbol, Any}(
          :key_params => old_params[:key_params],
          :entry_params => filter(kv -> kv.first != :key_params, old_params)
        )
      end
    elseif type_name in ("Polyhedron", "Cone", "PolyhedralComplex", "PolyhedralFan",
                         "SubdivisionOfPoints", "LinearProgram", "MixedIntegerLinearProgram")
      if dict[:_type] isa AbstractDict &&
         haskey(dict[:_type], :params) &&
         dict[:_type][:params] isa AbstractDict &&
         haskey(dict[:_type][:params], :pm_params)
        pm_params = dict[:_type][:params][:pm_params]
        if pm_params isa AbstractDict && haskey(pm_params, :params)
          inner = pm_params[:params]
          # The Dict upgrade (above) runs on pm_params first (it's a heterogeneous dict),
          # so property types may already be wrapped in entry_params.
          prop_types = if haskey(inner, :entry_params)
            inner[:entry_params]
          else
            filter(kv -> kv.first != :key_params, inner)
          end
          for (k, v) in prop_types
            dict[:_type][:params][k] = v
          end
        end
        delete!(dict[:_type][:params], :pm_params)
      end
    elseif type_name == "PhylogeneticModel"
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
      if dict[:data] isa AbstractDict && haskey(dict[:data], :model_parameter_name)
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
      if dict[:data] isa AbstractDict && haskey(dict[:data], :varnames_group_based)
        dict[:data][:model_parameter] = dict[:data][:varnames_group_based]
        delete!(dict[:data], :varnames_group_based)
      end
    elseif type_name in ("ZZLat", "ZZLatWithIsom")
      p = get(get(dict, :_type, Dict{Symbol, Any}()), :params, nothing)
      if p isa AbstractDict && haskey(p, :basis) &&
         p[:basis] isa AbstractDict && get(p[:basis], :name, nothing) == "MatElem"
        p[:basis] = p[:basis][:params]
      end
    elseif type_name == "QuadSpaceWithIsom"
      p = get(get(dict, :_type, Dict{Symbol, Any}()), :params, nothing)
      if p isa AbstractDict && haskey(p, :isom) &&
         p[:isom] isa AbstractDict && get(p[:isom], :name, nothing) == "MatElem"
        p[:isom] = p[:isom][:params]
      end
    end
    return dict
  end
))
