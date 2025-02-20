function upgrade_NamedTuple(d::Dict)
  for (i, entry) in enumerate(d[:data])
    if entry isa Dict
      func = Symbol(string("upgrade_", entry[:_type][:name]))
      upgraded_entry = @eval $func($dict)
    end
  end
end

function upgrade_QSMModel(d::Dict)
  upgraded_dict = d
  upgraded_dict[:_type] = Dict(
    :name => d[:_type],
    :params => Dict(
      :hs_model => d[:data][:hs_model][:_type],
      :genus_ci => d[:data][:genus_ci][:_type],
      :degree_of_Kbar_of_tv_restricted_to_ci => d[:data][:degree_of_Kbar_of_tv_restricted_to_ci][:_type]
    ))
  upgraded_dict[:data][:hs_model] = d[:data][:hs_model][:data]
  upgraded_dict[:data][:genus_ci] = d[:data][:genus_ci][:data]
  upgraded_dict[:data][:degree_of_Kbar_of_tv_restricted_to_ci] = d[:data][:degree_of_Kbar_of_tv_restricted_to_ci][:data]

  return upgraded_dict
end

push!(upgrade_scripts_set, UpgradeScript(
  v"1.3.0",
  function upgrade_1_3_0(s::UpgradeState, dict::Dict)
    if haskey(dict, :_type)
      if dict[:_type] isa Dict && haskey(dict[:_type], :name)
        type_name = dict[:_type][:name]
      else
        type_name = dict[:_type]
      end
      upgraded_dict = dict
      if type_name in ["PolyRing", "FreeAssociativeAlgebra", "MPolyRing"]
        upgraded_dict[:_type] = Dict(
          :name => dict[:_type],
          :params => dict[:data][:base_ring]
        )
        upgraded_dict[:data] = Dict(
          :symbols => dict[:data][:symbols],
        )
      elseif type_name == "MatSpace"
        upgraded_dict[:_type] = Dict(
          :name => type_name,
          :params => upgraded_dict[:data][:base_ring]
        )
        upgraded_dict[:data] = Dict(
          :ncols => upgraded_dict[:data][:ncols],
          :nrows => upgraded_dict[:data][:nrows]
        )
      elseif type_name == "fqPolyRepField"
        upgraded_dict[:_type] = Dict(
          :name => dict[:_type],
          :params => dict[:data][:def_pol][:_type][:params]
        )
        upgraded_dict[:data] = dict[:data][:def_pol][:data]
      elseif type_name == "MPolyDecRing"
        upgraded_dict[:_type] = Dict(
          :name => d[:_type],
          :params => Dict(
            :ring => d[:data][:ring],
            :grading_group => d[:data][:grading][:_type][:params]
          )
        )
        upgraded_dict[:data] = d[:data][:grading][:data]
      elseif type_name == "FqField"
        if dict[:data] isa Dict
          upgraded_dict[:_type] = Dict(
            :name => dict[:_type],
            :params => dict[:data][:def_pol][:_type][:params]
          )
          upgraded_dict[:data] = dict[:data][:def_pol][:data]
        end
      elseif type_name == "Tuple"
        upgraded_subtypes = Dict[]
        for (i, subtype) in enumerate(dict[:_type][:params])
          push!(upgraded_subtypes, upgrade_1_3_0(s, Dict(
            :_type => subtype,
            :data => dict[:data][i]
          )))
        end
        upgraded_dict[:_type][:params] = [subdict[:_type] for subdict in upgraded_subtypes]
        upgraded_dict[:data] = [subdict[:data] for subdict in upgraded_subtypes]
      elseif type_name == "ZZLat"
        upgraded_dict[:_type] = Dict(
          :name => type_name,
          :params => Dict(
            :basis => dict[:data][:basis][:_type],
            :ambient_space => dict[:data][:ambient_space]
          )
        )
        upgraded_dict[:data] = dict[:data][:basis][:data]
      elseif type_name in [
        "PolyRingElem", "MPolyRingElem", "NormalToricVariety", "FinGenAbGroup",
        "MPolyIdeal", "Dict", "MatElem"
        ]
        # do nothing
        upgraded_dict = dict
      else
        error("$type_name doesn't have upgrade")
      end        
    elseif haskey(dict, :data) && dict[:data] isa Dict
      upgraded_dict[:data] = upgrade_1_3_0(s, dict[:data])
    end

    if haskey(dict, :_refs)
      upgraded_refs = Dict()
      for (k, v) in dict[:_refs]
        upgraded_refs[k] = upgrade_1_3_0(s, v)
      end
      upgraded_dict[:_refs] = upgraded_refs
    end

    return upgraded_dict
  end
))
