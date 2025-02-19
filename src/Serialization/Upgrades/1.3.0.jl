function upgrade_NamedTuple(d::Dict)
  for (i, entry) in enumerate(d[:data])
    if entry isa Dict
      func = Symbol(string("upgrade_", entry[:_type][:name]))
      upgraded_entry = @eval $func($dict)
    end
  end
end

function upgrade_ring(d::Dict)
  upgraded_dict = d
  upgraded_dict[:_type] = Dict(
    :name => d[:_type],
    :params => d[:data][:base_ring]
  )
  upgraded_dict[:data] = Dict(
    :symbols => d[:data][:symbols],
  )
  return upgraded_dict
end


function upgrade_Dict(d::Dict)
  return d
end

function upgrade_fqPolyRepField(d::Dict)
  upgraded_dict = d
  upgraded_dict[:_type] = Dict(
    :name => d[:_type],
    :params => d[:data][:def_pol][:_type][:params]
  )
  upgraded_dict[:data] = d[:data][:def_pol][:data]
  return upgraded_dict
end

function upgrade_FqField(d::Dict)
  upgraded_dict = d
  if d[:data] isa Dict
    upgraded_dict[:_type] = Dict(
      :name => d[:_type],
      :params => d[:data][:def_pol][:_type][:params]
    )
    upgraded_dict[:data] = d[:data][:def_pol][:data]
  end
  return upgraded_dict
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

  println(json(upgraded_dict[:_type], 2))
  return upgraded_dict
end

function upgrade_MPolyDecRing(d::Dict)
  upgraded_dict = d
  upgraded_dict[:_type] = Dict(
    :name => d[:_type],
    :params => Dict(
      :ring => d[:data][:ring],
      :grading_group => d[:data][:grading][:_type][:params]
    )
  )
  upgraded_dict[:data] = d[:data][:grading][:data]
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

      if type_name in ["PolyRing", "FreeAssociativeAlgebra", "MPolyRing"]
        upgraded_dict = upgrade_ring(dict)
      elseif type_name in [
        "PolyRingElem", "MPolyRingElem", "NormalToricVariety", "FinGenAbGroup"]
        # do nothing
        upgraded_dict = dict
      else
        func = Symbol(string("upgrade_", type_name))
        upgraded_dict = @eval $func($dict)
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
