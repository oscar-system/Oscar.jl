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
      elseif type_name in ["AbsSimpleNumField", "Hecke.RelSimpleNumField"]
        upgraded_dict[:_type] = Dict(
          :name => dict[:_type],
          :params => dict[:data][:def_pol][:_type][:params]
        )
        upgraded_dict[:data] = Dict(
          :var => dict[:data][:var],
          :def_pol => dict[:data][:def_pol][:data]
        )
      elseif type_name in ["AbsNonSimpleNumField", "Hecke.RelNonSimpleNumField"]
        upgraded_dict[:_type] = Dict(
          :name => dict[:_type],
          :params => dict[:data][:def_pols][:_type][:params]
        )
        upgraded_dict[:data] = Dict(
          :vars => dict[:data][:vars],
          :def_pols => dict[:data][:def_pols][:data]
        )

      elseif type_name == "EmbeddedNumField"
        upgraded_dict[:_type] = Dict(
          :name => dict[:_type],
          :params => Dict(
            :embedding => dict[:data][:embedding],
            :num_field => dict[:data][:num_field]
        )
        upgraded_dict[:data] = []
      elseif type_name == "FqField"
        if dict[:data] isa Dict
          upgraded_dict[:_type] = Dict(
            :name => dict[:_type],
            :params => dict[:data][:def_pol][:_type][:params]
          )
          upgraded_dict[:data] = dict[:data][:def_pol][:data]
        end
      elseif type_name == "Dict"
        d = Dict()
        for (k, v) in dict[:_type][:params]
          if k == :key_type
            d[k] = v
          else
            d[k] = upgrade_1_3_0(s, Dict(
              :_type => dict[:_type][:params][k],
              :data => dict[:data][k]
            ))
          end
        end
        upgraded_dict = Dict(
          :_type => Dict(:name => "Dict", :params=>Dict()),
          :data => Dict()
        )
        for (k, v) in d
          if k == :key_type
            upgraded_dict[:_type][:params][k] = v
          else
            upgraded_dict[:_type][:params][k] = v[:_type]
            upgraded_dict[:data][k] = v[:data]
          end
        end
      elseif type_name in ["Vector", "Set"]
        subtype = dict[:_type][:params]
        upgraded_entries = []
        for entry in dict[:data]
          push!(upgraded_entries, upgrade_1_3_0(s, Dict(
            :_type => subtype,
            :data => entry
          )))
        end
        upgraded_dict[:data] = [e[:data] for e in upgraded_entries]
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
        "MPolyIdeal", "MatElem", "String", "Base.Int", "Bool", "Graph{Undirected}",
        "Graph{Directed}", "Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric}"
        ]
        # do nothing
        upgraded_dict = dict
      elseif type_name == "Polyhedron"
        if !(dict[:_type][:params] == "QQField")
          upgraded_subdict = upgrade_1_3_0(s, dict[:data])
          upgraded_subdict[:_type][:params][:key_type] = "Symbol"
          upgraded_dict[:_type] = Dict(
            :name => type_name,
            :params => upgraded_subdict[:_type]
          )
          upgraded_dict[:data] = upgraded_subdict[:data]

          upgraded_dict[:_type][:params][:params][:_polymake_type] = dict[:_type][:params][:params][:_type]
          upgraded_dict[:data][:_polymake_type] = dict[:data][:_type]
        else
          upgraded_dict = dict
        end
      elseif type_name in [
        "Hecke.RelSimpleNumFieldEmbedding", "Hecke.RelNonSimpleNumFieldEmbedding"
        ]
        upgraded_dict[:_type] = Dict(
          :name => type_name,
          :params => Dict(
            :num_field => dict[:data][:num_field],
            :base_field_emb => dict[:data][:base_field_emb]
          )
        )
        upgraded_dict[:data] = dict[:data][:data][:data]
      elseif type_name in  [
        "Hecke.AbsSimpleNumFieldEmbedding", "Hecke.AbsNonSimpleNumFieldEmbedding"
        ]
        upgraded_dict[:_type] = Dict(
          :name => type_name,
          :params => dict[:data][:num_field]
        )
        upgraded_dict[:data] = dict[:data][:data][:data]
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
