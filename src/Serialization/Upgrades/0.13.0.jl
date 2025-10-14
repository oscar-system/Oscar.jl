################################################################################
# Upgrade Summary
# This upgrade adds the param key word to types, it also changes any backref types
# to plain strings.

##### PLEASE NOTE #####
# This upgrade script does not cover all cases, only what is necessary for the
# upgrade tests to pass. As well as any tests that appear now in the data test
# section

# used to replace Ints in the file with strings
function upgrade_terms(terms)
  upgraded_terms = []
  if terms isa String
    return terms
  end
  if terms isa Int
    return string(terms)
  end
  for term in terms
    if term[1] isa Int
      upgraded_term = [string(term[1]), upgrade_terms(term[2])]
    elseif term[1] isa Vector{Int}
      upgraded_term = [
        [string(t) for t in term[1]],
        upgrade_terms(term[2])
      ]
    else
      upgraded_term = term
    end
    if term[2] isa Int
      upgraded_term[2] = string(term[2])
    end
    push!(upgraded_terms, upgraded_term)
  end
  return upgraded_terms
end

push!(upgrade_scripts_set, UpgradeScript(
  v"0.13.0",
  function upgrade_0_13_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    upgraded_dict = dict
    if !haskey(dict, :type)
      # this section deals will PolyRings and MPolyRings
      if haskey(dict, :base_ring) && dict[:base_ring] isa AbstractDict
        if dict[:base_ring][:type] == "#backref"
          upgraded_dict[:base_ring] = dict[:base_ring][:id]
        else
          if haskey(dict[:base_ring], :data)
            # should have only one key, since any base_ring that
            # doesn't use an id has one parameter
            key = first(keys(dict[:base_ring][:data]))
            upgraded_dict[:base_ring][:data] = string(dict[:base_ring][:data][key])
          end
          upgraded_dict[:base_ring][:_type] = dict[:base_ring][:type]
        end
        if haskey(dict, :symbols)
          if dict[:symbols] isa AbstractDict
            upgraded_dict[:symbols] = dict[:symbols][:data][:vector]
          end
        end
        return upgraded_dict
      end
      return upgraded_dict
    end

    if dict[:type] == "#backref"
      return dict[:id]
    end

    if dict[:type] == "Nemo.fpFieldElem"
      dict[:_type] = "fpFieldElem"
    end

    if dict[:type] isa String
      type_string = dict[:type]
    else
      type_string = dict[:type][:name]
    end

    # type has already been updated
    if haskey(dict, :_type) && dict[:_type] isa AbstractDict
      return dict
    end
    if contains(type_string, "Vector")
      # this section wasn't necessary for upgrading the folder of surfaces
      # which was our primary focus, this may be updated upon request
      error("The upgrade script needs an update")
    elseif contains(type_string, "PolyRingElem")
      if dict[:data] isa AbstractDict && haskey(dict[:data], :parents)
        # this currently only handles one parent
        upgraded_parent = dict[:data][:parents][end]
        params = upgraded_parent[:id]
        upgraded_dict[:data] = upgrade_terms(dict[:data][:terms])
      else
        # this section wasn't necessary for upgrading the folder of surfaces
        # which was our primary focus, this may be updated upon request
        error("The upgrade script needs an update")
      end

      upgraded_dict[:_type] = Dict{Symbol, Any}(:name => type_string, :params => params)
    elseif contains(type_string, "NamedTuple")
      upgraded_tuple = upgrade_0_13_0(s, dict[:data][:content])

      params = Dict{Symbol, Any}(
        :tuple_params => upgraded_tuple[:_type][:params],
        :names => dict[:data][:keys][:data][:content]
      )
      upgraded_dict[:_type] = Dict{Symbol, Any}(:name => type_string, :params => params)
      upgraded_dict[:data] = upgraded_tuple[:data]

    elseif contains(type_string, "Nemo.fpField")
      if dict[:data] isa String
        return dict
      elseif dict[:data] isa Int
        upgraded_dict[:data] = string(dict[:data])
        return upgraded_dict
      end
      upgraded_dict[:data] = string(dict[:data][:characteristic])
    elseif contains(type_string, "MPolyIdeal")
      upgraded_parent = upgrade_0_13_0(s, dict[:data][:parent])[:id]
      upgraded_dict[:_type] = Dict{Symbol, Any}(:name => type_string, :params => upgraded_parent)
      upgraded_gens = []
      for gen in dict[:data][:gens][:data][:vector]
        push!(upgraded_gens, upgrade_0_13_0(s, gen)[:data])
      end
      upgraded_dict[:data] = upgraded_gens
    elseif contains(type_string, "MPolyRing")
      if dict[:data][:base_ring] isa AbstractDict
        if dict[:data][:base_ring][:type] == "#backref"
          upgraded_dict[:data][:base_ring] = dict[:data][:base_ring][:id]
        end
      end
      # has already been updated
      if !(dict[:data][:symbols] isa AbstractDict)
        return dict
      end
      upgraded_dict[:data][:symbols] = dict[:data][:symbols][:data][:vector]
          # update the data section for specific types
    elseif contains(type_string, "PolyRing")
      upgraded_dict[:data] = upgrade_0_13_0(s, dict[:data])
    elseif contains(type_string, "Tuple")
      params = []
      entry_data = []

      for (i, field_type) in enumerate(dict[:data][:field_types])
        if !(dict[:data][:content][i] isa String)
          upgraded_entry = upgrade_0_13_0(s, dict[:data][:content][i])
          push!(params, upgraded_entry[:_type])
          push!(entry_data, upgraded_entry[:data])
        else
          push!(params, field_type)
          push!(entry_data, dict[:data][:content][i])
        end
      end
      upgraded_dict[:_type] = Dict{Symbol, Any}(:name => type_string, :params => params)
      upgraded_dict[:data] = entry_data
    elseif contains(type_string, "Matrix")
      params = dict[:data][:matrix][1][:data][:entry_type]
      upgraded_dict[:_type] = Dict{Symbol, Any}(:name => "Matrix", :params => params)
      matrix_data = []
      for v in dict[:data][:matrix]
        push!(matrix_data, v[:data][:vector])
      end
      upgraded_dict[:data] = matrix_data
    end

    if contains(type_string, "Field")
      if dict[:data] isa AbstractDict && haskey(dict[:data], :def_pol)
        upgraded_dict[:data][:def_pol] = upgrade_0_13_0(s, dict[:data][:def_pol])
      end
    end

    if haskey(upgraded_dict, :refs)
      upgraded_refs = Dict{Symbol, Any}()
      for (k, v) in upgraded_dict[:refs]
        upgraded_refs[k] = upgrade_0_13_0(s, v)
      end
      upgraded_dict[:_refs] = upgraded_refs
    end

    if haskey(upgraded_dict, :type) && !haskey(upgraded_dict, :_type)
      upgraded_dict[:_type] = upgraded_dict[:type]
    end
    return upgraded_dict
  end
))
