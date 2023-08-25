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
  function upgrade_0_13_0(s::UpgradeState, dict::Dict)
    upgraded_dict = dict
    if !haskey(dict, :type)
      # this section deals will PolyRings and MPolyRings
      if haskey(dict, :base_ring) && dict[:base_ring] isa Dict
        if dict[:base_ring][:type] == "#backref"
          upgraded_dict[:base_ring] = dict[:base_ring][:id]
        elseif haskey(dict[:base_ring], :data)
          # should have only one key, since any base_ring that
          # doesn't use an id has one paramater
          key = first(keys(dict[:base_ring][:data]))
          upgraded_dict[:base_ring][:data] = string(dict[:base_ring][:data][key])
        end
        if haskey(dict, :symbols)
          if dict[:symbols] isa Dict
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
      dict[:type] = "fpFieldElem"
    end
    T = decode_type(dict[:type])

    # type has already been updated
    if dict[:type] isa Dict
      return dict
    end
    if T <: Vector
      # this section wasn't necessary for upgrading the folder of surfaces
      # which was our primary focus, this may be updated upon request
      throw(Error("The upgrade script needs an update"))
    elseif T <: Union{PolyRingElem, MPolyRingElem}
      if dict[:data] isa Dict && haskey(dict[:data], :parents)
        upgraded_parents = Any[
          p[:id] for p in dict[:data][:parents]
            ]
        params = upgraded_parents
        upgraded_dict[:data] = upgrade_terms(dict[:data][:terms])
      else
        # this section wasn't necessary for upgrading the folder of surfaces
        # which was our primary focus, this may be updated upon request
        throw(Error("The upgrade script needs an update"))
      end

      upgraded_dict[:type] = Dict(:name => "PolyRingElem", :params => params)
    elseif T <: NamedTuple
      upgraded_tuple = upgrade_0_13_0(s, dict[:data][:content])
      
      params = Dict(
        :tuple_params => upgraded_tuple[:type][:params],
        :names => dict[:data][:keys][:data][:content]
      )
      upgraded_dict[:type] = Dict(:name => "NamedTuple", :params => params)
      upgraded_dict[:data] = upgraded_tuple[:data]

    elseif T <: Nemo.fpField
      if dict[:data] isa String
        return dict
      elseif dict[:data] isa Int
        upgraded_dict[:data] = string(dict[:data])
        return upgraded_dict
      end
      upgraded_dict[:data] = string(dict[:data][:characteristic])
    elseif T <: MPolyIdeal
      upgraded_parents = upgrade_0_13_0(s, dict[:data][:parent])
      upgraded_dict[:type] = Dict(:name => "MPolyIdeal", :params => upgraded_parents)
      upgraded_gens = []
      for gen in dict[:data][:gens][:data][:vector]
        push!(upgraded_gens, upgrade_0_13_0(s, gen)[:data])
      end
      upgraded_dict[:data] = upgraded_gens
    elseif T <: MPolyRing
      if dict[:data][:base_ring] isa Dict
        if dict[:data][:base_ring][:type] == "#backref"
          upgraded_dict[:data][:base_ring] = dict[:data][:base_ring][:id]
        end
      end
      # has already been updated
      if !(dict[:data][:symbols] isa Dict)
        return dict
      end
      upgraded_dict[:data][:symbols] = dict[:data][:symbols][:data][:vector]
    elseif T <: Tuple
      params = []
      entry_data = []
      
      for (i, field_type) in enumerate(dict[:data][:field_types])
        U = decode_type(field_type)
        if type_needs_params(U)
          upgraded_entry = upgrade_0_13_0(s, dict[:data][:content][i])
          push!(params, upgraded_entry[:type])
          push!(entry_data, upgraded_entry[:data])
        else
          push!(params, field_type)
          push!(entry_data, dict[:data][:content][i])
        end
      end
      upgraded_dict[:type] = Dict(:name => "Tuple", :params => params)
      upgraded_dict[:data] = entry_data
    elseif  T <: Matrix
      params = dict[:data][:matrix][1][:data][:entry_type]
      upgraded_dict[:type] = Dict(:name => "Matrix", :params => params)
      matrix_data = []
      for v in dict[:data][:matrix]
        push!(matrix_data, v[:data][:vector])
      end
      upgraded_dict[:data] = matrix_data
    end

    # update the data section for specific types
    if T <: PolyRing
      upgraded_dict[:data] = upgrade_0_13_0(s, dict[:data])
    end
    
    if T <: Field
      if dict[:data] isa Dict &&haskey(dict[:data], :def_pol)
        upgraded_dict[:data][:def_pol] = upgrade_0_13_0(s, dict[:data][:def_pol])
      end
    end
    
    if haskey(upgraded_dict, :refs)
      upgraded_refs = Dict()
      for (k, v) in upgraded_dict[:refs]
        upgraded_refs[k] = upgrade_0_13_0(s, v)
      end
      upgraded_dict[:refs] = upgraded_refs
    end
    return upgraded_dict
  end
))
