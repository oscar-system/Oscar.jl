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
    push!(upgraded_terms, upgraded_term)
  end
  return upgraded_terms
end




push!(upgrade_scripts_set, UpgradeScript(
  v"0.13.0",
  function upgrade_0_13_0(refs::Dict, dict::Dict)
    upgraded_dict = dict
    if !haskey(dict, :type)
      if haskey(dict, :base_ring) && dict[:base_ring] isa Dict
        if dict[:base_ring][:type] == "#backref"
          upgraded_dict[:base_ring] = dict[:base_ring][:id]
        elseif haskey(dict[:base_ring], :data)
          # should have only one key, since any base_ring that
          # doesn't use an id has one one paramater
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

    if dict[:type] isa Dict
      #println(json(dict, 2))
    end

    T = decode_type(dict[:type])
    # Types that now have params that haven't already been updated
    if type_needs_params(T)
      # type has already been updated
      if dict[:type] isa Dict
        return dict
      end
      if T <: Vector
        params = Dict(:params => Dict())
      elseif T <: Union{PolyRingElem}
        if dict[:data] isa Dict && haskey(dict[:data], :parents)
          upgraded_parents = Any[
            p[:id] for p in dict[:data][:parents]
              ]
          params = upgraded_parents
          upgraded_dict[:data] = upgrade_terms(dict[:data][:terms])
        else
          params = "this needs to be filled in"
        end
      end

      upgraded_dict[:type] = Dict(
        :name => dict[:type],
        :params => params
      )
    end

    if T <: PolyRing
      upgraded_dict[:data] = upgrade_0_13_0(refs, dict[:data])
    elseif T <: Field
      if haskey(dict[:data], :def_pol)
        upgraded_dict[:data][:def_pol] = upgrade_0_13_0(refs, dict[:data][:def_pol])
      end
    end
    
    if haskey(upgraded_dict, :refs)
      upgraded_refs = Dict()
      for (k, v) in upgraded_dict[:refs]
        upgraded_refs[k] = upgrade_0_13_0(refs, v)
      end
      upgraded_dict[:refs] = upgraded_refs
    end
    return upgraded_dict
  end
))
