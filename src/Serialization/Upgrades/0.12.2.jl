################################################################################
# Upgrade Summary
# The purpose of this upgrade is to make polynomial serialization more efficient.
# the terms of a polynomial are now stored as
#    terms = [(exponent, coeff), ..., (exponent, coeff)]
# where exponent is either an int or a vector of ints, and coeff is either a string
# or is given as a vector of terms. This results in a nested structure, the parents
# for the coefficients are found in the parents vector. The parent of the
# coefficient can be determined based on the depth level of the given coefficient
# This new format also allows us to improve serialization of vectors of ring
# elements, since the parents vector can be brought to the root, and we use
# the vector entries to store the terms

##### PLEASE NOTE #####
# This upgrade script does not cover all cases, only what is necessary for the
# upgrade tests to pass. Univariate polynomial and multivariate serialization
# has changed drastically as well as any types that depend on such types,
# and covering all cases would take a large amount of time and effort,
# and we feel at this time not worth it. However, if absolute necessary
# an upgrade for a user can be made upon request.


push!(upgrade_scripts_set, UpgradeScript(
  v"0.12.2",
  function upgrade_0_12_2(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    s.nested_level += 1
    # moves down tree to point where type exists in dict
    # since we are only doing updates based on certain types
    # no :type key implies the dict is data
    if !haskey(dict, :type)
      return upgrade_data(upgrade_0_12_2, s, dict)
    end

    upgraded_dict = dict
    # we only upgrade polynomials and MPolyIdeals
    if upgraded_dict[:type] == "MPolyIdeal"
      parent = dict[:data][:parent]
      s.id_to_dict[Symbol(parent[:id])] = parent
      upgraded_gens = []
      for gen in dict[:data][:gens][:data][:vector]
        push!(upgraded_gens, upgrade_0_12_2(s, gen))
      end
    elseif contains(upgraded_dict[:type], "PolyRingElem")
      if !haskey(upgraded_dict[:data], :parent)
        return upgraded_dict
      end

      # we first create the parent list
      parent_dict = upgraded_dict[:data][:parent]
      parents = []
      if haskey(parent_dict, :id)
        if parent_dict[:type] == "#backref"
          parent_dict = s.id_to_dict[Symbol(parent_dict[:id])]
        else
          s.id_to_dict[Symbol(parent_dict[:id])] = parent_dict
        end
        push!(parents, parent_dict)
      end

      while haskey(parent_dict[:data], :base_ring)
        parent_dict = parent_dict[:data][:base_ring]
        if !haskey(parent_dict, :id)
          break
        end
        if parent_dict[:type] == "#backref"
          parent_dict = s.id_to_dict[Symbol(parent_dict[:id])]
        else
          s.id_to_dict[Symbol(parent_dict[:id])] = parent_dict
        end
        push!(parents, parent_dict)
        if haskey(parent_dict[:data], :def_pol)
          parent_dict = parent_dict[:data][:def_pol][:data][:parent]
          push!(parents, parent_dict)
        end
      end

      local_refs = Dict()
      upgraded_parents = []

      # using the parent list we attach the refs 
      # to the root of the file
      for parent in parents
        id = parent[:id]
        
        pushfirst!(upgraded_parents, Dict(:id => id, :type => "#backref"))
        if haskey(local_refs, Symbol(id))
          continue
        end
        if haskey(parent[:data], :def_pol)
          upgraded_def_pol = upgrade_0_12_2(s, parent[:data][:def_pol])
          parent[:data][:def_pol] = upgraded_def_pol
        elseif haskey(parent[:data], :base_ring)
          base_dict = parent[:data][:base_ring]

          if haskey(base_dict, :id)
            parent[:data][:base_ring] = Dict(:id => base_dict[:id], :type => "#backref")
          end
        end
        local_refs[Symbol(id)] = parent
      end
      upgraded_dict[:data][:parents] = upgraded_parents
      terms = []
              
      if upgraded_dict[:type] == "MPolyRingElem"
        # convert the terms array to the new format
        for term in dict[:data][:terms]
          if term[:coeff][:type] == "#backref"
            term[:coeff] = s.id_to_dict[Symbol(term[:coeff][:id])]
          else
            s.id_to_dict[Symbol(term[:coeff][:id])] = term[:coeff]
          end
          push!(terms, [term[:exponent][:data][:vector],
                        term[:coeff][:data][:data]])
        end
      else
        # we convert the vector of coefficients
        # to the new terms format
        exponent = 0
        for coeff in upgraded_dict[:data][:coeffs][:data][:vector]
          if !isa(coeff, AbstractDict)
            upgraded_coeff = coeff
          elseif haskey(coeff[:data], :polynomial)
            upgraded_coeff = upgrade_0_12_2(s, coeff[:data][:polynomial])[:data][:terms]
          else
            upgraded_coeff = coeff[:data][:data]
          end
          term = (exponent, upgraded_coeff)
          push!(terms, term)
          exponent += 1
        end
      end
      upgraded_dict[:data] = Dict(:terms => terms, 
                                  :parents => upgraded_parents)
      merge!(s.id_to_dict, local_refs)
    end

    s.nested_level -= 1
    if !isempty(s.id_to_dict) && s.nested_level == 0
      upgraded_dict[:refs] = s.id_to_dict
    end

    return upgraded_dict
  end
))
