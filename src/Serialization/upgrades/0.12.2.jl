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


push!(upgrade_scripts, UpgradeScript(
  v"0.12.2",
  function upgrade_0_12_2(s::DeserializerState, dict::Dict)
    # moves down tree to point where type exists in dict
    # since we are only doing updates based on certain types
    # no :type key implies the dict is data
    if !haskey(dict, :type)
      return upgrade_data(upgrade_0_12_2, s, dict)
    end
    
    upgraded_dict = dict
    # we only upgrade univariate polynomials
    if contains(upgraded_dict[:type], "PolyRingElem")
      if !haskey(upgraded_dict[:data], :parent)
        return upgraded_dict
      end

      # we first create the parent list
      parent_dict = upgraded_dict[:data][:parent]
      parents = []
      if haskey(parent_dict, :id)
        push!(parents, parent_dict)
      end
      
      while haskey(parent_dict[:data], :base_ring)
        parent_dict = parent_dict[:data][:base_ring]

        if !haskey(parent_dict, :id)
          break
        end
        
        push!(parents, parent_dict)
        if haskey(parent_dict[:data], :def_pol)
          parent_dict = parent_dict[:data][:def_pol][:data][:parent]
          push!(parents, parent_dict)
        end
      end
      
      refs = Dict()
      upgraded_parents = []

      # using the parent list we make the refs to be attached
      # at the root of the file
      for parent in parents
        id = parent[:id]
        
        pushfirst!(upgraded_parents, Dict(:id => id, :type => "#backref"))
        if haskey(refs, Symbol(id))
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
        refs[Symbol(id)] = parent
      end
      
      upgraded_dict[:data][:parents] = upgraded_parents
      terms = []
      exponent = 0

      # we convert the vector of coefficients
      # to the new terms format
      for coeff in upgraded_dict[:data][:coeffs][:data][:vector]
        if !isa(coeff, Dict)
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
      upgraded_dict[:data] = Dict(:terms => terms, 
                                  :parents => upgraded_parents)
      merge!(s.refs, refs)
    end

    if !isempty(s.refs)
      upgraded_dict[:refs] = s.refs
    end

    return upgraded_dict
  end
))
