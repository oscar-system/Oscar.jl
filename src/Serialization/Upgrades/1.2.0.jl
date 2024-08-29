################################################################################
# Upgrade Summary
# 
# `FreeAssAlgebra`: We renamed `FreeAssAlgebra` to `FreeAssociativeAlgebra`, 
# as well as `FreeAssAlgElem` to `FreeAssociativeAlgebraElem` and this upgrade script takes care of this renaming.
  
push!(upgrade_scripts_set, UpgradeScript(
  v"1.2.0",
  function upgrade_1_2_0(s::UpgradeState, dict::Dict)
    upgraded_dict = dict

    renamings = Dict{String,String}([
      ("FreeAssAlgebra", "FreeAssociativeAlgebra"),
      ("FreeAssAlgElem", "FreeAssociativeAlgebraElem"),
      ("FreeAssAlgIdeal", "FreeAssociativeAlgebraIdeal"),
    ])

    function upgrade_type(d::String)
      return get(renamings, d, d)
    end
    function upgrade_type(d::Dict)
      upg_d = d
      upg_d[:name] = get(renamings, d[:name], d[:name])
      if d[:params] isa Dict && haskey(d[:params], :_type)
        upg_d[:params][:_type] = upgrade_type(d[:params][:_type])
      elseif d[:params] isa Vector
        upg_d[:params] = [upgrade_type(v) for v in d[:params]]
      end
      return upg_d
    end

    if haskey(dict, :_type)
      upgraded_dict[:_type] = upgrade_type(dict[:_type])
    end
        
    if haskey(dict, :data) && dict[:data] isa Dict
      upgraded_dict[:data] = upgrade_1_2_0(s, dict[:data])
    end

    if haskey(dict, :_refs)
      upgraded_refs = Dict()
      for (k, v) in dict[:_refs]
        upgraded_refs[k] = upgrade_1_2_0(s, v)
      end
      upgraded_dict[:_refs] = upgraded_refs
    end

    return upgraded_dict
  end
))
