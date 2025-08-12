################################################################################
# Upgrade Summary
# 
# `FreeAssAlgebra`: We renamed `FreeAssAlgebra` to `FreeAssociativeAlgebra`, 
# as well as `FreeAssAlgElem` to `FreeAssociativeAlgebraElem` and this upgrade script takes care of this renaming.
  
push!(upgrade_scripts_set, UpgradeScript(
  v"1.2.0",
  function upgrade_1_2_0(s::UpgradeState, dict::Dict)
    renamings = Dict{String,String}([
      ("FreeAssAlgebra", "FreeAssociativeAlgebra"),
      ("FreeAssAlgElem", "FreeAssociativeAlgebraElem"),
      ("FreeAssAlgIdeal", "FreeAssociativeAlgebraIdeal"),
    ])

    upgraded_dict = upgrade_types(dict, renamings)

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
