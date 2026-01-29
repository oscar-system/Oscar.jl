################################################################################
# Upgrade Summary
# 
# `FreeAssAlgebra`: We renamed `FreeAssAlgebra` to `FreeAssociativeAlgebra`, 
# as well as `FreeAssAlgElem` to `FreeAssociativeAlgebraElem` and this upgrade script takes care of this renaming.
  
push!(upgrade_scripts_set, UpgradeScript(
  v"1.2.0",
  function upgrade_1_2_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    renamings = Dict{String,String}([
      ("FreeAssAlgebra", "FreeAssociativeAlgebra"),
      ("FreeAssAlgElem", "FreeAssociativeAlgebraElem"),
      ("FreeAssAlgIdeal", "FreeAssociativeAlgebraIdeal"),
    ])

    upgraded_dict = rename_types(dict, renamings)

    if haskey(dict, :data) && dict[:data] isa AbstractDict
      upgraded_dict[:data] = upgrade_1_2_0(s, dict[:data])
    end

    return upgraded_dict
  end
))
