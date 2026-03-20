################################################################################
# Upgrade Summary
# 
# `FreeAssAlgebra`: We renamed `FreeAssAlgebra` to `FreeAssociativeAlgebra`, 
# as well as `FreeAssAlgElem` to `FreeAssociativeAlgebraElem` and this upgrade script takes care of this renaming.
  
push!(upgrade_scripts_set, UpgradeScript(
  v"1.7.0",
  function upgrade_1_7_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    renamings = Dict{String,String}([
      ("MatrixGroup", "MatGroup"),
      ("MatrixGroupElem", "MatGroupElem"),
    ])

    upgraded_dict = rename_types(dict, renamings)

    return upgraded_dict
  end
))
