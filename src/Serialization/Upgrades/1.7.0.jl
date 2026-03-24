################################################################################
# Upgrade Summary
# 
# We renamed `MatrixGroup` to `MatGroup`, as well as `MatrixGroupElem` to `MatGroupElem`
# in https://github.com/oscar-system/Oscar.jl/pull/5704 and this upgrade script takes care of this renaming.
  
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
