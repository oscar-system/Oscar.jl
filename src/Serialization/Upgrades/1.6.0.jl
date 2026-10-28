push!(upgrade_scripts_set, UpgradeScript(
  v"1.6.0",
  function upgrade_1_6_0(s::UpgradeState, dict::AbstractDict{String, Any})
    # recurse upgrade on containers
    upgrade_recursive(upgrade_1_6_0, s, dict)

    # Upgrades 
    if dict["_type"] == "PhylogeneticTree"
      dict["_type"] = Dict{String, Any}(
        "name" => "PhylogeneticTree",
        "params" => Dict{String, Any}(
          "_type" => dict["data"]["_type"] == "graph::PhylogeneticTree<Float>" ? "Floats" : "QQField"
        )
      )
      n_vertices = length(dict["data"]["LABELS"])
      dict["data"] = Dict{String, Any}(
        "pm_tree" => dict["data"],
        "vertex_perm" => string.(collect(1:n_vertices))
      )
    end
    return dict
  end
))
