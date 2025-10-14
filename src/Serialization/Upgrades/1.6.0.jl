push!(upgrade_scripts_set, UpgradeScript(
  v"1.6.0",
  function upgrade_1_6_0(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    # recurse upgrade on containers
    upgrade_containers(upgrade_1_6_0, s, dict)

    # Upgrades 
    if dict[:_type] == "PhylogeneticTree"
      dict[:_type] = Dict{Symbol, Any}(
        :name => "PhylogeneticTree",
        :params => Dict{Symbol, Any}(
          :_type => dict[:data][:_type] == "graph::PhylogeneticTree<Float>" ? "Floats" : "QQField"
        )
      )
      n_vertices = length(dict[:data][:LABELS])
      dict[:data] = Dict{Symbol, Any}(
        :pm_tree => dict[:data],
        :vertex_perm => string.(collect(1:n_vertices))
      )
    end
    return dict
  end
))
