################################################################################
# Upgrade Summary
# added for proof of concept to be removed later

push!(upgrade_scripts, UpgradeScript(
    v"0.11.4",
    function (s::DeserializerState, dict::Dict)
        println("Nothing to upgrade")
        return dict
    end
))
