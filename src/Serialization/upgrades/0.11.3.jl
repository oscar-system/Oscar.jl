################################################################################
# Upgrade Summary
# added for proof of concept to be removed later

upgrade_script = UpgradeScript(
    v"0.11.3",
    function (s::DeserializerState, dict::Dict)
        println("Nothing to upgrade")
        return dict
    end
)

push!(upgrade_scripts, upgrade_script)


