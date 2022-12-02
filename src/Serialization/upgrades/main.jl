struct UpgradeScript
    version::VersionNumber
    script::Function
end

function version(upgrade_script::UpgradeScript)
    return upgrade_script.version
end

const upgrade_scripts = Vector{UpgradeScript}()

include("0.11.2.jl")

sort!(upgrade_scripts; by=version)
