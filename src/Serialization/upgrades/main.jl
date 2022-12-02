struct UpgradeScript
    version::VersionNumber
    script::Function
end

function version(upgrade_script::UpgradeScript)
    return upgrade_script.version
end

function script(upgrade_script::UpgradeScript)
    return upgrade_script.script
end

(u_s::UpgradeScript)(s::DeserializerState,
                     dict::Dict{Symbol, Any}) = script(u_s)(s, dict)

const upgrade_scripts = Vector{UpgradeScript}()

include("0.11.2.jl")

sort!(upgrade_scripts; by=version)
