# Upgrade scripts are used to upgrade file formats from their current format
# to the next chronological file format. I.e. to upgrade to the latest file format
# a file will go through all possible upgrade scripts in succesive order until
# it has been upgraded to the lastest format.
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
include("0.11.3.jl")

sort!(upgrade_scripts; by=version)
