# Upgrade scripts contain code to upgrade serialization files from their
# current format to the next file format version. Thus to upgrade a file to
# the latest format, one needs to apply all possible upgrade scripts in
# successive order as determined by their version until finally the file has
# been upgraded to the lastest format.
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

# The list of all available upgrade scripts
const upgrade_scripts = Vector{UpgradeScript}()

include("0.11.2.jl")
include("0.11.3.jl")

sort!(upgrade_scripts; by=version)
