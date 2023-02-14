# Upgrade scripts contain code to upgrade files that were serialized with an older
# format to the next file format version. Thus to upgrade a file to
# the latest format, one needs to apply all possible upgrade scripts in
# successive order as determined by their version until finally the file has
# been upgraded to the lastest format.

struct UpgradeScript
    version::VersionNumber # version to be upgraded to
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

include("0.11.3.jl")
include("0.11.4.jl")

sort!(upgrade_scripts; by=version)

################################################################################
# Loading with upgrade checks on dict

# Finds the first version where an upgrade can be applied and then incrementally
# upgrades to the current version
function upgrade(dict::Dict{Symbol, Any}, dict_version::VersionNumber)
    upgraded_dict = dict
    for upgrade_script in upgrade_scripts
        script_version = version(upgrade_script)

        if dict_version < script_version
            # TODO: use a macro from Hecke that will allow user to surpress
            # such a message
            @info("upgrading serialized data....", maxlog=1)
            s = DeserializerState()
            upgraded_dict = upgrade_script(s, upgraded_dict)
        end
    end

    upgraded_dict[:_ns] = oscarSerializationVersion
    return upgraded_dict
end
