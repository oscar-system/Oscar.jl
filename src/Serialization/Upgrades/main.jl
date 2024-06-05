# Upgrade scripts contain code to upgrade files that were serialized with an older
# format to the next file format version. Thus to upgrade a file to
# the latest format, one needs to apply all possible upgrade scripts in
# successive order as determined by their version until finally the file has
# been upgraded to the latest format.

# Warning: The object is not saved to the new format, that is left to the user

"""
    UpgradeScript(version::VersionNumber, script::Function)

Any upgrade scripts schould be created using the `UpgradeScript`
constructor and then pushed to the `upgrade_scripts_set`. The
name of the function is not particularly important however one
should be assigned so that it can be used for recursion if
necessary for upgrade. The upgrde script should be independent
from any Oscar code and should relie solely on `julia` core
functionality so to avoid conflicts with any future Oscar code
deprecations.

# Example
```
push!(upgrade_scripts_set, UpgradeScript(
  v"0.13.0",
  function upgrade_0_13_0(s::UpgradeState, dict::Dict)
      ...
  end
))
```
"""
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

struct UpgradeState
  id_to_dict::Dict{Symbol, Any}
end

function UpgradeState()
  return UpgradeState(Dict{Symbol, Any}())
end

(u_s::UpgradeScript)(s::UpgradeState,
                     dict::Dict{Symbol, Any}) = script(u_s)(s, dict)

# The list of all available upgrade scripts
upgrade_scripts_set = Set{UpgradeScript}()

"""
    upgrade_data(upgrade::Function, s::UpgradeState, dict::Dict)

`upgrade_data` is a helper function that provides functionality for
recursing on the tree structure. It is independent of any particular
file format version and can be used in any upgrade script.
"""
function upgrade_data(upgrade::Function, s::UpgradeState, dict::Dict)
  # file comes from polymake
  haskey(dict, :_ns) && haskey(dict[:_ns], :polymake) && return dict
  
  upgraded_dict = Dict{Symbol, Any}()
  for (key, dict_value) in dict
    if dict_value isa String || dict_value isa Int64 || dict_value isa Bool
      upgraded_dict[key] = dict_value
    elseif dict_value isa Dict
      upgraded_dict[key] = upgrade(s, dict_value)
    else  # not a string or a dictionary, so must be a vector
      new_value = []
      for v in dict_value
        if v isa String
          push!(new_value, v)
        else
          push!(new_value, upgrade(s, v))
        end
      end
      upgraded_dict[key] = new_value
    end
  end
  return upgraded_dict
end

include("0.11.3.jl")
include("0.12.0.jl")
include("0.12.2.jl")
include("0.13.0.jl")
include("0.15.0.jl")
include("1.1.0.jl")

upgrade_scripts = collect(upgrade_scripts_set)
sort!(upgrade_scripts; by=version)

################################################################################
# Loading with upgrade checks on dict
const backref_sym = Symbol("#backref")

@doc raw"""
    upgrade(format_version::VersionNumber, dict::Dict)

Finds the first version where an upgrade can be applied and then incrementally
upgrades to each intermediate version until the structure of the current version
has been achieved.
"""
function upgrade(format_version::VersionNumber, dict::Dict)
  upgraded_dict = dict
  for upgrade_script in upgrade_scripts
    script_version = version(upgrade_script)
    if format_version < script_version
      # TODO: use a macro from Hecke that will allow user to suppress
      # such a message
      @info("upgrading serialized data....",
            maxlog=1)

      upgrade_state = UpgradeState()
      # upgrading large files needs a work around since the new load
      # uses JSON3 which is read only 
      upgraded_dict = upgrade_script(upgrade_state, upgraded_dict)
    end
  end
  upgraded_dict[:_ns] = get_oscar_serialization_version()
  return upgraded_dict
end
