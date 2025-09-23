# Upgrade scripts contain code to upgrade files that were serialized with an older
# format to the next file format version. Thus to upgrade a file to
# the latest format, one needs to apply all possible upgrade scripts in
# successive order as determined by their version until finally the file has
# been upgraded to the latest format.

# Warning: The object is not saved to the new format, that is left to the user

@doc raw"""
    UpgradeScript(version::VersionNumber, script::Function)

Any upgrade scripts should be created using the `UpgradeScript`
constructor and then pushed to the `upgrade_scripts_set`. The
name of the function is not particularly important however one
should be assigned so that it can be used for recursion if
necessary for upgrade. The upgrade script should be independent
from any Oscar code and should rely solely on `julia` core
functionality so to avoid conflicts with any future Oscar code
deprecations.

# Examples
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

mutable struct UpgradeState
  id_to_dict::Dict{Symbol, Any}
  nested_level::Int
end

function UpgradeState()
  return UpgradeState(Dict{Symbol, Any}(), 0)
end

(u_s::UpgradeScript)(s::UpgradeState,
                     dict::Dict{Symbol}) = script(u_s)(s, dict)

# The list of all available upgrade scripts
const upgrade_scripts_set = Set{UpgradeScript}()

"""
    upgrade_data(upgrade::Function, s::UpgradeState, dict::Dict)

`upgrade_data` is a helper function that provides functionality for
recursing on the tree structure. Used for upgrades up to 0.12.2.
"""
function upgrade_data(upgrade::Function, s::UpgradeState, dict::Dict)
  s.nested_level += 1
  # file comes from polymake
  haskey(dict, :_ns) && haskey(dict[:_ns], :polymake) && return dict
  
  upgraded_dict = Dict{Symbol, Any}()
  for (key, dict_value) in dict
    if dict_value isa String || dict_value isa Int64 || dict_value isa Bool
      upgraded_dict[key] = dict_value
    elseif dict_value isa Dict
      s.nested_level += 1
      upgraded_dict[key] = upgrade(s, dict_value)
      s.nested_level -= 1
    else  # not a string or a dictionary, so must be a vector
      new_value = []
      for v in dict_value
        if v isa String
          push!(new_value, v)
        else
          s.nested_level += 1
          push!(new_value, upgrade(s, v))
          s.nested_level -= 1
        end
      end
      upgraded_dict[key] = new_value
    end
  end
  s.nested_level -= 1
  return upgraded_dict
end

"""
    rename_types(dict::Dict, renamings::Dict{String, String})

Provides functionality for recursing on the tree structure of `dict`
and replace all type names that occur as keys in renamings with the values. 
"""
function rename_types(dict::Dict, renamings::Dict{String, String})
  function upgrade_type(d::String)
    return get(renamings, d, d)
  end

  function upgrade_type(v::Vector)
    return map(upgrade_type, v)
  end
  
  function upgrade_type(d::Dict)
    upg_d = d

    if haskey(d, :name)
      upg_d[:name] = get(renamings, d[:name], d[:name])
    else
      upg_d[:_type] = get(renamings, d[:_type], d[:_type])
      return upg_d
    end
    
    if d[:params] isa Dict
      if haskey(d[:params], :_type)
        upg_d[:params][:_type] = upgrade_type(d[:params][:_type])
      else
        for (k, v) in d[:params]
          upg_d[:params][k] = upgrade_type(d[:params][k])
        end
      end
    elseif d[:params] isa Vector
      upg_d[:params] = upgrade_type(d[:params])
    end
    return upg_d
  end

  if haskey(dict, :_type)
    dict[:_type] = upgrade_type(dict[:_type])
  end
  return dict
end


function upgrade_containers(upgrade::Function, s::UpgradeState, dict::Dict)
  # all containers have a Dict for their type description
  # with a name and a params key
  dict[:_type] isa String && return dict
  if dict[:data] isa Vector{String}
    ref_entry = get(s.id_to_dict, Symbol(dict[:data][1]), nothing)
    if !isnothing(ref_entry)
      ref_entry = upgrade(s, ref_entry)
      dict[:_type][:params] = ref_entry[:_type]
    end
  end
  
  if dict[:_type][:name] in ["Vector", "Set", "Matrix", "MultiDimArray"]
    subtype = dict[:_type][:params]
    upgraded_entries = Dict[]
    upgraded_entry = nothing
    for entry in dict[:data]
      upgraded_entry = upgrade(s, Dict(:_type => subtype, :data => entry))
      push!(upgraded_entries, upgraded_entry)
    end
    if !isnothing(upgraded_entry)
      dict[:_type][:params] = upgraded_entry[:_type]
    end
    dict[:data] = [u_e[:data] for u_e in upgraded_entries]
  elseif dict[:_type][:name] == "NamedTuple"
    upgraded_entries = Dict[]
    upgraded_entry = nothing
    for (type, entry) in zip(dict[:_type][:params], dict[:data])
      upgraded_entry = upgrade(s, Dict(:_type => type, :data => entry))
      push!(upgraded_entries, upgraded_entry)
    end
    dict[:_type][:params] = [u_e[:_type] for u_e in upgraded_entries]
    dict[:data] = [u_e[:data] for u_e in upgraded_entries]

  elseif dict[:_type][:name] == "NamedTuple"
    upgraded_entries = Dict[]
    upgraded_entry = nothing
    for (type, entry) in zip(dict[:_type][:params][:tuple_params], dict[:data])
      upgraded_entry = upgrade(s, Dict(:_type => type, :data => entry))
      push!(upgraded_entries, upgraded_entry)
    end
    dict[:_type][:params][:tuple_params] = [u_e[:_type] for u_e in upgraded_entries]
    dict[:data] = [u_e[:data] for u_e in upgraded_entries]
  elseif dict[:_type][:name] == "Dict"
    key_params = dict[:_type][:params][:key_params]

    if haskey(dict[:_type][:params], :value_params)
      value_params = dict[:_type][:params][:value_params]
      upgraded_entry = nothing
      upgraded_pairs = Tuple[]
      for (k, v) in dict[:data]
        upgraded_v = upgrade(s, Dict(:_type => value_params, :data => v))
        upgraded_k = upgrade(s, Dict(:_type => key_params, :data => k))
        push!(upgraded_pairs, (upgraded_k, upgraded_v))
      end

      if key_params in ["Symbol", "Int", "String"]
        dict[:data] = Dict()
        for (upgraded_k, upgraded_v) in upgraded_pairs
          dict[:data][upgraded_k[:data]] = upgraded_v[:data]
        end
      else
        dict[:data] = map(x -> [x[1][:data], x[2][:data]], upgraded_pairs)
      end

      if !isempty(upgraded_pairs)
        first_pair = first(upgraded_pairs)
        dict[:_type][:params][:key_params] = first_pair[1][:_type]
        dict[:_type][:params][:value_params] = first_pair[2][:_type]
      end
    else
      for entry in dict[:data]
        upgraded_entry = upgrade(s, Dict(:_type => dict[:_type][:params][entry.first],
                                         :data => entry.second))
        dict[:data][entry.first] = upgraded_entry[:data]
      end
    end
  end
end

# following order is the order in which upgrades will happen

include("0.11.3.jl")
include("0.12.0.jl")
include("0.12.2.jl")
include("0.13.0.jl")
include("0.15.0.jl")
include("1.1.0.jl")
include("1.2.0.jl")
include("1.3.0.jl")
include("1.4.0.jl")
include("1.6.0.jl")

const upgrade_scripts = collect(upgrade_scripts_set)
sort!(upgrade_scripts; by=version)

################################################################################
# Loading with upgrade checks on dict
const backref_sym = Symbol("#backref")

@doc raw"""
    upgrade(format_version::VersionNumber, dict::Dict)

Find the first version where an upgrade can be applied and then incrementally
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
      @debug("upgrading serialized data....",
             maxlog=1)

      upgrade_state = UpgradeState()
      # upgrading large files needs a work around since the new load
      # uses JSON3 which is read only
      upgraded_dict = upgrade_script(upgrade_state, upgraded_dict)
      if script_version > v"0.13.0"
        if haskey(upgraded_dict, :_refs)
          upgraded_refs = Dict{Symbol, Any}()
          for (k, v) in upgraded_dict[:_refs]
            upgraded_refs[k] = upgrade_script(upgrade_state, v)
          end
          upgraded_dict[:_refs] = upgraded_refs
        end
      end
    end
  end
  upgraded_dict[:_ns] = get_oscar_serialization_version()
  return upgraded_dict
end
