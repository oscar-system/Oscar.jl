# Upgrades are to change the data format of
# the types `FPGroup`, `SubFPGroup`, `PcGroup`, `SubPcGroup`
# such that the underlying GAP objects are no longer part of
# the serialized data.
# After the upgrade, some serialized GAP objects may be no longer needed
# (and no longer reachable) for the deserialization; the upgrade code does
# not remove these unnecessary data.
# Distinguish two cases:
# 1. The object in question is the only object in the file;
#    In this case, it does not occur in the `_refs` part.
# 2. The objects in question are part of bigger data,
#    such as a (nested) container.
#    In this case, all changes refer either to the `_refs` part
#    or to `params` vectors of `_type` entries.
#    In this case, we run over the `_refs` and adjust their contents.
#    At the same time, we collect the necessary replacements that must be
#    applied to the `_type` entries in the file.
#    Finally, we run the recursion in order to adjust the `params` vectors
#    in question.
#
push!(upgrade_scripts_set, UpgradeScript(
  v"1.8.0+2",
  function upgrade_1_8_0_2(s::UpgradeState, dict::AbstractDict{Symbol, Any})

    # We want to run this upgrade only on the top level,
    # and only if there is a `_refs` entry.
    haskey(dict, :_refs) || return dict

    # We are on the top level, prepare the global data for the recursion.
    # We need the `_refs` and the connections between the Oscar groups
    # and the corresponding GAP groups.

    # Step 1: Deal with one top level object of one of our types.
    refs = dict[:_refs]
    if dict[:_type] isa AbstractDict && haskey(dict[:_type], :name) &&
       dict[:_type][:name] in ["FPGroup", "SubFPGroup", "PcGroup", "SubPcGroup"]
      refs = copy(refs)
      refs[Symbol(dict[:id])] = Dict{Symbol, Any}(:_type => dict[:_type], :data => dict[:data])
      toplevel = true
    else
      toplevel = false
    end

    # Step 2: Adjust the `_refs` part, and collect `params` replacements
    #         (non-recursively).

    # Collect pointers from GAP groups to corresponding OSCAR groups.
    pointers = Dict{String, String}()
    for key in keys(refs)
      ref = refs[key]
      if ref[:_type] isa AbstractDict && haskey(ref[:_type], :name) &&
         ref[:_type][:name] in ["FPGroup", "SubFPGroup", "PcGroup", "SubPcGroup"]
        pointers[ref[:_type][:params]] = String(key)
      end
    end

    # Rewrite the refs between GAP groups to refs between Oscar groups.
    replace_type = Dict()
    pairs = collect(refs)

    # auxiliary functions
    insert_oscar_free_group = function(pointer, name)
      # no Oscar free group was stored
      oscaruuid = string(uuid4())
      gapdata = refs[Symbol(pointer)][:data]
      oscarobj = Dict{Symbol, Any}(
                   :_type => "FPGroup",
                   :data => Dict{Symbol, Any}(
                     :names => gapdata[:names],
                     :rep => Oscar.Serialization._word_filters[
                               Symbol(gapdata[:wfilt])]))

      dict[:_refs][Symbol(oscaruuid)] = oscarobj
      push!(pairs, Symbol(oscaruuid) => oscarobj)
      pointers[pointer] = oscaruuid

      replace_type[Dict{Symbol, Any}(
                     :name => name,
                     :params => oscaruuid)] =
                   Dict{Symbol, Any}(
                     :name => name,
                     :params => oscaruuid)
    end

    insert_oscar_full_group = function(pointer, name)
      # no Oscar full group was stored
      oscaruuid = string(uuid4())
      gapdata = refs[Symbol(pointer)][:data]

      if haskey(gapdata, :relators)
        # subgroup of a f.p. group, and this points to a free group
        oscarobj = Dict{Symbol, Any}(
                     :_type => Dict{Symbol, Any}(
                       :name => "FPGroup",
                       :params => pointer),
                     :data => Dict{Symbol, Any}(
                       :relators => gapdata[:relators]))

        dict[:_refs][Symbol(oscaruuid)] = oscarobj
        push!(pairs, Symbol(oscaruuid) => oscarobj)
        pointers[pointer] = oscaruuid

        replace_type[Dict{Symbol, Any}(
                       :name => name,
                       :params => oscaruuid)] =
                     Dict{Symbol, Any}(
                       :name => name,
                       :params => oscaruuid)

        insert_oscar_free_group(gapdata[:freeGroup], "FPGroup")
      elseif haskey(gapdata, :relord)
        # subgroup of a pc group
        oscarobj = Dict{Symbol, Any}(
                     :_type => "PcGroup",
                     :data => Dict{Symbol, Any}(
                       :relord => gapdata[:relord],
                       :power_rels => gapdata[:power_rels],
                       :comm_rels => gapdata[:comm_rels]))

        dict[:_refs][Symbol(oscaruuid)] = oscarobj
        push!(pairs, Symbol(oscaruuid) => oscarobj)
        pointers[pointer] = oscaruuid

        replace_type[Dict{Symbol, Any}(
                       :name => name,
                       :params => oscaruuid)] =
                     Dict{Symbol, Any}(
                       :name => name,
                       :params => oscaruuid)
      else
        # subgroup of a free group
       insert_oscar_free_group(pointer, "SubFPGroup")
      end
    end

    for (key, ref) in pairs
      if ref[:_type] isa AbstractDict && haskey(ref[:_type], :name)
        name = ref[:_type][:name]
        if name == "FPGroup"
          oldparams = ref[:_type][:params]
          gapref = refs[Symbol(oldparams)]
          data = gapref[:data]

          if haskey(data, :names)
            # full free group:
            # transfer the names and the filter identifier.
            ref[:data] = Dict{Symbol, Any}(
                           :names => data[:names],
                           :rep => Oscar.Serialization._word_filters[
                                     Symbol(data[:wfilt])])

            replace_type[Dict{Symbol, Any}(:name => name,
                              :params => oldparams)] = name

            # no ref is needed anymore
            ref[:_type] = "FPGroup"
          elseif haskey(gapref, :_type) && gapref[:_type] isa AbstractDict &&
                 haskey(gapref[:_type], :name) && gapref[:_type][:name] == "GapObj"
            # full f.p. group:
            # transfer the relators
            @assert haskey(data, :relators)
            ref[:data] = Dict{Symbol, Any}(
                           :relators => data[:relators])

            # lift the ref to the free group from GAP to OSCAR
            pointer = gapref[:_type][:params]

            if !haskey(pointers, pointer)
              insert_oscar_free_group(pointer, name)
            end

            replace_type[Dict{Symbol, Any}(:name => name,
                              :params => oldparams)] =
                         Dict{Symbol, Any}(:name => name,
                              :params => pointers[pointer])

            ref[:_type][:params] = pointers[pointer]
          end
        elseif name == "PcGroup"
          oldparams = ref[:_type][:params]
          gapref = refs[Symbol(oldparams)]
          data = gapref[:data]

          # transfer the defining data of the collector
          @assert haskey(data, :relord)
          ref[:data] = Dict{Symbol, Any}(
                         :relord => data[:relord],
                         :power_rels => data[:power_rels],
                         :comm_rels => data[:comm_rels])

          replace_type[Dict{Symbol, Any}(:name => name,
                            :params => oldparams)] = name

          # no ref is needed anymore
          ref[:_type] = "PcGroup"
        elseif name == "SubPcGroup" || name == "SubFPGroup"
          oldparams = ref[:_type][:params]
          gapref = refs[Symbol(oldparams)]
          if haskey(gapref, :_type) && gapref[:_type] isa AbstractDict &&
             haskey(gapref[:_type], :name) && gapref[:_type][:name] == "GapObj"

            data = gapref[:data]

            # transfer generators
            ref[:data] = Dict{Symbol, Any}(
                           :gens => data[:gens])

            # lift the ref to the full group from GAP to OSCAR
            pointer = gapref[:_type][:params]
            if !haskey(pointers, pointer)
              insert_oscar_full_group(pointer, name)
            end

            replace_type[Dict{Symbol, Any}(:name => name,
                              :params => oldparams)] =
                         Dict{Symbol, Any}(:name => name,
                              :params => pointers[pointer])

            ref[:_type][:params] = pointers[pointer]
          end
        end
      end
    end

    # Step 3: Adjust `params` recursively,
    #         using the precomputed replacements list.

    if toplevel
      # Replace the top level `:_type`.
      dict[:_type] = replace_type[dict[:_type]]

      # Replace the top level `:data`.
      dict[:data] = refs[Symbol(dict[:id])][:data]
    end

    # Define a function that has access to `replace_type`.
    recursion = function(s::UpgradeState, dict::AbstractDict{Symbol, Any})
      if haskey(dict, :_type) && haskey(replace_type, dict[:_type])
        dict[:_type] = replace_type[dict[:_type]]
      end
      return dict
    end

    # recurse upgrade on containers
    upgrade_recursive(recursion, s, dict)

    return dict
  end
))
