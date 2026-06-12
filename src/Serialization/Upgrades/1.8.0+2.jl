# Upgrades are to change the data format of
# `FPGroup`, `SubFPGroup`, `PcGroup`, `SubPcGroup`
# such that the underlying GAP objects are no longer part of
# the serialized data.
# After the upgrade, some serialized GAP objects may be no longer needed
# (and no longer reachable) for the deserialization; the upgrade code does
# not remove these unnecessary data.
#
push!(upgrade_scripts_set, UpgradeScript(
  v"1.8.0+2",
  function upgrade_1_8_0_2(s::UpgradeState, dict::AbstractDict{Symbol, Any})
    haskey(dict, :_refs) || return dict
    refs = dict[:_refs]

    # Record the connections from Oscar groups to GAP groups.
    pointers = Dict{String, String}()  # from GAP group to the corr. OSCAR group
    for key in keys(refs)
      ref = refs[key]
      if ref[:_type] isa AbstractDict && haskey(dict[:_type], :name) &&
         ref[:_type][:name] in ["FPGroup", "SubFPGroup", "PcGroup", "SubPcGroup"]
        pointers[ref[:_type][:params]] = String(key)
      end
    end

    # Rewrite the refs between GAP groups to refs between Oscar groups.
    for (key, ref) in refs
      if ref[:_type] isa AbstractDict && haskey(dict[:_type], :name)
        name = ref[:_type][:name]
        if name == "FPGroup"
          gapref = refs[Symbol(ref[:_type][:params])]
          data = gapref[:data]

          if haskey(data, :names)
            # full free group:
            # transfer the names and the filter identifier.
            ref[:data] = Dict{Symbol, Any}(
                           :names => data[:names],
                           :rep => Oscar.Serialization._word_filters[
                                     Symbol(data[:wfilt])])

            # no ref is needed anymore
            ref[:_type] = "FPGroup"
          else
            # full f.p. group:
            # transfer the relators
            @assert haskey(data, :relators)
            ref[:data] = Dict{Symbol, Any}(
                           :relators => data[:relators])

            # lift the ref to the full group from GAP to OSCAR
            pointer = gapref[:_type][:params]
            ref[:_type][:params] = pointers[pointer]
          end
        elseif name == "PcGroup"
          gapref = refs[Symbol(ref[:_type][:params])]
          data = gapref[:data]

          # transfer the defining data of the collector
          @assert haskey(data, :relord)
          ref[:data] = Dict{Symbol, Any}(
                         :relord => data[:relord],
                         :power_rels => data[:power_rels],
                         :comm_rels => data[:comm_rels])

          # no ref is needed anymore
          ref[:_type] = "PcGroup"
        elseif name == "SubPcGroup" || name == "SubFPGroup"
          gapref = refs[Symbol(ref[:_type][:params])]
          data = gapref[:data]

          # transfer generators
          ref[:data] = Dict{Symbol, Any}(
                         :gens => data[:gens])

          # lift the ref to the full group from GAP to OSCAR
          pointer = gapref[:_type][:params]
          ref[:_type][:params] = pointers[pointer]
        end
      end
    end
    return dict
  end
))

