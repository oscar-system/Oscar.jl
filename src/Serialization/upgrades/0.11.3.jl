################################################################################
# Upgrade Summary
# When `is_basic_serialization_type(T) == true` the type `Vector{T}` the entries
# are no longer serialized as dicts, The type is encoded on the Vector at the same level as
# the entry data. The entries no longer being dictionaries prevents the use of backrefs.


push!(upgrade_scripts, UpgradeScript(
    v"0.11.3", # version this script upgrades to
    function upgrade_0_11_3(s::DeserializerState, dict::Dict)
        # file comes from polymake
        haskey(dict, :_ns) && haskey(dict[:_ns], :polymake) && return dict

        # moves down tree to point where type exists in dict
        # since we are only doing updates based on certain types
        # no :type key implies the dict is data
        if !haskey(dict, :type)
            upgraded_dict = Dict{Symbol, Any}()
            for (key, value) in dict
                if value isa String
                    upgraded_dict[key] = value
                elseif value isa Dict{Symbol, Any}
                    upgraded_dict[key] = upgrade_0_11_3(s, value)
                else
                    values = []
                    for v in value
                        if v isa String
                            push!(values, v)
                        else
                            push!(values, upgrade_0_11_3(s, v))
                        end
                    end
                    upgraded_dict[key] = values
                end
            end
            return upgraded_dict
        end

        if dict[:type] == string(backref_sym)
            backref = s.objs[UUID(dict[:id])]

            # data no longer uses backref and is part of the update
            backref isa String && return backref

            return dict
        end

        # this handles types like FlintRationalField that have no associated data
        !haskey(dict, :data) && return dict

        # apply upgrade to vector entries
        # we only do this here since upgrade_0_11_13 is only defined for dicts
        # we could make this it's own function outside of this script
        if dict[:type] == "Vector"
            upgraded_vector = []
            entry_type = nothing

            for entry in dict[:data][:vector]
                result = upgrade_0_11_3(s, entry)

                # store values in state that are vector entries that aren't backrefs
                if entry[:type] != string(backref_sym)
                    if haskey(entry, :id)
                        s.objs[UUID(entry[:id])] = result
                    end
                    entry_type = entry[:type]
                end

                push!(upgraded_vector, result)
            end

            # no update at this level
            upgraded_vector[1] isa Dict && return Dict(
                :type => "Vector",
                :id => dict[:id],
                :data => Dict(
                    :vector => upgraded_vector
                )
            )

            return Dict(
                :type => "Vector",
                :id => dict[:id],
                :data => Dict(
                    :vector => upgraded_vector,
                    :entry_type => entry_type
                )
            )
        end

        U = decodeType(dict[:type])

        # Upgrades fmpq serialization
        if is_basic_serialization_type(U)
            if U === fmpq
                num = dict[:data][:num][:data]
                den = dict[:data][:den][:data]
                return "$num//$den"
            else
                return dict[:data]
            end
        end

        upgraded_data = upgrade_0_11_3(s, dict[:data])
        upgraded_dict = Dict(
            :type => dict[:type],
            :data => upgraded_data,
            :id => dict[:id]
        )

        # adds updated object in case it is referenced somewhere
        if haskey(dict, :id)
            s.objs[UUID(dict[:id])] = upgraded_dict
        end

        return upgraded_dict
    end
))
