################################################################################
# Upgrade Summary
# Types of the form `Vector{T}` are only serialized with backref when
# `is_basic_serialization_type(T) == false`. If
# `is_basic_serialization_type(T) == true` the type `Vector{T}` does not
# serialize `T` in the entries, but instead as a part of the `Vector` serialization.
# Instances of the form `Tuple{T}` with `is_basic_serialization_type(T) == true` 
# are no longer serialized with backref, however these types were already storing
# there types outside of the JSON array of entries

push!(upgrade_scripts, UpgradeScript(
    v"0.11.3", # version this script upgrades to
    function (s::DeserializerState, dict::Dict)
        # file comes from polymake
        if haskey(dict, :_ns)
            haskey(dict[:_ns], :polymake) && return dict
        end

        # moves down tree to point where type exists in dict
        # since we are only doing updates based on certain types
        if !haskey(dict, :type)
            upgraded_dict = Dict{Symbol, Any}()
            for (key, value) in dict
                if value isa String
                    upgraded_dict[key] = value
                else
                    # recursive call unnamed function with var"#self"
                    upgraded_dict[key] = var"#self#"(s, value)
                end
            end
            return upgraded_dict
        end

        if dict[:type] == string(backref_sym)
            backref = s.objs[UUID(dict[:id])]
            backref_type = decodeType(backref[:type])

            # recursive call unnamed function with var"#self"
            # replace backrefs for types that no longer support backrefs
            is_basic_serialization_type(backref_type) && return var"#self#"(s, backref)
            return dict
        end

        # this is one of the base cases
        !haskey(dict, :data) && return dict

        if dict[:type] == "Vector"
            upgraded_vector = []
            entry_type = nothing

            for entry in dict[:data][:vector]
                # recursive call unnamed function with var"#self"
                result = var"#self#"(s, entry)

                # store values that aren't backrefs in state
                if entry[:type] != string(backref_sym)
                    s.objs[UUID(entry[:id])] = entry
                    entry_type = entry[:type]
                end

                push!(upgraded_vector, result)
            end
            
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

        if is_basic_serialization_type(U)
            if U == fmpq
                num = dict[:data][:num][:data]
                den = dict[:data][:den][:data]
                return "$num//$den"
            else
                return dict[:data]
            end
        end

        upgraded_data = var"#self#"(s, dict[:data])

        return Dict(
            :type => dict[:type],
            :data => upgraded_data,
            :id => dict[:id]
        )
    end
))

