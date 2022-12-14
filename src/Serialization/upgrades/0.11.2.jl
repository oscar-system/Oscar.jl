################################################################################
# Upgrade Summary
# Instances of the types in this union are not serialized with backref,
# and when serialized as entry of a `Vector{T} where T <: OscarBasicType` then
# the type is not serialized in the entry, but only in the vector.
# Other containers types are also affected by this upgrade, mainly
# container types that use Vector serialization. 
# Instances of the form `Tuple{T} where T <: OscarBasicType` are also no longer
# serialized with backref, however these types were already storing there types
# outside of the JSON array of entries

upgrade_script = UpgradeScript(
    v"0.11.2",
    function (s::DeserializerState, dict::Dict) 
        if haskey(dict, :_ns)
            haskey(dict[:_ns], :polymake) && return dict
        end
    
        if !haskey(dict, :type)
            upgraded_dict = Dict{Symbol, Any}()
            for (key, value) in dict
                if value isa String
                    upgraded_dict[key] = value
                else
                    upgraded_dict[key] = var"#self#"(s, value)
                end
            end
            return upgraded_dict
        end

        !haskey(dict, :data) && return dict
        
        if dict[:type] == string(backref_sym)
            backref = s.objs[UUID(dict[:id])]
            backref_type = decodeType(backref[:type])

            backref_type <: OscarBasicType && return var"#self#"(s, backref)
            return dict
        end

        if dict[:type] == "Vector"
            upgraded_vector = []
            entry_type = nothing
            
            for entry in dict[:data][:vector]
                entry_type = entry[:type]
                push!(upgraded_vector, var"#self#"(s, entry))
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

        if U <: OscarBasicType
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
)

push!(upgrade_scripts, upgrade_script)

