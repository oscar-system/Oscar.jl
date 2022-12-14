################################################################################
# Upgrade Summary
# Defines an OscarBasictype, and any Vector{OscarBasictype} will not have
# backrefs in any entry. Instead all entries of the array will have their
# representative only and the entry_type is given for all entries once on
# the save level as the :vector key
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

