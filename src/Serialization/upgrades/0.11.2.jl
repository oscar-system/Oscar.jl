
function upgrade(s::DeserializerState, dict::Dict)
    if dict[:type] == string(backref_sym)
        backref = s.objs[UUID(dict[:id])]
        backref_type = decodeType(backref[:type])

        backref_type <: OscarBasicType && return upgrade(s, backref)
        return dict
    end

    if dict[:type] == "Vector"
        upgraded_vector = []
        entry_type = nothing
        
        for entry in dict[:data][:vector]
            entry_type = entry[:type]
            push!(upgraded_vector, upgrade(s, entry))
        end

        upgraded_vector[1] isa Dict && return upgraded_vector
        
        return Dict(
            :vector => upgraded_vector,
            :entry_type => entry_type
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
    
    return upgrade(s, dict[:data])
end


    
