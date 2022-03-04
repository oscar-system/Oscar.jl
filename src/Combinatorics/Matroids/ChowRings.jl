export chow_ring, select

function chow_ring(M::Matroid, extended=false)
    Flats = flats(M)
    number_flats = length(Flats)
    n = length(M.groundset)
    if !is_loopless(M)
        error("Matroid has loops")
    end
    if number_flats<2
        error("matroid has to few flats")
    end
    proper_flats = Flats[2:number_flats-1]
    var_names = [replace(string("x_",S), "["=>"{", "]"=>"}", ", "=>",") for S in proper_flats] # create variable names, indexed by sets  
    if extended 
        var_names = [var_names; [replace(string("h_",S), "["=>"{", "]"=>"}", ", "=>",") for S in [proper_flats;[Flats[number_flats]]]]]
    end

    ring, vars = GradedPolynomialRing(QQ, var_names)
    #ring, vars = PolynomialRing(ZZ, var_names) $which ring do we want here
    I = linear_relations(ring, proper_flats, vars, M)
    J = quadratic_relations(ring, proper_flats, vars)
    Ex = []
    if extended
        Ex = relations_extended_ring(ring, proper_flats, vars)
    end
    chow_modulus = ideal(ring, vcat(I, J, Ex))
    #chow_ring = ResidueRing(ring, chow_modulus)# is this a better choice?
    chow_ring, projection = quo(ring, chow_modulus)
    return chow_ring
end

function linear_relations(ring, proper_flats, vars, M)
    alpha = 0
    relations = []
    for i in M.groundset
        poly = ring(0)
        for index in findall(issubset([i],F) for F in proper_flats)
            poly+= vars[index]
        end
        if i==M.groundset[1]
            alpha = poly
        else
            push!(relations, alpha-poly)
        end
    end
    return relations
end

function quadratic_relations(ring, proper_flats, indeterminates)
    relations = Vector()
    s = size(proper_flats)[1]
    for i in 1:s
        F = proper_flats[i]
        for j in 1:i-1
            G = proper_flats[j]
            if !issubset(F,G) && !issubset(G,F) 
                push!(relations, indeterminates[i]*indeterminates[j])
            end
        end
    end
    return relations
end

function relations_extended_ring(ring, proper_flats, vars)
    relations = Vector()
    s = size(proper_flats)[1]

    # h_E = alpha = -x_E
    poly = ring(0)
    for index in findall(issubset([1],F) for F in proper_flats)
        poly+= vars[index]
    end
    push!(relations, poly-vars[2*s+1])

    # h_F = h_E - sum x_G where G is a proper flat containing F
    for i in 1:s
        F = proper_flats[i]
        poly = ring(0)
        for index in findall(issubset(F,G) for G in proper_flats)
            poly+= vars[index]
        end
        push!(relations, poly+vars[s+i]-vars[2*s+1]) #add h_F and subtract h_E
    end
    return relations
end

function select(include::Union{AbstractVector,Set},exclude::Union{AbstractVector,Set},set::Union{AbstractVector,Set})
    all = []
    for e in set
        all = union(all,e)
    end
    return findall(s->issubset(include,s)&&issubset(s,setdiff(all,exclude)),set);
end
