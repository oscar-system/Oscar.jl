###################################
# (1) Call cohomCalg to compute line bundle cohomolgy
###################################

function command_string(v::AbstractNormalToricVariety, c::Vector{fmpz})
    # Initialize a list of strings, which we will eventually join to form the command string
    string_list = Vector{String}()
    
    # Define helper function
    joincomma(list) = join([string(x) for x in list], ",")
    
    # Add information about grading of Cox ring to string_list
    divisors = gens(torusinvariant_divisor_group(v))
    for i in 1:length(divisors)
        tmp = joincomma(map_from_weil_divisors_to_class_group(v)(divisors[i]).coeff)
        push!(string_list, "vertex x$i|GLSM:($(tmp))")
    end
    
    # Add information about the Stanley-Reisner ideal to string_list
    current_coordinate_names = [string(x) for x in Hecke.gens(cox_ring(v))]
    new_coordinate_names = ["x$i" for i = 1:length(current_coordinate_names)]
    set_coordinate_names(v, new_coordinate_names)
    generators = [string(g) for g in gens(stanley_reisner_ideal(v))]
    push!(string_list, "srideal [" * joincomma(generators) * "]")
    set_coordinate_names(v, current_coordinate_names)
    
    # Add line bundle information to string_list
    push!(string_list, "ambientcohom O(" * joincomma(c) * ");")
    
    # Join and return
    return join(string_list, ";")
end
command_string(v::AbstractNormalToricVariety) = command_string(v, [fmpz(0) for i in 1:rank(picard_group(v))])



###################################
# (2) Compute the denominators of the rationoms
###################################

function contributing_denominators(variety::AbstractNormalToricVariety)
    
    # execute cohomCalg
    out = Pipe()
    err = Pipe()
    process = run(pipeline(ignorestatus(`$(cohomCalg_jll.cohomcalg()) --verbose2 --in=$(command_string(variety))`), stdout=out, stderr=err))
    close(out.in)
    close(err.in)
    
    # error during execution?
    stderr = read(err, String)
    code = process.exitcode
    if code != 0
        error("cohomCalg encountered the error " * stderr)
    end
    
    # process result
    stdout = read(out, String)
    start = findfirst("Final list of contributing monomials with factors:", stdout)[1]
    diff = findfirst("Verbose Level 1:", stdout)[1]
    output_string_reduced = [strip(s) for s in split(SubString(stdout, start, diff),"\n")]
    output_string_reduced = output_string_reduced[4:length(output_string_reduced)-3]
    
    # ambiguous monomial contributions found during execution?
    if last(output_string_reduced) != "There are no ambiguous contribution monomials to the variety."
        throw(ArgumentError("cohomCalg experienced ambiguous monomial contributions."))
    end
    
    # Monomials can in principal come with contributing factor different from 1.
    # In this case, the word "factor ..." occurs. Here, we do not care about these factors.
    # To process the result efficiently, we remove this information.
    for i in 1:length(output_string_reduced)
        if occursin("factor", output_string_reduced[i])
            output_string_reduced[i] = strip(SubString(output_string_reduced[i], 1, findfirst("factor", output_string_reduced[i])[1]-1))
        end
    end
    
    # at what positions in the output_string does the list of contributions to cohomology group i end and that of i+1 begins?
    hlist = Int[]
    for i in 1:length(output_string_reduced)-1
        if output_string_reduced[i][1] == 'H'
            push!(hlist, i)
        end
    end
    
    # form list of the contributing monomials
    contributing_monomials = [[""] for i in 1:dim(variety)+1]
    for i in 1:dim(variety)+1
        if ((i <= dim(variety)) && (hlist[i]+1 !== hlist[i+1]))
            contributing_monomials[i] = [output_string_reduced[j] for j in hlist[i]+1:hlist[i+1]-1]
        elseif ((i == dim(variety)+1) && (hlist[i]+1 !== length(output_string_reduced)))
            contributing_monomials[i] = [output_string_reduced[j] for j in hlist[i]+1:length(output_string_reduced)-1]
        end
    end
    
    # we execute cohomCalg with homogeneous variable names "x1", "x2" and so on
    # -> translate back into actually chosen variable names
    gens = Hecke.gens(cox_ring(variety))
    for i in 1:length(contributing_monomials)
        for j in 1:length(contributing_monomials[i])
            m = contributing_monomials[i][j]
            present_variables = [occursin("x$k", m) for k in 1:ngens(cox_ring(variety))]
            contributing_monomials[i][j] = string(prod(gens[k]^present_variables[k] for k in 1:ngens(cox_ring(variety))))
        end
    end
    
    # return result
    return contributing_monomials
    
end
export contributing_denominators



######################################
# (3) Turn denominators into shifted cones (-> vanishing sets)
######################################

function turn_denominator_into_polyhedron(variety::AbstractNormalToricVariety, monom::String)
    
    # (1) which variables appear in the monom?
    present_variables = [occursin(coordinate_names(variety)[i], monom) for i in 1:ngens(cox_ring(variety))]
    
    # (2) compute generators of the semigroup
    weights = [k.coeff for k in Oscar._cox_ring_weights(variety)]
    gens = vcat(unique([(-1)^Int(present_variables[i])*weights[i] for i in 1:length(present_variables)]))
    
    # (3) compute offset
    offset = zero(parent(weights[1]))
    for i in 1:length(present_variables)
        if present_variables[i]
            offset -= weights[i]
        end
    end
    return convex_hull(offset, gens)
    
end
export turn_denominator_into_polyhedron
