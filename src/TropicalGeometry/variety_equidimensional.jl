################################################################################
#
#  Tropicalization of equi-dimensional ideals
#
#  WARNING: assumes without test that `I` is equi-dimensional
#
################################################################################

function tropical_variety_equidimensional(I::MPolyIdeal, nu::TropicalSemiringMap{QQField,Nothing,<:Union{typeof(min),typeof(max)}}; weighted_polyhedral_complex_only::Bool=false)
    return tropical_variety_equidimensional_singular(I,nu; weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
end

# trivial valuation case
function tropical_variety_equidimensional_singular(I::MPolyIdeal, nu::TropicalSemiringMap{QQField,Nothing,<:Union{typeof(min),typeof(max)}}; weighted_polyhedral_complex_only::Bool=false)
    R = base_ring(I)
    singularCommand = join(["ring r=0,("*join(string.(symbols(R)),",")*"),dp;",
                            "ideal I = "*join(string.(gens(I)), ",")*";",
                            "if (!defined(tropicalVariety)) { LIB \"tropical.lib\"; };",
                            "fan TropI = tropicalVariety(I);",
                            "string TropIString = string(TropI);"])
    Singular.call_interpreter(singularCommand)
    TropIString = Singular.lookup_library_symbol("Top", "TropIString")
    Sigma = gfan_fan_string_to_oscar_complex(TropIString,convention(nu)==max,false)
    TropI = compute_weights_and_construct_tropical_variety(Sigma,I,nu)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropI,:algebraic_ideal,I)
        set_attribute!(TropI,:tropical_semiring_map,nu)
    end
    return TropI
end

function tropical_variety_equidimensional_singular(I::MPolyIdeal, nu::TropicalSemiringMap{QQField,ZZRingElem,<:Union{typeof(min),typeof(max)}}; weighted_polyhedral_complex_only::Bool=false)
    R = base_ring(I)
    singularCommand = join(["ring r=0,("*join(string.(symbols(R)),",")*"),dp;",
                            "ideal I = "*join(string.(gens(I)), ",")*";",
                            "if (!defined(tropicalVariety)) { LIB \"tropical.lib\"; };",
                            "fan TropI = tropicalVariety(I,number("*string(uniformizer(nu))*"));",
                            "string TropIString = string(TropI);"])
    Singular.call_interpreter(singularCommand)
    TropIString = Singular.lookup_library_symbol("Top", "TropIString")
    Sigma = gfan_fan_string_to_oscar_complex(TropIString,convention(nu)==max,true)
    TropI = compute_weights_and_construct_tropical_variety(Sigma,I,nu)
    if !weighted_polyhedral_complex_only
        set_attribute!(TropI,:algebraic_ideal,I)
        set_attribute!(TropI,:tropical_semiring_map,nu)
    end
    return TropI
end


# Takes the string of a gfan fan (= Singular fan) and converts it to an OSCAR polyhedral complex
# if negateFan==true, return the negative of the fan (because Singular is max-convention only, and OSCAR is flexible)
# if dehomogenizeFan==true, the gfan fan is a homogenized polyhedral complex, and we need to dehomogenize it
function gfan_fan_string_to_oscar_complex(input_string::String, negateFan::Bool=false, dehomogenizeFan::Bool=false)

    # Extracting the RAYS, LINEALITY_SPACE and MAXIMAL_CONES sections
    stringsParsed = Vector{SubString{String}}[]
    for regexp in [r"RAYS\n([\s\S]*?)\nN_RAYS", r"LINEALITY_SPACE\n([\s\S]*?)\nORTH_LINEALITY_SPACE", r"MAXIMAL_CONES\n([\s\S]*)"]
        sectionOfInterest = match(regexp, input_string).captures[1]
        linesOfInterest = split(sectionOfInterest, "\n")
        linesOfInterestFiltered = [split(line, r"\s+#")[1] for line in linesOfInterest if !isempty(line)]
        push!(stringsParsed, linesOfInterestFiltered)
    end

    # Convert Rays and ORTH_LINEALITY_SPACE to matrices
    # and negate if necessary
    rayGenerators = matrix(QQ,[parse.(Int, split(line)) for line in stringsParsed[1]])
    linealityGenerators = matrix(QQ,[parse.(Int, split(line)) for line in stringsParsed[2]])
    if negateFan
        rayGenerators *= -1
    end

    # Convert MAXIMAL_CONES to a Vector{Vector{Int}}
    coneIncidences = [parse.(Int, split(replace(line, r"[{}]" => ""), r"\s+")) .+ 1 for line in stringsParsed[3]]

    if dehomogenizeFan
        # if the singular fan is a homogenized polyhedral complex,
        # identify which rays have a non-zero first entry and actually represent vertices
        # then strip the first coordinate of the rays and lineality generators
        rayIndices = findall(iszero, rayGenerators[:,1])
        rayGenerators = rayGenerators[:,2:end]
        linealityGenerators = linealityGenerators[:,2:end]
        # in some cases, the first unit vector can be a lineality generator
        # hence we need to filter out potential zero rows from linealityGenerators
        linealityGenerators = linealityGenerators[findall(i->!iszero(linealityGenerators[i,:]),1:nrows(linealityGenerators)),:]
        return polyhedral_complex(IncidenceMatrix(coneIncidences), rayGenerators, rayIndices, linealityGenerators)
    else
        # if the singular fan is an honest polyhedral fan
        # we need to add the origin as a vertex and add it to all polyhedra
        rayIndices = collect(1:nrows(rayGenerators))
        rayGenerators = vcat(rayGenerators, zero_matrix(QQ,1,ncols(rayGenerators)))
        originIndex = nrows(rayGenerators)
        coneIncidences = [ vcat(incidence,[originIndex]) for incidence in coneIncidences ]
        return polyhedral_complex(IncidenceMatrix(coneIncidences), rayGenerators, rayIndices, linealityGenerators)
    end

end


# this function takes the tropical variety as a polyhedral complex,
# computes the weights of the tropical variety, and constructs the tropical variety
function compute_weights_and_construct_tropical_variety(Sigma::PolyhedralComplex,::MPolyIdeal,nu::TropicalSemiringMap)
    # TODO: compute multiplicities
    mults = ones(ZZRingElem,n_maximal_polyhedra(Sigma))
    TropI = TropicalVariety{typeof(convention(nu)),true}(Sigma,mults)
    return TropI
end
