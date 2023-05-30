include("Helper.jl")
include("GroebnerWalkUtilities.jl")
include("FractalWalkUtilities.jl")
include("GenericWalkUtilities.jl")
include("TranWalkUtilities.jl")


###############################################################
# Implementation of different variants of the Groebner Walk.
# The Groebner Walk is proposed by Collart, Kalkbrener & Mall (1997).
###############################################################

#=
Compute a reduced Groebner basis w.r.t. to a monomial order by converting it using the Groebner Walk.
The Groebner Walk is proposed by Collart, Kalkbrener & Mall (1997).
One can choose a strategy of:
Standard Walk (:standard) computes the Walk like it´s presented in Cox, Little & O´Shea (2005).
Generic Walk (:generic) computes the Walk like it´s presented in Fukuda, Jensen, Lauritzen & Thomas (2005).
Pertubed Walk (:pertubed, with p = degree of the pertubation) computes the Walk like it´s presented in Amrhein, Gloor & Küchlin (1997).
Tran´s Walk (:tran) computes the Walk like it´s presented in Tran (2000).
Fractal Walk (:fractalcombined) computes the Walk like it´s presented in Amrhein & Gloor (1998) with multiple extensions. The target monomial order has to be lex. This version uses the Buchberger Algorithm to skip weight vectors with entries bigger than Int32.
Fractal Walk (:fractal) computes the Walk like it´s presented in Amrhein & Gloor (1998). Pertubes only the target vector.

#Arguments
*`G::Oscar.IdealGens`: ideal one wants to compute a Groebner basis for.
*`startOrder::Symbol=:degrevlex`: monomial order to begin the conversion.
*`targetOrder::Symbol=:lex`: monomial order one wants to compute a Groebner basis for.
*`walktype::Symbol=standard`: strategy of the Groebner Walk. One can choose a strategy of:
    - `standard`: Standard Walk,
    - `pertubed`: Pertubed Walk,
    - `tran`: Tran´s Walk,
    - `generic`: Generic Walk,
    - `fractal`: standard-version of the Fractal Walk,
    - `fractalcombined`: combined Version of the Fractal Walk. Target monomial order needs to be lex,
*`pertubationDegree::Int=2`: pertubationdegree for the Pertubed Walk.
*'infoLevel::Int=0':
    -'0': no printout,
    -'1': intermediate weight vectors,
    -'2': information about the Groebner basis.
=#
function groebnerwalk(
    I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)),
    startOrder::Symbol = :degrevlex,
    targetOrder::Symbol = :lex,
    walktype::Symbol = :standard,
    pertubationDegree::Int = 2,
    infoLevel::Int = 0,
)
    R = change_order(
        base_ring(I),
        ordering_as_matrix(startOrder, nvars(base_ring(I))),
    )
    Gb = groebner_basis(
        ideal(R, [change_ring(x, R) for x in gens(I)]),
        complete_reduction = true,
    )

    return groebnerwalk(
        Gb,
        ordering_as_matrix(startOrder, nvars(R)),
        ordering_as_matrix(targetOrder, nvars(R)),
        walktype,
        pertubationDegree,
        infoLevel,
    )
end

#=
Computes a reduced Groebner basis w.r.t. the monomial order T by converting the reduced Groebner basis G w.r.t. the monomial order S using the Groebner Walk.
One can choose a strategy of:
Standard Walk (:standard) computes the Walk like it´s presented in Cox et al. (2005).
Generic Walk (:generic) computes the Walk like it´s presented in Fukuda et al. (2005).
Pertubed Walk (:pertubed, with p = degree of the pertubation) computes the Walk like it´s presented in Amrhein et al. (1997).
Tran´s Walk (:tran) computes the Walk like it´s presented in Tran (2000).
Fractal Walk (:fractal) computes the Walk like it´s presented in Amrhein & Gloor (1998). Pertubes only the target vector.
Fractal Walk (:fractalcombined) computes the Walk like it´s presented in Amrhein & Gloor (1998) with multiple extensions. The target monomial order has to be lex. This version uses the Buchberger Algorithm to skip weightvectors with entries bigger than Int32.

#Arguments
*`G::Oscar.IdealGens`: Groebner basis to convert to the Groebner basis w.r.t. the target order T. G needs to be a Groebner basis w.r.t. the start order S.
*`S::Matrix{Int}`: start monomial order w.r.t. the Groebner basis G. Note that S has to be a nxn-matrix with rank(S)=n and its first row needs to have positive entries.
*`T::Matrix{Int}`: target monomial order one wants to compute a Groebner basis for. Note that T has to be a nxn-matrix with rank(T)=n and its first row needs to have positive entries.
*`walktype::Symbol=standard`: strategy of the Groebner Walk. One can choose a strategy of:
    - `standard`: Standard Walk (default),
    - `pertubed`: Pertubed Walk,
    - `tran`: Tran´s Walk,
    - `generic`: Generic Walk,
    - `fractal`: standard-version of the Fractal Walk,
    - `fractalcombined`: combined version of the Fractal Walk. The target monomial order needs to be lex,
*`p::Int=2`: pertubationdegree for the pertubed Walk.
*'infoLevel::Int=0':
    -'0': no printout,
    -'1': intermediate weight vectors,
    -'2': information about the Groebner basis.
=#
function groebnerwalk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    walktype::Symbol = :standard,
    p::Int = 2,
    infoLevel::Int = 1,
) 
    if walktype == :standard
        walk = (x) -> standard_walk(x, S, T, infoLevel)
    elseif walktype == :generic
        walk = (x) -> generic_walk(x, S, T, infoLevel)
    elseif walktype == :pertubed
        walk = (x) -> pertubed_walk(x, S, T, p, infoLevel)
    elseif walktype == :fractal
        walk = (x) -> fractal_walk(x, S, T, infoLevel)
    elseif walktype == :fractal_start_order
        walk = (x) -> fractal_walk_start_order(x, S, T, infoLevel)
#    elseif walktype == :fractal_lex
 #       walk = (x) -> fractal_walk_lex(x, S, T, infoLevel)
 #   elseif walktype == :fractal_look_ahead
 #       walk = (x) -> fractal_walk_look_ahead(x, S, T, infoLevel)
    elseif walktype == :tran
        walk = (x) -> tran_walk(x, S, T, infoLevel)
    elseif walktype == :fractal_combined
        walk = (x) -> fractal_walk_combined(x, S, T, infoLevel)
    end

    delete_counter()
    !check_order_M(S, T, G) && throw(
        error(
            "The matrices representing the monomial order have to be nxn-matrices with full rank.",
        ),
    )

   # R = singular_poly_ring(G.Ox, matrix_ordering(G.Ox, S))
   # i = Singular.Ideal(R, [R(x) for x = G])
    Gb = walk(G)

    if infoLevel >= 1
        println("Cones crossed: ", delete_counter())
    end
    return Oscar.IdealGens(gens(Gb), matrix_ordering(base_ring(G),T))
end

###########################################
# Counter for the steps in the Fractal Walk.
###########################################
counter = 0
function delete_counter()
    global counter
    temp = counter
    counter = 0
    return temp
end
function getcounter()
    global counter
    return counter
end
function raise_counter()
    global counter = getcounter() + 1
end
###############################################################
# Implementation of the standard walk.
###############################################################

function standard_walk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    infoLevel::Int,
)
    if infoLevel >= 1
        println("standard_walk results")
        println("Crossed Cones in: ")
    end

    Gb = standard_walk(G, S, T, S[1, :], T[1, :], infoLevel)

    return Gb
end

function standard_walk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    currweight::Vector{Int},
    tarweight::Vector{Int},
    infoLevel::Int,
)
    while true
        G = standard_step(G, currweight, T)
        if infoLevel >= 1
            println(currweight)
            if infoLevel == 2
                println(G)
            end
        end
        raise_counter()
        if currweight == tarweight
            return G
        else
            currweight = next_weight(G, currweight, tarweight)
        end
    end
end

###############################################################
# The standard step is used for the strategies standard and pertubed.
###############################################################

function standard_step(G::Oscar.IdealGens, w::Vector{Int}, T::Matrix{Int})
    R = base_ring(G)
    ordAlt = G.ord
    Gw = 0

    #check if no entry of w is bigger than Int32. If it´s bigger multiply it by 0.1 and round.
    if !checkInt32(w)
        Gw = ideal(initials(G, w))
        w, b = truncw(G, w, Gw)
        if !b
            throw(
                error(
                    "Some entries of the intermediate weight-vector $w are bigger than int32",
                ),
            )
        end
    else
        Gw = ideal(initials(G, w))
    end
    ordNew = create_order(R,w,T)
    H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)
    #H = liftGW2(G, R, Gw, H, Rn)
    H = lift(G, ordAlt, H, ordNew)
    return interreduce_walk(H)
end

###############################################################
# Generic-version of the Groebner Walk.
###############################################################

function generic_walk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    infoLevel::Int,
)
    Lm = [leading_term(g, ordering=matrix_ordering(base_ring(G), S)) for g in G]
    v = next_gamma(G, Lm, [0], S, T)
    ordNew=matrix_ordering(base_ring(G),T)
    if infoLevel >= 1
        println("generic_walk results")
        println("Crossed Cones with: ")
    end
    while !isempty(v)
        G, Lm = generic_step(G, Lm, v, ordNew)
        raise_counter()

        if infoLevel >= 1
            println(v)
            if infoLevel == 2
                println(G)
            end
        end
        G=Oscar.IdealGens(G, ordNew)
        v = next_gamma(G, Lm, v, S, T)
    end
    return G
end

function generic_step(
    G::Oscar.IdealGens,
    Lm::Vector{T},
    v::Vector{Int},
    ord::MonomialOrdering
) where {T<:MPolyRingElem}
    facet_Generators = facet_initials(G, Lm, v)
    H = groebner_basis(ideal(facet_Generators), ordering = ord, complete_reduction = true, algorithm = :buchberger)
    H, Lm = lift_generic(gens(G), Lm, gens(H), ord)
    G = interreduce(H, Lm, ord)
    return G, Lm
end

###############################################################
#Pertubed-version of the Groebner Walk.
###############################################################

function pertubed_walk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    p::Int,
    infoLevel::Int,
)
    if infoLevel >= 1
        println("pertubed_walk results")
        println("Crossed Cones in: ")
    end

    currweight = pertubed_vector(G, S, p)

    while true
        tarweight = pertubed_vector(G, T, p)
        Tn = add_weight_vector(tarweight, T)
        G = standard_walk(G, S, Tn, currweight, tarweight, infoLevel)
        if same_cone(G, T)
            return G
        else
            p = p - 1
            currweight = tarweight
            S = Tn
        end
    end
end

###############################################################
# The Fractal Walk
###############################################################

##########################################
# global weightvectors
##########################################
pTargetWeights = []
pStartWeights = []
firstStepMode = false

###############################################################
# Combined version of the extensions of the Fractal Walk.
# This version
# - checks if the starting weight vector represents the monomial order and pertubes it if necessary.
# - analyses the Groebner basis Gw of the initialforms and uses the Buchberger-algorithm if the generators of Gw are binomial or less.
# - skips a step in top level in the last step.
# - checks if an entry of an intermediate weight vector is bigger than int32. In case of that the Buchberger-Algorithm is used to compute the Groebner basis of the ideal of the initialforms.
###############################################################
function fractal_walk_combined(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    infoLevel::Int,
)
    if infoLevel >= 1
        println("fractal_walk_combined results")
        println("Crossed Cones in: ")
    end

    global pTargetWeights =
        [pertubed_vector(G, T, i) for i = 1:nvars(base_ring(G))]
    return fractal_walk_combined(G, S, T, S[1, :], pTargetWeights, 1, infoLevel)
end

function fractal_walk_combined(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    currweight::Vector{Int},
    pTargetWeights::Vector{Vector{Int}},
    p::Int,
    infoLevel::Int,
)
R = base_ring(G)
G.isGB = true
w = currweight
ordAlt = G.ord

    # Handling the weight of the start order.
    if (p == 1)
        if !ismonomial(initials(G, w))
            global pStartWeights = [pertubed_vector(G, S, i) for i = 1:nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        w = pStartWeights[p]
    else
        w = currweight
    end

    # main loop
    while true
        t = next_weightfr(G, w, pTargetWeights[p])

        # Handling the final step in the current depth.
        # Next_weightfr may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
        # -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weightfr returns 1 in the last step on every local path.
        if t == 1 && p != 1
            if same_cone(G, T)
                if infoLevel >= 1
                    println("depth $p: in cone ", currweight, ".")
                end

                # Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
                if !inCone(G, T, pTargetWeights, p)
                    global pTargetWeights =
                        [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                    if infoLevel >= 1
                        println("depth $p: not in cone ",pTargetWeights[p], ".")
                    end
                end
                return G
            end
        elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
            if inCone(G, T, pTargetWeights, p)
                if infoLevel >= 1
                    println("depth $p: in cone ", pTargetWeights[p], ".")
                end
                return G
            end
            global pTargetWeights =
                [pertubed_vector(G, T, i) for i = 1:nvars(R)]
            if infoLevel >= 1
                println("depth $p: not in cone ",pTargetWeights[p], ".")
            end
            continue
        end

        # skip a step for target monomial order lex.
        if t == 1 && p == 1
            if infoLevel >= 1
                println("depth $p: recursive call in ", pTargetWeights[p])
            end
            return fractal_walk_combined(
                G,
                S,
                T,
                w,
                pTargetWeights,
                p + 1,
                infoLevel,
            )
        else
            w = w + t * (pTargetWeights[p] - w)
            w = convert_bounding_vector(w)
            Gw = ideal(initials(G, w))

            # handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
            if !checkInt32(w)
                w, b = truncw(G, w, Gw)
                if !b
                    ordNew = create_order(R,T)
                    w = T[1, :]
                    G = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)
    
                    if !inCone(G, T, pTargetWeights, p)
                        global pTargetWeights =
                            [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                        if infoLevel >= 1
                            println("depth $p: not in cone ",pTargetWeights[p], ".")
                        end
                    end
                    return G
                end
            end
            ordNew = create_order(R,w,T)

            # converting the Groebner basis
            if (p == nvars(R) || isbinomial(gens(Gw)))
                H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

                if infoLevel >= 1
                    println("depth $p: conversion in ", w, ".")
                end
                raise_counter()
            else
                if infoLevel >= 1
                    println("depth $p: recursive call in $w.")
                end
                H = fractal_walk_combined(
                    Oscar.IdealGens(R,gens(Gw), ordAlt),
                    S,
                    T,
                    deepcopy(currweight),
                    pTargetWeights,
                    p + 1,
                    infoLevel,
                )
                global firstStepMode = false
            end
        end
        #H = liftGW2(G, R, Gw, H, Rn)
        H = lift_fractal_walk(G, H, ordNew)
        G = interreduce_walk(H)
        ordAlt = ordNew
        currweight = w
    end
end

###############################################################
# Plain version of the Fractal Walk.
# This version checks if an entry of an intermediate weight vector is bigger than int32. In case of that the Buchberger-Algorithm is used to compute the Groebner basis of the ideal of the initialforms.
###############################################################

function fractal_walk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    infoLevel::Int,
)
    if infoLevel >= 1
        println("FractalWalk_standard results")
        println("Crossed Cones in: ")
    end

    global pTargetWeights =
        [pertubed_vector(G, T, i) for i = 1:nvars(base_ring(G))]
    return fractal_recursiv(G, S, T, S[1, :], pTargetWeights, 1, infoLevel)
end

function fractal_recursiv(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    currweight::Vector{Int},
    pTargetWeights::Vector{Vector{Int}},
    p::Int,
    infoLevel::Int,
)
    R = base_ring(G)
    G.isGB = true
    w = currweight
    ordAlt = G.ord
    while true
        t = next_weightfr(G, w, pTargetWeights[p])

        # Handling the final step in the current depth.
        # Next_weightfr may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
        # -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weightfr returns 1 in the last step on every local path.        if t == 1 && p != 1
        if t == 1 && p != 1
            if same_cone(G, T)
                if infoLevel >= 1
                    println("depth $p: in cone ", currweight, ".")
                end

                # Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
                if !inCone(G, T, pTargetWeights, p)
                    global pTargetWeights =
                        [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                    if infoLevel >= 1
                        println("depth $p: not in cone ",pTargetWeights[p], ".")
                    end
                end
                return G
            end
        elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
            if inCone(G, T, pTargetWeights, p)
                if infoLevel >= 1
                    println("depth $p: in cone ",pTargetWeights[p], ".")
                end
                return G
            end
            global pTargetWeights =
                [pertubed_vector(G, T, i) for i = 1:nvars(R)]
            if infoLevel >= 1
                println("depth $p: not in cone ",pTargetWeights[p], ".")
            end
            continue
        end

        w = w + t * (pTargetWeights[p] - w)
        w = convert_bounding_vector(w)
        Gw = ideal(initials(G, w))

        # Handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
        if !checkInt32(w)
            w, b = truncw(G, w, Gw)
            if !b
                ordNew = create_order(R,T)
                w = T[1, :]
                G = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

                if !inCone(G, T, pTargetWeights, p)
                    global pTargetWeights =
                        [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                    if infoLevel >= 1
                        println("depth $p: not in cone ",pTargetWeights[p], ".")
                    end
                end
                return G
            end
        end
        ordNew = create_order(R,w,T)
        # Converting the Groebner basis
        if p == nvars(R)
            H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

            if infoLevel >= 1
                println("depth $p: conversion in ", w, ".")
            end
            raise_counter()
        else
            if infoLevel >= 1
                println("depth $p: recursive call in $w.")
            end
            H = fractal_recursiv(
                Oscar.IdealGens(R,gens(Gw), ordAlt),
                S,
                T,
                deepcopy(currweight),
                pTargetWeights,
                p + 1,
                infoLevel,
            )
        end
        #H = liftGW2(G, R, Gw, H, Rn)
        H = lift_fractal_walk(G, H, ordNew)
        G = interreduce_walk(H)
        ordAlt = ordNew
        currweight = w
    end
end

###############################################################
# Extends the plain Fractal Walk by checking the start order.
###############################################################

function fractal_walk_start_order(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    infoLevel::Int,
)
    if infoLevel >= 1
        println("fractal_walk_withStartorder results")
        println("Crossed Cones in: ")
    end

    global pTargetWeights =
        [pertubed_vector(G, T, i) for i = 1:nvars(Singular.base_ring(G))]
    return fractal_walk_recursiv_startorder(
        G,
        S,
        T,
        S[1, :],
        pTargetWeights,
        1,
        infoLevel,
    )
end

function fractal_walk_recursiv_startorder(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    currweight::Vector{Int},
    pTargetWeights::Vector{Vector{Int}},
    p::Int,
    infoLevel::Int,
)
    R = base_ring(G)
    G.isGB = true
    ordAlt = G.ord

    # Handling the starting weight.
    if (p == 1)
        if !ismonomial(initials(G, currweight))
            global pStartWeights = [pertubed_vector(G, S, i) for i = 1:nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        w = pStartWeights[p]
    else
        w = currweight
    end

    while true
        t = next_weightfr(G, w, pTargetWeights[p])

        # Handling the final step in the current depth.
        # Next_weightfr may return 0 if the target vector does not lie in the cone of T while G already defines the Groebner basis w.r.t. T.
        # -> Checking if G is already a Groebner basis w.r.t. T solves this problem and reduces computational effort since next_weightfr returns 1 in the last step on every local path.        if t == 1 && p != 1
        if t == 1 && p != 1
            if same_cone(G, T)
                if infoLevel >= 1
                    println("depth $p: in cone ", currweight, ".")
                end

                # Check if a target weight of pTargetWeights[p] and pTargetWeights[p-1] lies in the wrong cone.
                if !inCone(G, T, pTargetWeights, p)
                    global pTargetWeights =
                        [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                    if infoLevel >= 1
                        println("depth $p: not in cone ",pTargetWeights[p], ".")
                    end
                end
                return G
            end
        elseif t == [0] # The Groebner basis w.r.t. the target weight and T is already computed.
            if inCone(G, T, pTargetWeights, p)
                if infoLevel >= 1
                    println("depth $p: in cone ",pTargetWeights[p], ".")
                end
                return G
            end
            global pTargetWeights =
                [pertubed_vector(G, T, i) for i = 1:nvars(R)]
            if infoLevel >= 1
                println("depth $p: not in cone ",pTargetWeights[p], ".")
            end
            continue
        end

        w = w + t * (pTargetWeights[p] - w)
        w = convert_bounding_vector(w)
        Gw = ideal(initials(G, w))

        # Handling the current weight with regards to Int32-entries. If an entry of w is bigger than Int32 use the Buchberger-algorithm.
      if !checkInt32(w)
            w, b = truncw(G, w, Gw)
            if !b
                ordNew = create_order(R,T)
                w = T[1, :]
                G = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

                if !inCone(G, T, pTargetWeights, p)
                    global pTargetWeights =
                        [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                    if infoLevel >= 1
                        println("depth $p: not in cone ",pTargetWeights[p], ".")
                    end
                end
                return G
            end
        end
        ordNew = create_order(R,w,T)

        # Converting the Groebner basis
        if p == nvars(R)
            H = groebner_basis(Gw, ordering = ordNew, complete_reduction = true, algorithm = :buchberger)

            if infoLevel >= 1
                println("depth $p: conversion in ", w, ".")
            end
            raise_counter()
        else
            if infoLevel >= 1
                println("depth $p: recursive call in $w.")
            end
            H = fractal_walk_recursiv_startorder(
                    Oscar.IdealGens(R,gens(Gw), ordAlt),
                    S,
                    T,
                    deepcopy(currweight),
                    pTargetWeights,
                    p + 1,
                    infoLevel,
                )
            global firstStepMode = false
        end
        #H = liftGW2(G, R, Gw, H, Rn)
        H = lift_fractal_walk(G, H, ordNew)
        G = interreduce_walk(H)
        ordAlt = ordNew
        currweight = w
    end
end


###############################################################
# Tran´s version of the Groebner Walk.
# Returns the intermediate Groebner basis if an entry of an intermediate weight vector is bigger than int32.
###############################################################

function tran_walk(
    G::Oscar.IdealGens,
    S::Matrix{Int},
    T::Matrix{Int},
    infoLevel::Int,
)
    if infoLevel >= 1
        println("tran_walk results")
        println("Crossed Cones in: ")
    end

    currweight = S[1, :]
    tarweight = T[1, :]
    R = base_ring(G)
    if !ismonomial(initials(R, Singular.gens(G), currweight))
        currweight = pertubed_vector(G, S, nvars(R))
    end

    while true
        w = next_weight(G, currweight, tarweight)

        # return the Groebner basis if an entry of w is bigger than int32.
        if !checkInt32(w)
            w, b = truncw(G, w, initials(R, gens(G), w))
            if !b
                return G
            end
        end
        Rn = change_order(R, w, T)
        if w == tarweight
            if same_cone(G, T)
                if infoLevel >= 1
                    println("Cones crossed: ", counter)
                end
                return G
            elseif inSeveralCones(initials(base_ring(G), gens(G), tarweight))
                tarweight = representation_vector(G, T)
                continue
            end
        end
        G = standard_step_without_int32_check(G, w, T)
        if infoLevel >= 1
            println(w)
            if infoLevel == 2
                println(G)
            end
        end
        R = Rn
        currweight = w
        raise_counter()
    end
end

###############################################################
# Standard step without checking of the entries of a given weight vector.
###############################################################

function standard_step_without_int32_check(
    G::Oscar.IdealGens,
    w::Vector{Int},
    T::Matrix{Int},
)
    R = base_ring(G)
    Rn = change_order(R, w, T)
    Gw = initials(Rn, gens(G), w)
    H = Singular.std(Singular.Ideal(Rn, Gw), complete_reduction = true)
    #H = liftGW2(G, R, Gw, H, Rn)
    H = lift(G, R, H, Rn)
    return interreduce_walk(H)
end
