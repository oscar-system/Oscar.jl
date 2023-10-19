

mutable struct MatroidRealizationSpace
    defining_ideal::Union{Ideal, NumFieldOrdIdl, Nothing}
    inequations::Union{Vector{Oscar.RingElem},Nothing}
    ambient_ring::Union{Oscar.MPolyRing, Ring, Nothing}
    realization_matrix::Union{Oscar.MatElem,Nothing}
    realizable::Union{Bool,Nothing} # Nothing if we're not sure yet if the matroid is realizable.
    F::AbstractAlgebra.Ring
    char::Union{Int,Nothing}
    q::Union{Int,Nothing}
    one_realization::Bool
end

function Base.show(io::IO, RS::MatroidRealizationSpace)
    if RS.realizable == false
        if RS.char == nothing && RS.q == nothing && RS.F == ZZ
            print(io, "The matroid is not realizable.")
        else
            print(io, "The matroid is not realizable over the specified field or characteristic.")
        end
    else
        if RS.one_realization
            println(io, "One realization is given by")
        else
            println(io, "The realizations are parametrized by")
        end
        # println isn't ideal as it prints the matrix as one big line
        display(RS.realization_matrix)
        println(io, "in the ", RS.ambient_ring)
        I = RS.defining_ideal
        if (typeof(I) isa NumFieldOrdIdl && I.gen != ZZ(0)) || (typeof(I) isa Ideal && !iszero(I))
            println(io, "within the vanishing set of the ideal\n", RS.defining_ideal)
        end
        if length(RS.inequations) > 0
            println(io, "avoiding the zero loci of the polynomials\n", RS.inequations)
        end
    end
end

# Constructors
function MatroidRealizationSpace(F::AbstractAlgebra.Ring, char::Union{Int,Nothing}, q::Union{Int,Nothing})
    return MatroidRealizationSpace(nothing, nothing, nothing, nothing, false, F, char, q, false)
end

function MatroidRealizationSpace(I::Union{Ideal, NumFieldOrdIdl, Nothing},
    ineqs::Union{Vector{<:Oscar.RingElem},Nothing},
    R::Union{Oscar.MPolyRing, Ring, Nothing}, mat::Union{Oscar.MatElem,Nothing},
    F::AbstractAlgebra.Ring, char::Union{Int,Nothing}, q::Union{Int,Nothing})
    return MatroidRealizationSpace(I, ineqs, R, mat, nothing, F, char, q, false)
end

@doc raw"""
    is_realizable(M; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)

* If char = nothing, then this function determines whether the matroid is realizable over some field.
  
* If `char == 0`, then this function determines whether the matroid is realizable over some field of
    characteristic 0.
  
* If char = p is prime, this function determines whether the matroid is realizable
    of characteristic p.
  
* If `char == p` and `q` is a power of `p`, this function determines whether the matroid is realizable over the
    finite field ``GF(q)``.
"""
function is_realizable(M::Matroid; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)::Bool
    RS = realization_space(M, char=char, q=q)
    return is_realizable(RS)
end

function is_realizable(RS::MatroidRealizationSpace)
    !isnothing(RS.realizable) && return RS.realizable

    if !(typeof(RS.ambient_ring) isa MPolyRing)
        RS.realizable = true
        return RS.realizable
    end
    for p in minimal_primes(RS.defining_ideal)
        component_non_trivial = true
        for f in RS.inequations
            if f in p
                component_non_trivial = false
                break
            end
        end
        if component_non_trivial
            RS.realizable = true
            return RS.realizable
        end
    end
    RS.realizable = false
    return RS.realizable
end


@doc raw"""
    defining_ideal(RS::MatroidRealizationSpace)

The ideal of the matroid realization space `RS`.
"""
defining_ideal(RS::MatroidRealizationSpace) = RS.defining_ideal

@doc raw"""
inequations(RS::MatroidRealizationSpace)

Generators of the localizing semigroup of `RS`.
These are the polynomials that need to be nonzero in any realization.
"""
inequations(RS::MatroidRealizationSpace) = RS.inequations

@doc raw"""
ambient_ring(RS::MatroidRealizationSpace)

The polynomial ring containing the ideal `defining_ideal(RS)` and the polynomials in `inequations(RS)`. 
"""
ambient_ring(RS::MatroidRealizationSpace) = RS.ambient_ring

@doc raw"""
realization_matrix(RS::MatroidRealizationSpace)

A matrix with entries in ambient_ring(RS) whose columns, when filled in with values satisfying equalities
from `defining_ideal(RS)` and inequations from `inequations(RS)`, form a realization for the matroid. 
"""
realization_matrix(RS::MatroidRealizationSpace) = RS.realization_matrix

function realization_space_matrix(M::Matroid, B::Vector{Int}, F::AbstractAlgebra.Ring)
    # prepare the combinatorial data
    
    circs = fundamental_circuits_of_basis(M,B)
    
    nonIdCols = setdiff(matroid_groundset(M),B)
    circs = [setdiff( c, nonIdCols ) for c in circs]

    rk = rank(M)
    n = length(M)

    # we start by computing the number of variables:
    numVars = 0
    unUsedRowsForOnes = collect(2:rk)
    for col in 1:(n-rk)
        for row in 1:rk
            circ = circs[col]
            if !(B[row] == minimum(circ)) && B[row] in circ
                if row in unUsedRowsForOnes
                    unUsedRowsForOnes = setdiff(unUsedRowsForOnes,[row])
                else
                    numVars += 1 
                end
            end
        end
    end

    
    if numVars > 0
        R, x = polynomial_ring(F, numVars)
    else
        R = F
        x = Vector{MPolyRingElem}()
    end

    unUsedRowsForOnes = collect(2:rk)
    
    # create the matrix and fill it with entries
    
    mat = zero_matrix(R,rk,n)
    
    for i in 1:rk
        mat[i,B[i]]=R(1)
    end   
    
    varCounter = 1
    
    for col in 1:(n-rk)
        for row in 1:rk
            circ = circs[col]
            c = nonIdCols[col]
            
            if B[row] == minimum(circ)
                mat[row,c] = R(1)
            elseif B[row] in circ 
                if row in unUsedRowsForOnes
                    mat[row,c] = R(1)
                    unUsedRowsForOnes = setdiff(unUsedRowsForOnes,[row])
                else
                    mat[row,c] = x[varCounter]
                    varCounter = varCounter+1
                end
            else
                mat[row, c] = R(0)
            end
        end
    end
    
    return (R, mat)
end


# Given a matroid M with a basis B, this functions computes for a i in groundset(M)\B a circuit in B \cup i.
# This is a function used in the construction of the matrix underlying the realization space computation.
function fundamental_circuits_of_basis(M::Matroid, B::Vector{Int})
    remaining_elts = setdiff(matroid_groundset(M),B)
    fund_circs = Vector{Vector{Int}}()
    for i in remaining_elts
         push!(fund_circs,fundamental_circuit(M,B,i))
     end
    return fund_circs
 end


@doc raw"""
    realization_space(M::Matroid; B::Union{GroundsetType,Nothing} = nothing,
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false,
    char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)::MatroidRealizationSpace
        
    This function returns the data for the coordinate ring of the matroid realization space of the matroid M
    as a `MatroidRealizationSpace`. This function has several optional parameters.
    
    * `B` is a basis of M that specifies which columns of `realization_matrix(M)` form the identity matrix.
      The default is `nothing`, in which case the basis is chosen for you.
    
    * `F` is a coefficient ring for the realization space. The default is `ZZ`.
      Other options are `QQ` or `GF(p)` for some prime p.
      
    * `char` specifies the characteristic of the coefficient ring, and is used to determine if the matroid
      is realizable over a field of this characteristic. The default is `nothing`.
    
    * `q` is an integer, and when char = p, this input is used to determine whether the matroid
      is realizable over the finite field ``GF(p^{q})``. The default is `nothing`.

    * `reduce` determines whether a reduced realization space is returned which means that the equations
      are used to eliminate variables as far as possible. The default is `true`.

    * `saturate` determines whether `defining_ideal(M)` should be saturated with respect to the semigroup
      generated by `inequations(M)`. The default is `false`. This can be rather slow for large instances.

# Examples
```jldoctest
julia> M = fano_matroid();

julia> RS = realization_space(M)
The realizations are parametrized by
[0   1   1   1   1   0   0]
[1   0   1   1   0   1   0]
[1   0   1   0   1   0   1]
in the Integer Ring
within the vanishing set of the ideal
2ZZ

julia> realization_space(non_fano_matroid())
The realizations are parametrized by
[1   1   0   0   1   1   0]
[0   1   1   1   1   0   0]
[0   1   1   0   0   1   1]
in the Integer Ring
avoiding the zero loci of the polynomials
RingElem[2]

julia> realization_space(pappus_matroid(),char=0)
The realizations are parametrized by
[1   0   1   0   x2   x2                 x2^2    1    0]
[0   1   1   0    1    1   -x1*x2 + x1 + x2^2    1    1]
[0   0   0   1   x2   x1                x1*x2   x1   x2]
in the Multivariate Polynomial Ring in x1, x2 over Rational Field
avoiding the zero loci of the polynomials
RingElem[x1 - x2, x2, x1, x2 - 1, x1 + x2^2 - x2, x1 - 1, x1*x2 - x1 - x2^2]

julia> realization_space(uniform_matroid(3,6))
The realizations are parametrized by
[1   0   0   1    1    1]
[0   1   0   1   x1   x3]
[0   0   1   1   x2   x4]
in the Multivariate Polynomial Ring in x1, x2, x3, x4 over Integer Ring
avoiding the zero loci of the polynomials
RingElem[x1*x4 - x2*x3, x2 - x4, x1 - x3, x1*x4 - x1 - x2*x3 + x2 + x3 - x4, x3 - x4, x4 - 1, x3 - 1, x3, x4, x1 - x2, x2 - 1, x1 - 1, x1, x2]

```
"""
function realization_space(M::Matroid; B::Union{GroundsetType,Nothing} = nothing, 
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false, simplify::Bool=true,
    char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)::MatroidRealizationSpace

    if char!=nothing && !isprime(char) && char!=0
        error("The characteristic has to be 0 or a prime number.")
    end

    #Construct the base ring as F_p if q=p^k
    if q!=nothing
        isprimepower, p, k = is_prime_power_with_data(q)
        if !isprimepower
            error("The given q has to be a prime power.")
        end
        if char!=nothing && char!=p
            error("The given characteristic doesn't match q.")
        else
            char = p
        end
    end

    if char == 0
        F = QQ
    elseif char != nothing
        F = GF(char)
    end

    rk = rank(M)
    n = length(M)

    goodM = isomorphic_matroid(M, [i for i in 1:n])

    Bs = bases(goodM)
    
    if !isnothing(B) 
        goodB = sort!(Int.([M.gs2num[j] for j in B]))
    else
        goodB = find_good_basis_heuristically(goodM)
    end
    polyR, mat = realization_space_matrix(goodM, goodB, F)

    
    eqs = Vector{RingElem}()
    ineqs = Vector{RingElem}()

    #need to catch the corner-case if there are no variables at all
    if !(typeof(polyR) isa MPolyRing)
        RS = MatroidRealizationSpace(ideal(polyR,0), ineqs, polyR, mat, F, char, q)
        RS.realizable = true
        return RS
    end
    
    for col in subsets(Vector(1:n),rk)
        
        col_det = det(mat[:,col])
        
        if total_degree( col_det ) <= 0 
            
            if col_det != 0 && col in Bs 
                if isunit( col_det ) 
                    continue
                end
            elseif col_det != 0 # and col is not a basis
                # determinant nonzero but set not a basis
                push!( eqs, col_det )
            elseif col in Bs 
                #determinant zero but set is a basis, i.e. M is not realizable
                return MatroidRealizationSpace(F, char, q)
            else
                continue
            end
            
        end
        
        if  col in Bs
            push!( ineqs, col_det )
        else
            push!( eqs, col_det )
        end
        
    end

    if q != nothing
        for i in 1:polyR.nvars
            push!(eqs, polyR[i]^q-polyR[i])
        end
    end

    def_ideal = ideal(polyR,eqs)
    def_ideal = ideal(groebner_basis(def_ideal))
    isone(def_ideal) && return MatroidRealizationSpace(F, char, q)

    # Unclear if we should use this. Can be too slow in some cases.
    #=
    if !iszero(def_ideal)
        for i in 1:length(ineqs)
            ineqs[i] = Oscar.reduce(ineqs[i], gens(def_ideal))
        end
    end
    (polyR(0) in ineqs) && return MatroidRealizationSpace(nothing, nothing, nothing, nothing, false, F, char, q)
    =#

    ineqs = gens_2_prime_divisors(ineqs)
    #=
    # Only saturate easy inequations as it's too slow otherwise.
    for i in 1:length(ineqs)
        if count(k -> k!=0, degrees(ineqs[i])) < 3 && sum(degrees(ineqs[i])) < 3
            def_ideal = saturation(def_ideal, ideal([ineqs[i]]))
        end
    end
    =#
    isone(def_ideal) && return MatroidRealizationSpace(F, char, q)
    
    RS = MatroidRealizationSpace(def_ideal, ineqs, polyR, mat, F, char, q)

    if simplify
        RS = reduce_realization_space(RS)
    end

    if saturate
        RS.defining_ideal = stepwise_saturation(def_ideal,ineqs)
        isone(RS.defining_ideal) && return MatroidRealizationSpace(F, char, q)
        RS.realizable = true
    end

    return RS
end


# A heuristic function that tries to find a sensible basis for the moduli space computation for which the defining ideal is not too complicated
function find_good_basis_heuristically(M::Matroid)
    bs = bases(M)
    cs = circuits(M)
    min_num_vars = length(cs)*rank(M)
    min_basis = bs[1]
    for bi in 1:length(bs)
        current_num_vars = 0
        for c in cs
            for e in c
                @inbounds if e in bs[bi]
                    current_num_vars += 1
                end
            end
        end
        if current_num_vars < min_num_vars
            min_num_vars = current_num_vars
            min_basis = bs[bi]
        end
    end
    return min_basis
end


# Return the prime divisors of f. 
function poly_2_prime_divisors(f::RingElem)
    return map(first, factor(f))
end

# Return the unique prime divisors of the elements of Sgen, again no exponents. 
function gens_2_prime_divisors(Sgens::Vector{<:RingElem})
    return unique!(vcat([poly_2_prime_divisors(f) for f in Sgens]...))
end


function stepwise_saturation(I::MPolyIdeal, Sgens::Vector{<:RingElem})
    for f in Sgens
        I = saturation(I, ideal([f]))
    end
    return I
end

@doc raw"""
    realization(M::Matroid; B::Union{GroundsetType,Nothing} = nothing, 
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false, 
    char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)::MatroidRealizationSpace
        
    This function tries to find one realization in the matroid realization space of the matroid `M`.
    The output is again a `MatroidRealizationSpace`. 

    If the matroid is only realizable over an extension of the prime field the extension field is specified
    as a splitting field of an irreducible polynomial. Every root of this polynomial gives an equivalent
    realization of the matroid. 

    This function has several optional parameters. Note that one must input either the characteristic or a specific
    field of definition for the realization. 
    
    * `B` is a basis of M that specifies which columns of `realization_matrix(M)` form the identity matrix.
      The default is `nothing`, in which case the basis is chosen for you. 
    
    * `F` is a coefficient ring for the realization space. The default is `ZZ`.
      Other options are `QQ` or `GF(p)` for some prime p.
      
    * `char` specifies the characteristic of the coefficient ring, and is used to determine if the matroid
      is realizable over a field of this characteristic. The default is `nothing`.
    
    * `q` is an integer, and when char = p, this input is used to determine whether the matroid 
      is realizable over the finite field ``GF(p^{q})``. The default is `nothing`.

    * `reduce` determines whether a reduced realization space is returned which means that the equations
      are used to eliminate variables as far as possible. The default is `true`.

    * `saturate` determines whether `defining_ideal(M)` should be saturated with respect to the semigroup
      generated by `inequations(M)`. The default is `false`. This can be rather slow for large instances.

# Examples
```jldoctest
julia> realization(pappus_matroid(), char=0)
One realization is given by
[1   0   1   0   2   2   4   1   0]
[0   1   1   0   1   1   1   1   1]
[0   0   0   1   2   3   6   3   2]
in the Rational Field

julia> realization(pappus_matroid(), q=4)
One realization is given by
[1   0   1   0        1        1    1        1    0]
[0   1   1   0   x1 + 1   x1 + 1   x1        1    1]
[0   0   0   1        1       x1   x1   x1 + 1   x1]
in the Multivariate Polynomial Ring in x1 over Galois field with characteristic 2
within the vanishing set of the ideal
ideal(x1^2 + x1 + 1)

julia> realization(uniform_matroid(3,6), char=5)
One realization is given by
[1   0   0   1   1   1]
[0   1   0   1   4   3]
[0   0   1   1   3   2]
in the Galois field with characteristic 5

```
"""
function realization(M::Matroid; B::Union{GroundsetType,Nothing} = nothing, 
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false, simplify::Bool=true,
    char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)

    RS = realization_space(M,B=B,F=F,saturate=saturate,simplify=simplify,char=char,q=q)

    return realization(RS)
end


@doc raw"""
realization(RS::MatroidRealizationSpace)

This function tries to find one realization in the matroid realization `RS`.
    The output is again a `MatroidRealizationSpace`.
"""
function realization(RS::MatroidRealizationSpace)
    # If the matroid is not realizable we stop
    if RS.char == nothing && RS.q == nothing && RS.F == ZZ
        error("A field or characteristic must be specified")
    end

    if RS.realizable == false
        return RS
    end

    # If the ambient ring is not a polynomial ring we can reduce we stop
    R = RS.ambient_ring
    !(typeof(R) isa MPolyRing) && return RS    
    Inew = RS.defining_ideal
    eqs = copy(gens(Inew))

    d = min(dim(Inew), nvars(R))
    if d == 0
        nvars(R) == 0 && return RS
        for p in gens(Inew)
            f = factor(p)
            if length(f) > 1
                push!(eqs,collect(f)[1][1])
            end
        end
    end
    ineqsnew = RS.inequations

    counter = 0
    base = 7
    if RS.char != nothing && RS.char > 0
        base = RS.char
    end
    upperbound = min(base^d, 10^3)
    while counter < upperbound
        eqsnew = copy(eqs)
        values = digits(counter, base=base,pad=d)
        counter += 1

        for i in 1:d
            push!(eqsnew, R[i]-values[i])
        end
        Inew = ideal(groebner_basis(ideal(R,eqsnew)))
        isone(Inew) && continue
        ineqsnew = copy(RS.inequations)
        for i in 1:length(ineqsnew)
            ineqsnew[i] = reduce(ineqsnew[i], gens(Inew))
        end
        (R(0) in ineqsnew) && continue
        Inew = stepwise_saturation(Inew,ineqsnew)
        isone(Inew) && continue
        break
    end
    counter == upperbound && d != 0 && return RS

    RSnew = MatroidRealizationSpace(Inew, ineqsnew, R, RS.realization_matrix, RS.F, RS.char, RS.q)
    RSnew = reduce_realization_space(RSnew)
    ineqsnew = RSnew.inequations
    if length(ineqsnew) > 0
        Inew = RSnew.defining_ideal
        Rnew = RSnew.ambient_ring
        ineqsnew = filter(p ->!isone(ideal(groebner_basis(Inew+ideal(Rnew,p)))),ineqsnew)
        RSnew = MatroidRealizationSpace(Inew, ineqsnew, Rnew, RSnew.realization_matrix, RSnew.F, RSnew.char, RSnew.q)
    end

    RSnew.one_realization = true

    return RSnew
end


#####################
# full reduction    #
#####################

# computes the coefficient of v in monomial m. 
function coefficient_monomial(v::RingElem, m::RingElem)
    isone(degree(m,v)) || return "The variable is not a degree 1 factor."
    mf = factor(m)
    u = unit(mf)
    not_v = [k^e for (k,e) in mf if k != v ]
    length(not_v) == 0 ? u : u*prod(not_v)
end

# computes the coefficient of v in f. 
function coefficient_v(v::RingElem, f::RingElem)
    isone(degree(f,v)) || return "degree of variable must be 1"
    withv = [term(f,i) for i in 1:length(f) if v in vars(monomial(f,i))]
    return sum([coefficient_monomial(v,m) for m in withv])
end


function find_solution_v(v::RingElem, Igens::Vector{<:RingElem}, 
                         Sgens::Vector{<:RingElem}, R::MPolyRing) 

    
    with_v_deg_1 = [g for g in Igens if isone(degree(g,v))] 
    length(with_v_deg_1) != 0 || return "can't isolate"

    for f in with_v_deg_1

        den = coefficient_v(v, f)
        fac_den = poly_2_prime_divisors(den)
        !issubset(fac_den, Sgens) && continue

        no_v = [term(f,i) for i in 1:length(f) if !(v in vars(monomial(f,i)))]
        iszero(length(no_v)) && continue
              
        h = R(-1)*sum(no_v)
        return h//den
        
    end
    return "can't solve for v"
end


# v is replaced by t in f
function sub_map(v::RingElem, t::RingElem, R::MPolyRing, xs::Vector{<:RingElem}) 
    xs_v = map(x -> x==v ? t : x, xs )    
    return hom(R,FractionField(R), a->a, xs_v)
end



# replace v by t in f, only return the numerator.
function sub_v(v::RingElem, t::RingElem, f::RingElem, R::Ring, xs::Vector{<:RingElem}) 
    m = sub_map(v,t,R,xs) 
    new_f = numerator(m(f))
    return new_f 
end


# removes factors that are in the semigroup generated by Sgens
function clean(f::RingElem, R::MPolyRing, Sgens::Vector{<:RingElem})   
    
    fFactors = factor(f)
    FactorsDict = Dict(fFactors) 
    cleanf_arr = [k^(FactorsDict[k]) for k in keys(FactorsDict) if !(k in Sgens) || is_unit(k)]
    
    length(cleanf_arr) > 0 ? prod(cleanf_arr) : unit(fFactors)
    
end

#variables in ideal
function ideal_vars(Igens::Vector{<:RingElem}) 
    return unique!(vcat([vars(gen) for gen in Igens]...))
end

function n_new_Sgens(x::RingElem, t::RingElem, Sgens::Vector{<:RingElem}, 
                     R::Ring, xs::Vector{<:RingElem}) 
    preSgens = unique!([sub_v(x, t, f, R, xs) for f in Sgens])
    return gens_2_prime_divisors(preSgens)
end

function n_new_Igens(x::RingElem, t::RingElem, Igens::Vector{<:RingElem}, 
                     Sgens::Vector{<:RingElem}, R::Ring, xs::Vector{<:RingElem}) 

    preIgens = unique!([clean(sub_v(x, t, f, R, xs), R, Sgens) for f in Igens])
    return filter(!iszero, preIgens)
end



function matrix_clear_den_in_col(X::Oscar.MatElem, c::Int)
    Xc = [denominator(f) for f in X[:, c]]
    t = lcm(Xc)
    result = multiply_column!(X, t, c)
    return result
end



function matrix_clear_den(X::Oscar.MatElem)
    rs, cs = size(X)
    for c in 1:cs
        X = matrix_clear_den_in_col(X, c)
    end
    return X
end


function reduce_ideal_one_step(MRS::MatroidRealizationSpace, 
                               elim::Vector{<:RingElem}, 
                               fullyReduced::Bool)

    Igens = gens(MRS.defining_ideal)
    Sgens = MRS.inequations
    R = MRS.ambient_ring
    FR = FractionField(R)
    xs = gens(R)
    X = MRS.realization_matrix
    nr, nc = size(X)
    
    Ivars = ideal_vars(Igens);

    for x in Ivars 
        t = find_solution_v(x, Igens, Sgens, R)
        t isa String && continue
        
        phi = sub_map(x, t, R, xs)
        
	    Sgens_new = n_new_Sgens(x, t, Sgens, R, xs);
        if length(Sgens_new) == 0
            Sgens_new = Vector{RingElem}()
        end
        Igens_new = n_new_Igens(x, t, Igens, Sgens_new, R, xs);
        push!(elim, x)
        
        phiX = matrix(FR, [phi(X[i,j]) for i in 1:nr, j in 1:nc  ] )
        nX_FR = matrix_clear_den(phiX)
        nX = matrix(R, [numerator(nX_FR[i,j])  for i in 1:nr, j in 1:nc ])
        
        GBnew = collect(groebner_basis(ideal(R, Igens_new)))         
        
        MRS_new = MatroidRealizationSpace(ideal(R, GBnew), Sgens_new, R, nX, MRS.F, MRS.char, MRS.q )

        
        return (MRS_new, elim, fullyReduced)
    end

    return (MRS, elim, true)

end


function reduce_realization_space(MRS::MatroidRealizationSpace,
                               elim::Vector{RingElem} = Vector{RingElem}(), 
                               fullyReduced::Bool = false) 
    
    
    #If there are no variables left, we don't reduce anything
    if !(typeof(MRS.ambient_ring) isa MPolyRing)
        return MRS
    end

    output = reduce_ideal_one_step(MRS, elim, fullyReduced)
    output isa String && return "Not Realizable 0 in Semigroup"
    (MRS, elim, fullyReduced) = output    
    
    !fullyReduced && return reduce_realization_space(MRS, elim, fullyReduced)
    
    
    R = MRS.ambient_ring
    xs = gens(R)
    cR = coefficient_ring(R)
    X = MRS.realization_matrix
    nr, nc = size(X)
    Igens = gens(MRS.defining_ideal)
    Sgens = MRS.inequations

    
    zero_elim = []        
    for i in 1:length(xs)
        if xs[i] in elim
            push!(zero_elim, 0)
        else
            push!(zero_elim, "x"*string(i) ) 
        end
    end
        
    xnew_str = Vector{String}(filter(x -> x!=0,  zero_elim))    
    
    
    
    if length(xnew_str) == 0
        phi = hom(R, cR, a->a, [cR(0) for i in 1:length(xs)])
        ambR = codomain(phi);
        if length(Igens) == 0
            Inew = ideal(ambR, ambR(0))
        else
            Inew = ideal(ambR, phi.(Igens)); 
        end
        
        if length(Sgens) == 0
            normal_Sgens = Vector{RingElem}()
        else             
            normal_Sgens = phi.(Sgens)
        end        
        
        
    else
        Rnew, xnew = polynomial_ring(coefficient_ring(R), length(xnew_str)) 
    
        zero_elim_var = []
        j=1
        for i in 1:length(zero_elim)
            if xs[i] in elim
                push!(zero_elim_var, Rnew(0))
            else
                push!(zero_elim_var, xnew[j] ) 
                j+=1
            end
        end
        
        phi = hom(R, Rnew, a->a, zero_elim_var)
        ambR = codomain(phi);
        if length(Igens) == 0
            Inew = ideal(ambR, ambR(0))
        else
            Inew = ideal(ambR, phi.(Igens)); 
        end
        
        if length(Sgens) == 0
            normal_Sgens = Vector{RingElem}()
        else             
            Sgens_new = phi.(Sgens)
            normal_Sgens = gens_2_prime_divisors([normal_form(g, Inew) for g in Sgens_new])
        end
    
    end

    if isone(Inew)
        return MatroidRealizationSpace(MRS.F, MRS.char, MRS.q)
    end
    
    Xnew = matrix(ambR, [phi(X[i,j]) for i in 1:nr, j in 1:nc])

    #Try to reduce the matrix one last time using the ideal and the inequations
    m, n = size(Xnew)
    A = base_ring(R)
    if length(gens(Inew)) > 0 && gens(Inew)[1] != 0
        for i in 1:m
            for j in 1:n
                if A == ZZ
                    Xnew[i,j] = mod(Xnew[i,j], gens(Inew)[1])
                else
                    Xnew[i,j] = reduce(Xnew[i,j], gens(Inew))
                end
            end
        end
    end

    for j in 1:n
        g = gcd(Xnew[:,j]...)
        prime_divisors = poly_2_prime_divisors(g)
        for f in prime_divisors
            if f in normal_Sgens
                for i in 1:m
                    Xnew[i,j] = Xnew[i,j]/f
                end
            end
        end
    end
    MRS_new = MatroidRealizationSpace(Inew, normal_Sgens, ambR, Xnew, MRS.F, MRS.char, MRS.q)

    return MRS_new
end


