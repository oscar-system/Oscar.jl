module PrimDec
using Oscar
import Markdown
import Singular

include("ZeroDec.jl")
#Note: primary_decomposition in now in src/Rings/mpoly.jl and uses Singular.LibPrimdec

######################################### PRIMARY DECOMPOSITION VIA GTZ #################################################################################
# This is an implementation of the GTZ Algorithm, primarily based on Chapter 4 of "A SINGULAR introduction to Commutative Algebra" by Greuel and Pfister.
#
# Main algorithms:
# decomp
# reductionToZero
#########################################################################################################################################################

######################################### decomp ####################################################################
#initializes decomp_internal, which is called recursively to compute the primary decomposition with the method of Algorithm 4.3.4
#This is an implementation of Algorithm 4.3.4.
#REQUIRES:  INPUT sideal I
#           I has a field as coefficient ring
#           basering is not a quotient ring
#           basering has a global ordering
#           basering has characteristic 0
#OUTPUT:    Primary decomposition of I with associated primes in a list 
@doc Markdown.doc"""
    decomp(I::Oscar.ideal; usefglm::Bool)

Computes the primary decomposition of an Ideal $I$ over a basefield of characteristig 0 that has a global ordering,
via GTZ. Returns the primary decomposition with associated primes in a list. If the additional boolean 'usefglm' is
set to 'true' then the FGLM - algorithm is used in zero-dimensional computations.
"""
function primary_decomposition(I::Oscar.MPolyIdeal; usefglm::Bool = false)
    Oscar.singular_assure(I)
    B = base_ring(I)
    trivial = Vector{Singular.sideal}(undef, 0)
    n = 1
    Is = I.gens.S
    J = changeOrderOfBasering(I.gens.S, :lex)
    check = oneideal(J.base_ring)
    PDint = decomp_internal(J, trivial, n, check, usefglm)
    
    #Map back to original Oscar-ring:
    R = Is.base_ring
    S = J.base_ring
    phi = Singular.AlgebraHomomorphism(S, R, gens(R))
    PDres = Vector{Oscar.MPolyIdeal}(undef, 0)
    
    for i in 1:length(PDint)
      push!(PDres, Oscar.MPolyIdeal(B, phi(PDint[i])))
    end

    return PDres
end

function decomp_internal(I::Singular.sideal, primary::Array{Singular.sideal}, count::Int64, check::Singular.sideal, usefglm::Bool = false)
###########################################check whether we are already done computing############################################################

    R = I.base_ring
    I = Singular.std(I, complete_reduction = true)
    copy_I = Singular.deepcopy(I)
    copy_I.isGB = true


    #check if check is contained in I
    iszero(Singular.reduce(check, copy_I)) && return primary
        

    if Singular.reduce(R(1), copy_I) == R(0)
        
        return primary
    end
    
    X = Singular.gens(R)
    
############################################ check whether we are already zerodimensional #########################################################    
    u = Singular.maximal_independent_set(copy_I)
    if u == []
        # if maximal independent set is empty, it is zerodimensional, hence we have to return the zerodimensional primary decomposition
        n_zerodecomp = zerodecomp(copy_I, true, usefglm)
        for i in 1:length(n_zerodecomp)
           push!(primary, n_zerodecomp[i])
        end
        return primary
    end


############################################ actual computation begins here ########################################################################

#prepare and then do dimension reduction (see reductionToZero)
    ureturn = Singular.deepcopy(u)
   
    x_minus_u = filter((x) -> !(x in u), X)
    #switch lexicographic ordering to K[x\u, u]
    switch_morphism, switch_morphism_inverse, T = switchLexOrdering(copy_I, x_minus_u, u)
    copy_I = switch_morphism(copy_I)

    u, G, hdash = reductionToZero(copy_I, ureturn, x_minus_u)
#change ring to K(u)[x\u]
    n = Singular.ngens(G)
    X = Singular.gens(R)
    x_minus_u = filter((x) -> !(x in u), X)

    T = G.base_ring
    K, u = Singular.FunctionField(Singular.QQ, string.(u))
    S, x_minus_u = Singular.PolynomialRing(K, string.(x_minus_u), ordering=:lex)

    ptr1 = Singular.libSingular.n_SetMap(base_ring(T).ptr, base_ring(S).ptr)
    Q = Singular.Ideal(S, S(0))
    transitionA = vcat(0, -1*collect(1:length(u)), collect(1:length(x_minus_u)))
    gi = gens(G)
    
    for i in 1:n
        ptr2 = Singular.libSingular.p_PermPoly(gi[i].ptr, Int32[transitionA...], T.ptr, S.ptr, ptr1, Int32[])

        p = S(ptr2)

        Q = addGenerator(Q, p)
    end
    Q = Singular.std(Q; complete_reduction = true)

    #compute primary decomposition of Q ( which is G over K(u)[x\u]). This is zero-dimensional by the dimension reduction step.
    qprimary = zerodecomp(Q, true, usefglm)

    #change ring back to K[x] and compute the intersection of all ideals in the primary decomposition and all associated primes with K[x]
    primary, check = map_back_and_intersect(S, T, length(u), length(x_minus_u), qprimary, primary, switch_morphism_inverse, check)

    #transform hdash so that it has the right parent
    #substitute by morphism

    hdash = switch_morphism_inverse(hdash)
    
    copy_I = switch_morphism_inverse(copy_I)
   
    #Add h (obtained from the reduction step) to the ideal and 
    I_split = addGenerator(I, hdash)

    #if h was in I before, we are done
    if Singular.reduce(hdash, Singular.std(copy_I)) == parent(hdash)(0)
        return primary
    end
    
    
    return decomp_internal(I_split, primary, count, check, usefglm)
    
end


######################################### reductionToZero ####################################################################
#
#This is an implementation of Algorithm 4.3.2.
#REQUIRES:  INPUT sideal I in K[x], u maximal independent subset of x wrt. I, x\u
#           I has a field as coefficient ring
#           basering is not a quotient ring
#           basering has lexicographic ordering w.r.t K[x\u. u] (x\u > u)
#           basering has characteristic 0
#OUTPUT:    u mapped to the correct ring
#           G a Groebner basis of I w.r.t lex:= [x\u , u] over K[]
#           h in K(u) s.t. IK(u)[x\u] intersecting K[x] = I:<h> = I:<h^infinity> 

function reductionToZero(I::Singular.sideal, u::Vector{Singular.spoly{U}}, x_minus_u::Vector{Singular.spoly{U}}) where U <: Singular.FieldElem
    ureturn = Singular.deepcopy(u)
    G = Singular.std(I)
    T = G.base_ring

    #map to K(u)[x\u]
    K, u = Singular.FunctionField(Singular.QQ, string.(u))
    S, x_minus_u = Singular.PolynomialRing(K, string.(x_minus_u))
    ptr1 = Singular.libSingular.n_SetMap(base_ring(T).ptr, base_ring(S).ptr)
    h = K(1)
    transitionA = vcat(0, -1*collect(1:length(u)), collect(1:length(x_minus_u)))
    #compute the product of the leadiing coefs of the gröbner basis of I, considered as polynomials in x\u with coeffs in K(u) to compute h
    for i in 1:Singular.ngens(G)
        ptr2 = Singular.libSingular.p_PermPoly(G[i].ptr, Int32[transitionA...], T.ptr, S.ptr, ptr1, Int32[]) #(poly.ptr, Int32[0, variablen auf parameter], Ring 1, Ring2. ptr1, Int32[parameter auf variable]) 
        p = S(ptr2)
        h = h*leading_coefficient(p)
    end
    #map h back to K[x\u, u]
    h = Singular.n_transExt_to_spoly(h, parent_ring = T)
    
    H = Singular.Ideal(T, h)

    m =saturationExponent(G,H)
    #return the output
    if m == 0
        return ureturn, G, h
    end
    h = h^m
    return ureturn, G, h
end




####################################### helpprocs, mainly for mapping forwards and backwards between different orderings #######################

#this is used to compute the mapping back to K[x] and intersecting the primary ideals and primes 
function map_back_and_intersect(S::Singular.PolyRing, T::Singular.PolyRing, lengthu::Int64, lengthx_minus_u::Int64, qprimary::Vector{W}, primary::Vector{Singular.sideal }, switch_morphism_inverse, check::Singular.sideal) where { U <: Singular.FieldElem, V <: Singular.spoly, W <: Singular.sideal}
    S = qprimary[1].base_ring
    for j in 1:length(qprimary)
        qprimary_gens = Singular.gens(qprimary[j])
        newgens = Vector{Singular.spoly}(undef,length(qprimary_gens))
        K = base_ring(S)
        #map all elements of the primary decomposition and their associated primes back to original K[u,x\u] ring
        ptr3 = Singular.libSingular.n_SetMap(base_ring(S).ptr, base_ring(T).ptr)
        h = K(1)
        
        for i in 1:length(qprimary_gens)
            h = h*leading_coefficient(qprimary_gens[i])

            qg = qprimary_gens[i]
            ptr4 = Singular.libSingular.p_PermPoly(qg.ptr, Int32[0,collect(lengthu+1:lengthu+lengthx_minus_u)...], S.ptr, T.ptr, ptr3, Int32[collect(1:lengthu)...])

            newgens[i] = T(ptr4)

        end
        h = Singular.n_transExt_to_spoly(h, parent_ring = T)
        H = Singular.Ideal(T, h)
        new_qprimary = Singular.Ideal(T, newgens...)
        #saturate with H to compute intersection with K[X]
        new_qprimary_intersect = Singular.saturation(new_qprimary, H)
        #switch back to standard lex ordering

        new_qprimary_intersect = switch_morphism_inverse(new_qprimary_intersect)
       
        check = Singular.intersection(check, new_qprimary_intersect)
        push!(primary, new_qprimary_intersect)
    end    

    return primary, check
end

#using the method of Greuel/Pfister 1.8.9 (see Example 1.8.15), computes saturationExponent
function  saturationExponent(G::Singular.sideal, h::Singular.sideal)
    J = Singular.quotient(G, h);
    I1 = Singular.deepcopy(G);
    k = 0;

    while Singular.ngens(Singular.std(Singular.reduce(J, Singular.std(I1))))!=0 
        k+=1;
        I1 = J;
        J = Singular.quotient(I1,h);
    end
    return k
end

#Switches basering of ideal I to a custom lexicographic ordering where the variables are given by the Array H,returns both the forwards and backwards map as an AlgebraHomomorphism, and the ring that is mapped into.
function customLexOrdering(I:: Singular.sideal, H::Vector{T}) where T <: Singular.spoly
    R = I.base_ring;
    Hstring = ["$i" for i in H]
    S, = Singular.PolynomialRing(R.base_ring, Hstring, ordering=:lex)

    index_H = Array{Singular.spoly}(undef, length(H))
    map_back = Array{Singular.spoly}(undef, length(H))
    for i in 1:length(index_H)

        index_H[Singular.var_index(H[i])] = Singular.gen(S, i)
        map_back[i] = H[i]
    end



    res1 = Singular.AlgebraHomomorphism(R, S, index_H);
    res2 = Singular.AlgebraHomomorphism(S, R, map_back)
    

    return res1, res2, S
end

#Calls customLexOrdering for our case.
function switchLexOrdering(I:: Singular.sideal, x_minus_u::Vector{T}, u::Vector{T}) where T <: Singular.spoly 
    H = vcat(u, x_minus_u);
    return customLexOrdering(I,H)
end

end

using .PrimDec
