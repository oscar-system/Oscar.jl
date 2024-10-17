## Functions visible on the outside
export get_monoid_algebra
export irreducible_res
export injective_res

# export mod_quotient
# export get_all_ass_primes
# export quotient_by_ideal
# export shifted_module
# export compute_shift
# export face_to_prime
# export get_faces_of_polyhedral_cone
# export get_bounding_hyperplanes
# export get_prime_of_face
# export get_polyhedral_cone
# export get_zonotope
# export mod_quotient_2
# export ZF_basis
# export coefficients_incl
# export mod_saturate
# export irreducible_hull
# export generators_W_H
# export rational_to_integer_vector
# export get_hyperplane_H_presentation

# # Data 
# export FaceQ 
# export IrrSum
# export IndecInj
# export IrrSum
# export MonoidAlgebra
# export InjMod
# export HyperplaneQ
# export MonoidAlgebraModule

#########################
# some composite types
#########################

struct FaceQ # face of semigroup  
    prime::Ideal 
    poly::Polyhedron  #face of Q corresponding to prime
end

struct HyperplaneQ # a hyperplane bounding the cone RR_{\geq 0}Q
    hyperplane::Polyhedron
    A::Matrix{Int}
    b::Vector{Int}
end

struct IrrHull # irreducible hull
    faces::Vector{FaceQ}
    vectors::Vector{Tuple{SubquoModuleElem,Vector{Int}}}
    lambda::AbstractAlgebra.Generic.MatSpaceElem
end

struct IndecInj # indecomposable injective
    prime::Ideal
    vector::Vector{Int}
end

struct IrrSum # irreducible sum
    mod::SubquoModule
    components::Vector{IndecInj}
end
struct IrrRes # irreducible resolution (including all computed data and the cochain complex)
    irrSums::Vector{IrrSum}
    cochainMaps::Vector{SubQuoHom}
    surjections::Vector{SubQuoHom}
    inclusions::Vector{SubQuoHom}
    cokernels::Vector{SubquoModule}
    cochainComplex::ComplexOfMorphisms # if sequence not exact return trivial cochain_complex (M0 -> M0)
end

struct MonoidAlgebra # monoid algebra with associated data
    algebra::Union{MPolyRing, MPolyQuoRing}
    cone::Polyhedron
    faces::Vector{FaceQ}
    hyperplanes::Vector{HyperplaneQ}
    # normal::Bool
    pointed::Bool #cone pointed
    zonotope::Tuple{Polyhedron,Vector{Int}}
end

struct MonoidAlgebraIdeal
    monoid_algebra::MonoidAlgebra
    ideal::Ideal
end

struct MonoidAlgebraModule
    baseAlgebra::MonoidAlgebra
    mod::SubquoModule
end

struct InjMod
    monoidAlgebra::MonoidAlgebra
    indecInjectives::Vector{IndecInj}
end

struct InjRes
    injMods::Vector{InjMod}
    cochainMaps::Vector{MatElem}
    upto::Int
    irrRes::IrrRes
    shift::Vector{Int}
end

# given a graded polynomial ring or a quotient, return the MonoidAlgebra if it is one
# INPUT:    MpolyRing or MPolyQuoRing k_Q
# OUTPUT:   MonoidAlgebra(...)
function get_monoid_algebra(k_Q::Union{MPolyRing, MPolyQuoRing})
    #check if monoid algebra
    @req is_zm_graded(k_Q) "Not ZZ^d-graded."
    gg_Q = grading_group(k_Q)
    @req is_free(gg_Q) && is_abelian(gg_Q) "Not a monoid algebra."

    # polyhedral cone RR_{\geq 0}Q
    C_Q = get_polyhedral_cone(k_Q) 
    P_Q = polyhedron(C_Q) # C_Q as polyhedron
    G_Q,c = get_zonotope(P_Q)
    return MonoidAlgebra(k_Q,P_Q,get_faces_of_polyhedral_cone(k_Q,G_Q,P_Q),get_bounding_hyperplanes(P_Q),is_pointed(C_Q),(G_Q,c))
end

function ideal(kQ::MonoidAlgebra,gens::Vector{T}) where T <: Union{MPolyRingElem,MPolyQuoRingElem}
    if length(gens) == 0
        return MonoidAlgebraIdeal(kQ,ideal(kQ.algebra,elem_type(base_ring(kQ.algebra))[]))
    end
    for g in gens 
        kQ.algebra == parent(g) || error("Base rings do not match.")
    end
    return MonoidAlgebraIdeal(kQ,ideal(kQ.algebra,gens))
end

function quotient_ring_as_module(I::MonoidAlgebraIdeal)
    return MonoidAlgebraModule(I.monoid_algebra,quotient_ring_as_module(I.ideal))
end

# given a face F of a the cone of a monoid algebra, return the prime k{Q\F}
function get_prime_of_face(k_Q::Union{MPolyRing, MPolyQuoRing},G_Q::Polyhedron,F::Polyhedron)
    gens_PF = []
    for lp in lattice_points(G_Q) 
        if dim(intersect(convex_hull(lp),F))== -1
            a_v = Vector{ZZRingElem}()
            for a_p in lp 
                push!(a_v,a_p)
            end
            push!(gens_PF,monomial_basis(k_Q,a_v)[1])
        end
    end
    return ideal(k_Q,gens_PF)
end


# given a QQ^d vector v, this function returns a vector w that lies on the ray through v and the origin
# INPUT:    rational d-vector 
# OUTPUT:   integer d-vector
function rational_to_integer_vector(v::Union{Vector{QQFieldElem},AbstractVector{Rational}})
    denominators = [denominator(x) for x in v]
    lcm_denominators = lcm(denominators)

    return [Int(numerator(x*lcm_denominators)) for x in v]
end

# given a monoid algebra, this function returns the corresponding polyhedral cone
# INPUT:    monoid algebra k[Q]
# OUTPUT:   polyhedral cone \RR_{\geq 0}Q
function get_polyhedral_cone(R::Union{MPolyDecRing,MPolyQuoRing})
    D = [degree(Vector{Int},g) for g in gens(R)]
    return positive_hull(D)
end

# (works more general for polyhedra)
# given a polyhedral cone, this function returns its faces
# INPUT:    polyhedral cone C
# OUTPUT:   set of faces of C
function get_faces_of_polyhedral_cone(k_Q:: Union{MPolyRing, MPolyQuoRing}, G_Q::Polyhedron,P::Polyhedron)
    P_faces = Vector{Polyhedron}()
    for i=0:dim(P) in 
    append!(P_faces, faces(P,i))
    end
    return [FaceQ(get_prime_of_face(k_Q,G_Q,F),F) for F in P_faces]
end

# given a polyhedral cone, this function returns the hyperplanes bounding it
# INPUT:    polyhedral cone C
# OUTPUT:   set of hyperplanes bounding C
function get_bounding_hyperplanes(P::Polyhedron)
    hyperplanes = Vector{HyperplaneQ}()
    for f in facets(Polyhedron,P)
        hyperplane = polyhedron(affine_hull(f)[1])
        A,b = get_hyperplane_H_presentation(hyperplane)
        push!(hyperplanes,HyperplaneQ(hyperplane,A,b))
    end
    return hyperplanes
end

# given a polyhedral cone, return the zonotope as in Lemma 3.10 in [HM2004]
# INPUT:    polyhedral cone C  
# OUTPUT:   zonotope, sum of primitive integer vector ong rays of C 
function get_zonotope(P::Polyhedron)
    d = ambient_dim(P)
    P_rays = [rational_to_integer_vector(Vector(r)) for r in rays(P)]

    c = zeros(Int,d)
    zonotope = convex_hull(zeros(Int,d))
    for r in P_rays
        zonotope = zonotope + convex_hull([zeros(Int,1,d);reshape(r,1,length(r))])
        c = c + r
    end
    return zonotope,c
end

# given a hyperplane, return the H-presentation of it
# INPUT:    hyperplane h
# OUTPUT:   matrix A, vector b corresponding to Ax \leq b which defines h
function get_hyperplane_H_presentation(h::Polyhedron)
    aff_hull = affine_hull(h).Obj.pm_polytope.AFFINE_HULL
    _M = Matrix{Rational}(aff_hull)
    M = hcat(map(row -> reshape(rational_to_integer_vector(row),1,:),eachrow(_M))...)
    A = [M[:,2:n_columns(M)];-M[:,2:n_columns(M)]]
    b = [M[:,1];-M[:1]]
    return A,b
end

# given an ideal of R and a ZZ^d-vector a, return the module R/I shifted by a
# INPUT:    ideal I,  ZZ^d-vector a
# OUTPUT:   R/I(-a) the module R/I shifted by a
function shifted_module(I::MonoidAlgebraIdeal,shift::Vector)
    S = I.monoid_algebra.algebra
    m_shift = monomial_basis(S,shift)[1]
    I_shifted = ideal(S,[m_shift])*I.ideal
    _M = quotient_ring_as_module(I_shifted)
    return MonoidAlgebraModule(I.monoid_algebra, sub(_M,[m_shift*_M[1]])[1])
end

# computes the generators of k{(a + H_+^°)\cap Q}
# Algorithm 3.11 in HM2004
# INPUT:    MonoidAlgebra kQ
#           HyperplaneQ H that bounding the polyhedral cone RR_{\geq 0}kQ
#           vector a in \ZZ^d
# OUTPUT:   finite set B such that (x^b : b \in B) equals k{(a + H_+^°)\cap Q}
function generators_W_H(kQ::MonoidAlgebra, H::HyperplaneQ, a::Vector{Int})
    F = intersect(H.hyperplane,kQ.cone)

    #get faces of Q intersecting F only at 0 in Q
    D = []
    for f in kQ.faces
        if dim(intersect(f.poly,F)) == 0
            push!(D,f.poly)
        end
    end

    B = []
    PaF = polyhedron(H.A,H.A*a+H.b) #a + RR h
    for d in D
        I = intersect(PaF,d) #(a + RR H)\cap RR_+D
        if dim(I) >= 0
            push!(B,lattice_points(I + kQ.zonotope[1]))
        end
    end

    #eliminate lattice points on (a + \RR H)
    B_red = []
    if length(B) > 0
        for b in B[1]
            if dim(intersect(convex_hull(b),PaF)) < 0
                push!(B_red,b)
            end
        end
    end
    return B_red
end

@doc raw"
    compute_bass_numbers(I,i)

    computes degrees of non-zero Bass numbers up to some given cohomological degree

    INPUT:  monomial ideal I
            integer i
    OUTPUT: list of ZZ^d graded degrees of Bass numbers up to cohomological degree i
"
function compute_bass_numbers(I::Ideal,i::Int)
    R_Q = base_ring(I)
    I_R = modulus(R_Q)
    if is_zero(I_R) #R_Q polynomial ring
        S = R_Q
    else # R_Q quotient ring
        S = base_ring(R_Q)
    end

    # residue_field
    I_m = ideal(S,gens(S))
    R_k,_ = quotient_by_ideal(I_m)

    # module R_Q/I as a subquotient of S^1
    M_mod = quotient_ring_as_module(quo(R_Q,I)[1])
    
    D = Vector{Vector{Int}}()
    for j=0:i in 
        E = Nothing
        try 
            E = ext(R_k,M_mod,j)
        catch e
            if isa(e,AssertionError)
                E = Nothing
            end
        else
            for g in gens(E) 
                push!(D,degree(Vector{Int},g))
            end
        end
    end
    return D
end

##check if points lie in P_Q
##INPUT:    list of points in ZZ^d
##          polyhedron P_Q
##OUTPUT:   true if all points lie in P_Q
function points_in_Q(bass::Vector{Vector{Int}},P_Q::Polyhedron)
    for b in bass
        if !is_subset(convex_hull(b),P_Q)
            return false
        end
    end
    return true
end

##computes a in ZZ^d such that all Bass numbers of M(a) lie in Q
##INPUT:    monomial ideal I
##          integer i up to which cohomological degree Bass numbers are considered
##OUTPUT:   ZZ^d degree of shift
function compute_shift(I::MonoidAlgebraIdeal,i::Int)
    n_bass = compute_bass_numbers(I.ideal,i)
    c = I.monoid_algebra.zonotope[2]
    j = 0
    while !points_in_Q(n_bass,I.monoid_algebra.cone)
        bass_ = [a_bass + c for a_bass in n_bass]
        n_bass = bass_
        j = j +1
    end
    return j*c
end

#computes the quotient of a module M by an ideal I, i.e. (0 :_M I)
#INPUT:     SubquoModule M
#           Ideal I
#OUTPUT:    SubquoModule (0 :_M I)
function mod_quotient(M::SubquoModule,I::Ideal)
    M_I,_ = quotient_by_ideal(I) #base_ring(I)/I as module
    m = M_I[1] #generator of M_I
    H = hom(M_I,M)

    Q_gens = [element_to_homomorphism(g)(m) for g in gens(H[1])]
    return sub(M,Q_gens)
end

#same as mod_quotient except that the result is generated by non-zero elements
function mod_quotient_2(M::SubquoModule,I::Ideal)
    T = elem_type(M)
    M_I,_ = quotient_by_ideal(I) #base_ring(I)/I as module
    m = M_I[1] #generator of M_I
    H = hom(M_I,M)

    Q_gens = Vector{T}()
    for g in gens(H[1])
        h_g = element_to_homomorphism(g)
        b_g = h_g(m)
        if !is_zero(b_g) && !(-b_g in Q_gens)
            push!(Q_gens,b_g)
        end
    end
    return sub(M,Q_gens)
end

#compute saturation of a module M by an ideal, i.e. (0 :_M I^{infty})
#INPUT:     Submodule M
#           Ideal I
#OUTPUT:    SubquoModule (0 :_M I^{infty})
function mod_saturate(M::SubquoModule,I::Ideal)
    M_sat,_ = mod_quotient_2(M,I)
    M_prev = M_sat #previous module quotient
    i = 2
    while true
        M_q,_ = mod_quotient_2(M,I^i)
        if M_prev == M_q
            break
        end
        M_sat = M_sat + M_q
        M_prev = M_q #update
        i = i +1
    end
    return M_sat
end

#Compute k[ZF]-basis as in Algorithm 3.3
#INPUT:     SubquoModule M = (0 :_M P_F) computed with mod_quotient(_,_)
#           prime ideal P_F = k{Q\F}
#OUTPUT:    k[ZF]-basis of localization (0 :_M P_F)[ZF]
function ZF_basis(M::SubquoModule,PF::Ideal)
    T = elem_type(M)
    N = M
    h_N = identity_map(N)
    B = Vector{T}()
    for g in gens(M)
        S = sub(N,[h_N(g)]) #submodule of N =(0 :_M P_F)/(y0,...,yn) generated by g
        if annihilator(S[1]) == PF && !is_zero(S[1])
        #if is_subset(PF,annihilator(S[1]))
            push!(B,g)
        end
        M_B,_ = sub(M,B)
        N,h_N = quo(M,M_B)
        if is_zero(N)
            break
        end
    end
    return filter(!is_zero,B) 
end

#Compute coefficients as in algorithm 3.5
#INPUT:     SubquoModule N
#           SubquoModule Mp = (0 :_M P_F)
#           homomorphism (inclusion) incl_Np: Mxy -> M
#           k[ZF]-basis Bp computed as in ZF_basis(_,_)
#           prime ideal pr
#OUTPUT:    k[ZF]-basis computed as in ZF_basis(_,_)
function coefficients_incl(N::SubquoModule,Np::SubquoModule,incl_Np::SubQuoHom,Bp::Vector,pr::Ideal)
    R_N = base_ring(N)
    if is_zero(modulus(R_N))
        R_poly = R_N
    else
        R_poly = base_ring(R_N)
    end
    T = elem_type(R_N)
    lambda_out = Vector{Vector{T}}()
    p_F = ideal_in_poly_ring(pr)
    i = 1
    for g in gens(N) 
        m_g = monomial_basis(R_N,degree(g))[1]
        N_g,incl_N_g = sub(N,[g;]) #submodule of M generated by g
        I_g,incl_g,incl_p = intersect(Np,N_g) # (0 :_(g) PF)
        lambda_g = Vector{T}()
        j = 1
        for b in Bp 
            lambda_gb = R_N() #zero in R_M
            if length(coordinates(ambient_representative(b))) > 1 #ugly workaround for the special case not considered in algorithm 3.6
                a_yg = degree(b) - degree(g)
                coeff_gb = coordinates(ambient_representative(b))[i] # coefficient of b on m_yg*g
                if is_zero(coeff_gb)
                    lambda_gb = 0*one(R_N)
                else
                    lambda_gb = monomial_basis(parent(coeff_gb),degree(coeff_gb)-a_yg)[1]*one(R_N)
                end
            else
                for g_I in filter(!is_zero,gens(I_g)) 
                    b_m = monomial_basis(R_N,degree(b))[1]
                    g_I_m = monomial_basis(R_N,degree(g_I))[1]
                    if equal_mod_ZF(b_m,g_I_m,p_F)
                        a_yg = degree(b) - degree(g)
                        m_yg = monomial_basis(R_poly,a_yg)[1] #monomial of degree a_yg
                        gg = preimage(incl_Np,m_yg*g) # m_yg*g as an element of Np !not unique
                        lambda_gb = m_g*coordinates(gg)[j] # coefficient of b on m_yg*g
                    end
                end
            end
            push!(lambda_g,lambda_gb) 
            j = j + 1
        end
        push!(lambda_out,lambda_g)
        i = i + 1
    end
    return matrix(R_N,lambda_out) 
end

#Returns ideal of underlying polynomial ring
#INPUT:     ideal I of quotient ring
#OUTPUT:    ideal I as an ideal of polynomial ring
function ideal_in_poly_ring(I::Ideal)
    R_I = base_ring(I)
    J = modulus(R_I)
    if is_zero(J)
        return I
    else
        S = base_ring(R_I)
        G = [monomial_basis(S,degree(g))[1] for g in gens(I)]
        return ideal(S,G)
    end
end

#get quotient ring S/I as module (S quotient ring) 
#alternative to quotient_ring_as_module(I), which gives wrong result when base_ring(I) = quotient ring 
#INPUT:     ideal I 
#OUTPUT:    quotient S/I as SubquoModule (S = base_ring(I))
function quotient_by_ideal(I::Ideal)
    R_I = base_ring(I)
    F = graded_free_module(R_I,1)
    f = F[1]
    V = [g*f for g in gens(I)]
    return quo(F,V)
end


@doc raw"""
    irreducible_hull(M::SubquoModule,P::Vector{FaceQ})

    Computes an irreducible hull of a $\ZZ^d$-graded $k[Q]$-module

    INPUT:    SubquoModule M over semigroup ring k[Q]
               List of prime ideals corresponding to the faces of Q
    OUTPUT:    effective irreducible hull W_bar and effective vector set lambda
"""
function irreducible_hull(M::SubquoModule,P::Vector{FaceQ})
    N = M
    W_F = Vector{FaceQ}()
    W_deg = Vector{Tuple{SubquoModuleElem,Vector{Int}}}()
    lambda = []
    for p in P 
        # Np,incl_Np = mod_quotient(N,p.prime) #compute quotient of N by p
        Np,incl_Np = mod_quotient_2(N,p.prime)
        if is_zero(Np)
            continue
        end
        Bp = ZF_basis(Np,p.prime) #compute ZF-basis of quotient
        Np,incl_Np_n = sub(Np,Bp)
        incl = incl_Np_n*incl_Np
        for b in Bp 
            push!(W_F,p)
            push!(W_deg,(b,degree(Vector{Int},b)))
        end
        lambda_p = coefficients_incl(N,Np,incl,Bp,p.prime)
        push!(lambda,lambda_p)
        N,_ = quo(M,mod_saturate(M,p.prime)) #update N
        if is_zero(N)
            break
        end
    end
    # return W_bar,foldl(hcat,lambda)
    return IrrHull(W_F,W_deg,foldl(hcat,lambda))
end

@doc raw"""
    irreducible_res(M::SubquoModule,P::Vector{FaceQ},P_Q::Polyhedron,G::Polyhedron,F::Vector{Polyhedron{QQFieldElem}},H::Vector{FaceQ})

    Computes a minimal irreducible resolution of a finitely generated $\ZZ^d$-graded module over some affine semigroup ring

    INPUT:  SubquoModule M over semigroup ring $k[Q]$
            list of primes P corresponding to faces
            semigroup Q as polyhedron P_Q
            list of faces F of Q 
            list of hyperplanes H bounding Q as tuples (polyhedron,A,b)
    OUTPUT: irreducible resolution of M given by a list of of k[Q]-modules and homomorphisms in between
"""
# function irreducible_res(M::SubquoModule,P::Vector{FaceQ},P_Q::Polyhedron,G::Polyhedron,F::Vector{Polyhedron{QQFieldElem}},H::Vector{FaceQ})
function irreducible_res(M::MonoidAlgebraModule)
    kQ = M.baseAlgebra
    Mi = M.mod # current module in resolution
    gi = identity_map(Mi) #initilalize
    res_Wi = Vector{IrrSum}()
    res_Mi = [Mi] #cokernels 
    res_hi = Vector{SubQuoHom}()
    res_fi = Vector{SubQuoHom}()
    res_gi = [gi] #quotient maps

    while !is_zero(Mi) #until cokernel Mi is zero
        W_bar = irreducible_hull(Mi,kQ.faces)
        irreducible_ideals = [] #list of irreducible ideals constituting an irreducible hull (sum)
        for k=1: length(W_bar.faces)
            B_i = []
            for h in kQ.hyperplanes
                if is_subset(W_bar.faces[k].poly,h.hyperplane)
                    B_h = generators_W_H(kQ,h,W_bar.vectors[k][2])
                    push!(B_i,B_h)
                end
            end

            G_W = Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
            for b in B_i
                for bb in b
                    a_v = Vector{ZZRingElem}()
                    for a in bb
                        push!(a_v,a)
                    end
                    push!(G_W, monomial_basis(kQ.algebra,a_v)[1])
                end
            end
            push!(irreducible_ideals,ideal(kQ.algebra,G_W))
        end
        
        faces_ = map(p -> p.prime,W_bar.faces)
        vectors = map(x -> x[2],W_bar.vectors)
        indec_injectives = map((x,y) -> IndecInj(x,y),faces_,vectors)

        irreducible_comp = map(I -> quotient_ring_as_module(I),irreducible_ideals)
        d_sum(x,y) = direct_sum(x,y,task=:none)
        Wi = foldl(d_sum,irreducible_comp)
        fi = hom(Mi,Wi,W_bar.lambda)
        hi = gi*fi
        push!(res_Wi,IrrSum(Wi,indec_injectives))
        push!(res_hi,hi)
        push!(res_fi,fi)
        
        Mi_,gi_ = quo(Wi,image(hi)[1]) #cokernel
        
        ### 
        Mi,ji = prune_with_map(Mi_)
        gi = gi_*inv(ji)
        ###
        
        if length(filter(is_zero,relations(Mi))) > 0 # fix modules with "zero" relations
            Mi,h = fix_module(Mi)
            gi = gi*h
        end
        push!(res_Mi,Mi)
        push!(res_gi,gi)
    end
    C = cochain_complex([identity_map(M.mod)]) # default value if sequence not exact
    try
        C = cochain_complex(res_hi)
    catch;
    end
    return IrrRes(res_Wi,res_hi,res_gi,res_fi,res_Mi,C)
end

# compute a minimal injective resolution up to cohomological degree i
# INPUT:    MonoidAlgebraIdeal I
#           positive integer i
# OUTPUT:   injective resolution of R/I up to cohomological degree I
function injective_res(I::MonoidAlgebraIdeal, i::Int)
    kQ = I.monoid_algebra
    a_shift = compute_shift(I,i)
    irrRes = irreducible_res(shifted_module(I,a_shift))
    inj_modules = Vector{InjMod}()
    irr_sums = irrRes.irrSums
    for component in irr_sums 
       shifted_comp = map(indec -> IndecInj(indec.prime,indec.vector - a_shift),component.components)
       push!(inj_modules,InjMod(kQ,shifted_comp))
    end
    cochain_maps = [matrix(cm) for cm in irrRes.cochainMaps]
    cochain_maps[1] = map(a_ij -> kQ.algebra(monomial_basis(kQ.algebra,degree(Vector{Int},a_ij)-a_shift)[1]),cochain_maps[1])
    return InjRes(inj_modules,cochain_maps,i,irrRes,a_shift)
end

##fix SubquoModule with one or more "zero"-relations
function fix_module(M::SubquoModule)
    R = base_ring(M)
    M_F = ambient_free_module(M)
    M_gens = gens(M)
    if length(M_gens) > 0
        M_rels = filter(!is_zero,relations(M))
        d = length(M_gens)
        M_fixed,_ = quo(sub(M_F,M_gens)[1],M_rels)
        return M_fixed,hom(M,M_fixed,identity_matrix(R,d))
    else
        return M,identity_map(M)
    end
end


@doc raw"""
    equal_mod_ZF(a::RingElem,b::RingElem,p_F::Ideal)

check if two ring elements are equal modulo ZF
INPUT:     RingElem a
           RingElem b
           prime ideal p_F 
OUTPUT:    true if a is congruent to b modulo ZF
"""
function equal_mod_ZF(a::RingElem,b::RingElem,p_F::Ideal)
    # @assert parent(a) == parent(b) && parent(a) == base_ring(I) "base ring and parents must match"
    if degree(a) == degree(b)
        return true 
    elseif !(div(a,b) == 0)
        return !ideal_membership(div(a,b),p_F)
    elseif !(div(b,a) == 0)
        return !ideal_membership(div(b,a),p_F)
    else
        return false
    end
end

##computes system of inequalities for polyhedron {r*a : r >= 0}, a \in ZZ^2
##INPUT:    vector a in ZZ^2
##OUTPUT:   matrix A, vector b such that Ax <= b for all x on hyperplane through a
function compute_inequalities(a::Vector{Int})
    num_variables = length(a)
    num_inequalities = 2 * num_variables

    A = zeros(Int, num_inequalities, num_variables)
    b = zeros(Int, num_inequalities)
    if !is_zero(a)
        A[1,1] = -1
        A[1,2] = 0
        A[2,1] = 0
        A[2,2] = -1
        A[3,1] = a[2]
        A[3,2] = -a[1]
        A[4,1] = -a[2]
        A[4,2] = a[1]
    else
        A[1,1] = -1
        A[1,2] = 0
        A[2,1] = 0
        A[2,2] = -1
        A[3,1] = 1
        A[3,2] = 0
        A[4,1] = 0
        A[4,2] = 1
    end
    return A, b
end
    
@doc raw"
    get_face_to_prime(P_F::Ideal)

    Given a homogeneous prime ideal, return the corresponding face F

    INPUT:    prime ideal p_F
    OUTPUT:   list [face F, matrix A, vector b] such that F is defined by Ax <= b
"
function get_face_to_prime(P_F::Ideal)
    R_Q = base_ring(P_F)

    #get generators of k[F] = k[Q]/P_F
    Q_F,_ = quo(R_Q,P_F)
    L_F = [degree(Vector{Int},g) for g in filter(!is_zero, gens(Q_F))]
    
    #get corresponding face defining it by a system of inequalities
    m = length(gens(grading_group(R_Q)))
    C = convex_hull(zeros(Int, m))

    for a in L_F
        A,b = compute_inequalities(a)
        C_a = polyhedron(A,b)
        C = C + C_a
        return FaceQ(P_F,C,A,b)
    end

    #return trivial A,b if L_F is empty
    A = zeros(Int, m,m)
    b = zeros(Int, m)
    return FaceQ(P_F,C,A,b)
end
    
##get all prime ideals P_F associated to M
##INPUT:    ideal I
##OUTPUT:   list of all minimal associated primes + the maximal ideal and the corresponding face F as a polyhedron  
function get_all_ass_primes(I)
    R_Q = base_ring(I)
    I_m = ideal(R_Q,gens(R_Q)) #maximal ideal

    MP = minimal_primes(I) #only minimal primes...
    return [get_face_to_prime(I_m); [get_face_to_prime(mp) for mp in MP]]
end
