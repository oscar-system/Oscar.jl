#injective resolution of a ring k[Q]
function injective_resolution(kQ::MonoidAlgebra,i::Integer)
    I_M = ideal(kQ,[])
    M = quotient_ring_as_module(I_M)
    return Oscar.injective_resolution(M,i)
end

# local cohomology of a ring at maximal ideal
function local_cohomology(kQ::MonoidAlgebra,i::Integer)
    I_M = ideal(kQ,[])
    M = quotient_ring_as_module(I_M)
    m = ideal(kQ,gens(kQ))
    return Oscar.local_cohomology(M,m,i)
end 

# get an ideal over k[Q] as a k[Q]-module
function ideal_as_module(I::MonoidAlgebraIdeal)
    kQ = base_ring(I)
    n = length(gens(I))
    
    I_M = ideal(kQ,[])
    M = quotient_ring_as_module(I_M)

    F,_ = direct_sum([M for i in 1:n])
    G = gens(I)
    return sub(F,[G[i]*F[i] for i in 1:n])[1]
end

# get a_invariant with respect to standard \ZZ grading (total degree)
# how to input some other admissable \ZZ-grading?
function a_invariant(kQ::MonoidAlgebra)
    @req is_normal(kQ) "k[Q] must be normal"
    #compute top local cohomology
    Hd = local_cohomology(kQ,dim(kQ))
    non_zero_sectors = filter(!is_zero, Hd.sectors)
    max = -Inf
    for s in non_zero_sectors
        p = s.sector #get polyhedron
        if lineality_dim(p) <=0
            V = vertices(p)
            i = maximum([sum(v) for v in V])
            if i > max
                max = i
            end
        end
    end
    return max
end

# homology is just local cohomology at zero ideal
function cohomology(M::SubquoModule, i::Integer)
    kQ = base_ring(M)
    Z = ideal(kQ,[])
    return Oscar.local_cohomology(M,Z,i)
end

#Example normal, but not Gorenstein
kQ = monoid_algebra([[1,0],[1,1],[1,2],[1,3]],QQ)
is_normal(kQ)

# top local cohomology
H2 = local_cohomology(kQ,2)
is_zero(H2)
s = filter(!is_zero, H2.sectors) # 10 sectors with non-zero local cohomology
# is this essentially just one sector? 

a = a_invariant(kQ)

# canonical module 
_C = ideal(kQ,[kQ[1],kQ[2]])

# ideal as kQ-module
C = ideal_as_module(_C)

# dualizing complex
D = Oscar.injective_resolution(C,3)

# check vanishing of cohomology module of the dualizing complex
krull_dim(kQ)
is_zero(homology(C,1))
is_zero(homology(C,2))
