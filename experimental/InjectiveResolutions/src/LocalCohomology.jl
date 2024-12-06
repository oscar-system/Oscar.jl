# compute local cohomology modules 
using Oscar

export sector_partition
export _local_cohomology_sectors

# data structures
struct Sector
    A::Vector{Vector{Int}}
    sector::Polyhedron
    indexVector::Vector{Int}
    vectorSpace::Generic.FreeModule{QQFieldElem}
    H::Generic.QuotientModule{QQFieldElem}
end

struct SectorPartition
    sectors::Vector{Sector}
    maps::Vector{Tuple{Sector,Sector,Generic.ModuleHomomorphism}}
end

# definition of monoid algebra as quotient of polynomial ring
S,(x,y,z) = graded_polynomial_ring(QQ,["x","y","z"]; weights = [[0,1],[1,1],[2,1]])
J = ideal(S,[x*z-y^2])
R_Q,phi = quo(S,J)
x,y,z = gens(R_Q)

# get MonoidAlgebra
kQ = get_monoid_algebra(R_Q)

# short example to compute a sector partition
p_0 = kQ.faces[1].prime
F_0 = kQ.faces[1].poly

p_1 = kQ.faces[2].prime
F_1 = kQ.faces[2].poly

p_2 = kQ.faces[3].prime
F_2 = kQ.faces[3].poly

J_0 = IndecInj(p_0,F_0,[2,4])
J_1 = IndecInj(p_1,F_1,[1,4])
J_2 = IndecInj(p_2,F_2,[4,3])
J = [J_0,J_1,J_2]

S = sector_partition(kQ,J...)

# example to a sector partition of a local cohomology
# first steps in an injective resolution of M = k[Q]/(x^4*y,x^2*z), where k[Q] = k[x,y,z]/(xz-y^2)
J_0 = [IndecInj(p_0,F_0,[1,4]), IndecInj(p_1,F_1,[0,4]), IndecInj(p_2,F_2,[1,2])]
J_1 = [IndecInj(p_0,F_0,[1,2]), IndecInj(p_0,F_0,[0,3]), IndecInj(p_0,F_0,[-1,3]), IndecInj(p_1,F_1,[-1,5]), IndecInj(p_2,F_2,[1,0])]
J_2 = [IndecInj(p_0,F_0,[-1,-1]), IndecInj(p_0,F_0,[-1,2]), IndecInj(p_0,F_0,[-2,2])]
J = [J_0...,J_1...,J_2...]
J0 = InjMod(kQ,J_0)
J1 = InjMod(kQ,J_1)
J2 = InjMod(kQ,J_2)
j,k = length(J_0),length(J_0) + length(J_1)

# monomial matrices corresponding to the maps J_0 -> J_1 and J_1 -> J_2
phi_01 = matrix(QQ, [1 1 1 0 0; 0 -1 0 1 0; -1 0 0 0 1])
psi_12 = matrix(QQ, [1 0 0; 0 1 1; -1 -1 -1; 0 1 1; 1 0 0])

I = ideal(kQ,[x,y,z])
H = _local_cohomology_sectors(J0,J1,J2,phi_01,psi_12,I)

# where is the local cohomology module non-zero?
H_nz = filter(h -> dim(h[2]) > 0, H) 

# example for complete computation of local cohomology modules
I_M = ideal(kQ,[x^2*z,x^4*y])
I = ideal(kQ,[x,y,z]) #maximal ideal
i = 2

function local_cohomology(I_M::MonoidAlgebraIdeal,I::MonoidAlgebraIdeal,i::Integer)
    kQ = I_M.monoidAlgebra.algebra
    inj_res = injective_res(I_M,i+1)
    # _phi = matrix(QQ, inj_res.cochainMaps[i])
    # _psi = matrix(QQ, inj_res.cochainMaps[i+1])

    _phi = matrix(QQ,[1 0 0; 0 1 1; -1 -1 -1; 0 1 1; 1 0 0])
    _psi = matrix(QQ,[0 0; -1 -1; 1 1])
    J0 = inj_res.injMods[i]
    J1 = inj_res.injMods[i+1]
    J2 = inj_res.injMods[i+2]

    J,phi,psi,(j,k) = apply_gamma(J0,J1,J2,_phi,_psi,I)
    H = sector_partition(kQ,phi,psi,j,k,J...)
    maps = maps_needed(H)
    return SectorPartition(H,maps)
end

function compute_taus(kQ::MonoidAlgebra, J::IndecInj...)
    # kQ = 
    # get linear functionals tau_i for every hyperplane bounding Q
    A,_ = halfspace_matrix_pair(facets(kQ.cone))
    H = []
    for i=1:length(kQ.hyperplanes) 
        push!(H,[kQ.hyperplanes[i].hyperplane,A[i,:]])
    end

    _tau = []
    for J_i in J in 
        tau_i = []
        for h in H
            if issubset(J_i.face.poly,h[1])
                a_ih = dot(h[2],J_i.vector)
                push!(tau_i,a_ih)
            else
                push!(tau_i,-Inf)
            end
        end
        push!(_tau,tau_i)
    end
    tau = []
    for i=1:length(kQ.hyperplanes) 
        push!(tau,[t[i] for t in _tau])
    end
   return tau
end

function get_halfspace_eq(P::Polyhedron)
    A,_ = halfspace_matrix_pair(facets(P))
    H = []
    for i=1:length(kQ.hyperplanes) 
        push!(H,[kQ.hyperplanes[i].hyperplane,A[i,:]])
    end
    return H
end

function sector_partition(kQ::MonoidAlgebra,phi::QQMatrix,psi::QQMatrix,j::Integer,k::Integer, J::IndecInj...)
    # check that kQ is compatible with the indecomposable injectives
    # ...

    r = length(J)
    # compute tau's
    tau = compute_taus(kQ, J...)

    #sort tau's
    sorted_taus = []
    for t in tau
        pi_t = sortperm(t)
        t_sorted = push!([],-Inf)
        append!(t_sorted,t[pi_t])
        push!(t_sorted,Inf)
        push!(sorted_taus,t_sorted)
    end

    # get vectors corresponding to...
    H = get_halfspace_eq(kQ.cone)

    # compute all stripes
    stripes = []
    for j=1:length(tau) # loop over all linear functionals tau_i
        stripes_j = []
        for i=1:r+1 # loop over all indices
            _A = Vector{Vector{QQFieldElem}}()
            b = Vector{QQFieldElem}()

            if sorted_taus[j][i + 1] == -Inf*one(QQ)
                push!(stripes_j,NaN) # corresponds to empty stripe
                continue
            elseif sorted_taus[j][i + 1] != Inf*one(QQ)
                push!(b,sorted_taus[j][i + 1]-1)
                push!(_A,H[j][2])
            end
            if sorted_taus[j][i] != -Inf*one(QQ)
                push!(b,-sorted_taus[j][i])
                push!(_A,-H[j][2])
            end

            if length(b) < 1
                push!(stripes_j,NaN) # corresponds to empty stripe
                continue
            end
            A = reduce(vcat,transpose(_A))
            push!(stripes_j, polyhedron(A,b))
        end
        push!(stripes,stripes_j)
    end

    # compute the sectors of the sector partition
    # sectors = Vector{Polyhedron}()
    # delta_A = []
    # S_A = []
    S_A = Vector{Sector}()
    for tuple in Iterators.product(ntuple(_ -> 1:r+1,length(tau))...)
        _delta = Vector{Polyhedron}()
        valid = true # is one strip empty? => the intersection is empty
        A_tuple = []
        for i in eachindex(tau)
            if stripes[i][tuple[i]] isa Polyhedron
                push!(_delta,stripes[i][tuple[i]])
            else # in this case the stripe is empty
                valid = false
                break
            end

            # compute the set of all indices j such that the current stripe lies in J_j 
            A_i = Vector{Int}()
            if sorted_taus[i][tuple[i]] == -Inf*one(QQ)
                A_i = findall(x -> x == -Inf*one(QQ),tau[i])
            else
                A_i = findall(x ->  x == -Inf*one(QQ) || (x <= sorted_taus[i][tuple[i]]),tau[i])
            end
            push!(A_tuple,A_i)
        end
        if length(_delta) > 0 && valid
            # compute the sector \Delta(tuple)
            poly_tuple = intersect(_delta...)
            if dim(poly_tuple) < 0
                continue
            elseif dim(poly_tuple) == 0 && length(lattice_points(poly_tuple)) == 0
                continue
            end
            # push!(sectors,poly_tuple)

            # compute the corresponding set A \subseteq {1,...,r}
            # i.e. the indices j such that poly_tuple lies in J_j
            A = intersect(A_tuple...)
            # if !(A in delta_A)
            #     push!(delta_A,A)
            # end

            H_A,A_ = _local_cohomology_sector(A,j,k,phi,psi)
            push!(S_A,Sector(A_,poly_tuple,collect(tuple),free_module(QQ,length(A)),H_A))

            # # add poly_tuple to S_A and add A if necessary
            # if A in [x[1] for x in S_A]
            #     index = findfirst(x -> x[1] == A, S_A)
            #     push!(S_A[index][2],(poly_tuple,tuple))
            # else
            #     H_S = free_module(QQ,length(A))
            #     push!(S_A,(A,[(poly_tuple,tuple)],H_S))
            # end
        end
    end
    return S_A
end

function _local_cohomology_sector(A::Vector{Int},j::Integer,k::Integer,phi::QQMatrix, psi::QQMatrix)
    # divide A into triple
    A_0 = [a for a in A if a <= j]
    A_1 = [a for a in A if j < a <= k]
    A_2 = [a for a in A if k < a]
    _A = [A_0,A_1,A_2]

    J_S0 = free_module(QQ,length(A_0))
    J_S1 = free_module(QQ,length(A_1))
    J_S2 = free_module(QQ,length(A_2))

    #compute the maps by deleting rows and columns in phi and psi
    #phi
    rows = []
    for i in A_0
        push!(rows, phi[i,:])
    end
    _phi_del = transpose(hcat(rows...))

    columns = []
    for i in A_1
        push!(columns, _phi_del[:,i-j])
    end
    if length(columns) > 0
        phi_del = matrix(QQ, hcat(columns...))
    else
        phi_del = matrix(QQ,zeros(QQ,length(A_0),length(A_1)))
    end

    #psi
    rows = []
    for i in A_1
        push!(rows, psi[i-j,:])
    end
    _psi_del = transpose(hcat(rows...))

    columns = []
    for i in A_2
        push!(columns, _psi_del[:,i-k])
    end
    if length(columns) > 0
        psi_del = matrix(QQ,hcat(columns...))
    else 
        psi_del = matrix(QQ,zeros(QQ,length(A_1),length(A_2)))
    end

    f0 = hom(J_S0,J_S1,phi_del)
    f1 = hom(J_S1,J_S2,psi_del)
    return quo(kernel(f1)[1],image(f0)[1])[1], _A
end

function apply_gamma(J0::InjMod,J1::InjMod,J2::InjMod,phi::QQMatrix, psi::QQMatrix,I::MonoidAlgebraIdeal)
    @req J0.monoidAlgebra == J1.monoidAlgebra == J2.monoidAlgebra == I.monoidAlgebra "Input not over same monoid algebra!"
    J_0 = J0.indecInjectives
    J_1 = J1.indecInjectives
    J_2 = J2.indecInjectives

    _J_0 = []
    rows_phi = []
    for i in eachindex(J_0)
        if issubset(I.ideal,J_0[i].face.prime) 
            push!(_J_0,J_0[i])
            push!(rows_phi,phi[i,:])
        end
    end
    _phi = transpose(hcat(rows_phi...))

    columns_phi = []
    rows_psi = []
    _J_1 = []
    for i in eachindex(J_1)
        if issubset(I.ideal,J_1[i].face.prime) 
            push!(_J_1,J_1[i])
            push!(columns_phi,_phi[:,i] )
            push!(rows_psi,psi[i,:])
        end
    end
    phi = matrix(QQ,hcat(columns_phi...))
    _psi = transpose(hcat(rows_psi...))

    columns_psi = []
    _J_2 = []
    for i in eachindex(J_2)
        if issubset(I.ideal,J_2[i].face.prime) 
            push!(_J_2,J_2[i])
            push!(columns_psi, _psi[:,i])
        end
    end
    psi = matrix(QQ,hcat(columns_psi...))

    return [_J_0...,_J_1...,_J_2...],phi,psi,(length(_J_0), length(_J_0) + length(_J_1))
end

function maps_needed(S_A::Vector{Sector})
    # # get a list of all sectors and the corresponding set A and vector space H_A
    # _S = []
    # for x in S_A
    #     _Sx = [(x[1][1],y,x[1][3]) for y in x[1][2]]
    #     append!(_S,_Sx)
    # end

    maps = []
    # for s1 in _S, s2 in _S # loop over all pairs
    for s1 in S_A, s2 in S_A
        if s1 == s2
            continue
        end
        if issubset(s2.A,s1.A) && s1.indexVector >= s2.indexVector # check if B \subseteq A and (l_1, ..., l_n) <= (l_1', ..., l_n')
            K = s2.sector + (-1)*s1.sector # minkowski sum (|Delta_B - \Delta_A)
            if dim(intersect(K,kQ.cone)) > -1 # check if Q \cap (|Delta_B - \Delta_A) \neq \emptyset
                # compute the map x^{B - A}: H_A -> H_B
                A = vcat(s1.A...)
                B = vcat(s2.A...)
                H_A = s1.vectorSpace
                H_B = s2.vectorSpace
                T = elem_type(H_A)
                V = Vector{T}()
                for i in eachindex(A)
                    if A[i] in B
                        j = findfirst(x -> x == A[i], B)
                        push!(V,H_B[j])
                    else
                        push!(V,0*zero(H_B))
                    end
                end

                # define map and push
                f_AB = hom(H_A,H_B,V)

                if is_zero(f_AB)
                    continue
                else
                    push!(maps,(s1,s2,f_AB))
                end
            end
        else 
            continue
        end
    end
    return maps
end

# (([[1], [2, 3], Int64[]], (Polyhedron in ambient dimension 2, (4, 8)), Vector space of dimension 3 over QQ), ([[1], Int64[], Int64[]], (Polytope in ambient dimension 2, (3, 3)), Vector space of dimension 1 over QQ))
# s1 = ([[1], [2, 3], Int64[]], (Polyhedron in ambient dimension 2, (4, 8)), Vector space of dimension 3 over QQ)
# -> _S[23]
# s2 = ([[1], Int64[], Int64[]], (Polytope in ambient dimension 2, (3, 3)), Vector space of dimension 1 over QQ)
# -> _S[11]

function _local_cohomology_sectors(J0::InjMod,J1::InjMod,J2::InjMod,_phi::QQMatrix, _psi::QQMatrix,I::MonoidAlgebraIdeal)
    J,phi,psi,(j,k) = apply_gamma(J0,J1,J2,_phi,_psi,I)
    S = sector_partition(I.monoidAlgebra, J...)

    S_A = []
    for s in S[1]
        # divide A into triple
        A_0 = [a for a in s[1] if a <= j]
        A_1 = [a for a in s[1] if j < a <= k]
        A_2 = [a for a in s[1] if k < a]
        A = [A_0,A_1,A_2]

        J_S0 = free_module(QQ,length(A_0))
        J_S1 = free_module(QQ,length(A_1))
        J_S2 = free_module(QQ,length(A_2))

        #compute the maps by deleting rows and columns in phi and psi
        #phi
        rows = []
        for i in A_0
            push!(rows, phi[i,:])
        end
        _phi_del = transpose(hcat(rows...))

        columns = []
        for i in A_1
            push!(columns, _phi_del[:,i-j])
        end
        if length(columns) > 0
            phi_del = matrix(QQ, hcat(columns...))
        else
            phi_del = matrix(QQ,zeros(QQ,length(A_0),length(A_1)))
        end

        #psi
        rows = []
        for i in A_1
            push!(rows, psi[i-j,:])
        end
        _psi_del = transpose(hcat(rows...))

        columns = []
        for i in A_2
            push!(columns, _psi_del[:,i-k])
        end
        if length(columns) > 0
            psi_del = matrix(QQ,hcat(columns...))
        else 
            psi_del = matrix(QQ,zeros(QQ,length(A_1),length(A_2)))
        end

        f0 = hom(J_S0,J_S1,phi_del)
        f1 = hom(J_S1,J_S2,psi_del)
        H_012 = quo(kernel(f1)[1],image(f0)[1])[1]
        push!(S_A,([A,s[2],s[3]], H_012))
    end
    return S_A
end

function _sector_partition(kQ::MonoidAlgebra, J::IndecInj...)
    # check that kQ is compatible with the indecomposable injectives
    # ...

    r = length(J)
    # compute tau's
    tau = compute_taus(kQ, J...)

    #sort tau's
    sorted_taus = []
    for t in tau
        pi_t = sortperm(t)
        t_sorted = push!([],-Inf)
        append!(t_sorted,t[pi_t])
        push!(t_sorted,Inf)
        push!(sorted_taus,t_sorted)
    end

    # get vectors corresponding to...
    H = get_halfspace_eq(kQ.cone)

    # compute all stripes
    stripes = []
    for j=1:length(tau) # loop over all linear functionals tau_i
        stripes_j = []
        for i=1:r+1 # loop over all indices
            _A = Vector{Vector{QQFieldElem}}()
            b = Vector{QQFieldElem}()

            if sorted_taus[j][i + 1] == -Inf*one(QQ)
                push!(stripes_j,NaN) # corresponds to empty stripe
                continue
            elseif sorted_taus[j][i + 1] != Inf*one(QQ)
                push!(b,sorted_taus[j][i + 1]-1)
                push!(_A,H[j][2])
            end
            if sorted_taus[j][i] != -Inf*one(QQ)
                push!(b,-sorted_taus[j][i])
                push!(_A,-H[j][2])
            end

            if length(b) < 1
                push!(stripes_j,NaN) # corresponds to empty stripe
                continue
            end
            A = reduce(vcat,transpose(_A))
            push!(stripes_j, polyhedron(A,b))
        end
        push!(stripes,stripes_j)
    end

    # compute the sectors of the sector partition
    sectors = Vector{Polyhedron}()
    delta_A = []
    S_A = []
    for tuple in Iterators.product(ntuple(_ -> 1:r+1,length(tau))...)
        _delta = Vector{Polyhedron}()
        valid = true # is one strip empty? => the intersection is empty
        A_tuple = []
        for i in eachindex(tau)
            if stripes[i][tuple[i]] isa Polyhedron
                push!(_delta,stripes[i][tuple[i]])
            else # in this case the stripe is empty
                valid = false
                break
            end

            # compute the set of all indices j such that the current stripe lies in J_j 
            A_i = Vector{Int}()
            if sorted_taus[i][tuple[i]] == -Inf*one(QQ)
                A_i = findall(x -> x == -Inf*one(QQ),tau[i])
            else
                A_i = findall(x ->  x == -Inf*one(QQ) || (x <= sorted_taus[i][tuple[i]]),tau[i])
            end
            push!(A_tuple,A_i)
        end
        if length(_delta) > 0 && valid
            # compute the sector \Delta(tuple)
            poly_tuple = intersect(_delta...)
            if dim(poly_tuple) < 0
                continue
            elseif dim(poly_tuple) == 0 && length(lattice_points(poly_tuple)) == 0
                continue
            end
            push!(sectors,poly_tuple)

            # compute the corresponding set A \subseteq {1,...,r}
            # i.e. the indices j such that poly_tuple lies in J_j
            A = intersect(A_tuple...)
            if !(A in delta_A)
                push!(delta_A,A)
            end

            # add poly_tuple to S_A and add A if necessary
            if A in [x[1] for x in S_A]
                index = findfirst(x -> x[1] == A, S_A)
                push!(S_A[index][2],(poly_tuple,tuple))
            else
                H_S = free_module(QQ,length(A))
                push!(S_A,(A,[(poly_tuple,tuple)],H_S))
            end
        end
    end
    return S_A, sectors, delta_A
end

function _maps_needed(S_A)
    # get a list of all sectors and the corresponding set A and vector space H_A
    _S = []
    for x in S_A
        _Sx = [(x[1][1],y,x[1][3]) for y in x[1][2]]
        append!(_S,_Sx)
    end

    maps = []
    for s1 in _S, s2 in _S # loop over all pairs
        if s1 == s2
            continue
        end
        if issubset(s2[1],s1[1]) && s1[2][2] >= s2[2][2] # check if B \subseteq A and (l_1, ..., l_n) <= (l_1', ..., l_n')
            K = s2[2][1] + (-1)*s1[2][1] # minkowski sum (|Delta_B - \Delta_A)
            if dim(intersect(K,kQ.cone)) > -1 # check if Q \cap (|Delta_B - \Delta_A) \neq \emptyset
                # compute the map x^{B - A}: H_A -> H_B
                A = vcat(s1[1]...)
                B = vcat(s2[1]...)
                H_A = s1[3]
                H_B = s2[3]
                T = elem_type(H_A)
                V = Vector{T}()
                for i in eachindex(A)
                    if A[i] in B
                        j = findfirst(x -> x == A[i], B)
                        push!(V,H_B[j])
                    else
                        push!(V,0*zero(H_B))
                    end
                end

                # define map and push
                f_AB = hom(H_A,H_B,V)

                if is_zero(f_AB)
                    continue
                else
                    push!(maps,(s1,s2,f_AB))
                end
            end
        else 
            continue
        end
    end
    return maps
end
#######################################
## first try in computing a sector partition
# get halfspaces defining Q, i.e. the linear functionals \tau_i
Q_halfspaces = facets(kQ.cone)
_H1 = kQ.hyperplanes[1].hyperplane
_H2 = kQ.hyperplanes[2].hyperplane

# get linear functional as vector
A,b = halfspace_matrix_pair(facets(kQ.cone))
H1 = [_H1,A[1,:]]
H2 = [_H2,A[2,:]]
H = [H1,H2]

# get faces of Q and set a fixed order
F = kQ.faces
F_empty = F[1].poly
F_1 = F[2].poly
F_2 = F[3].poly

# define indecomposable injectives
J_0 = [[2,4],F_empty]
J_1 = [[1,4],F_1]
J_2 = [[4,3],F_2]

J0 = IndecInj(F[1].prime,F_empty,[2,4])
J1 = IndecInj(F[2].prime,F_1,[1,4])
J2 = IndecInj(F[3].prime,F_2,[4,3])

J = [J_0,J_1,J_2]
J = [J0,J1,J2]

tau = []
for J_i in J
    tau_i = []
    for h in H
        if issubset(J_i[2],h[1])
            a_ih = dot(h[2],J_i[1])
            push!(tau_i,a_ih)
        else
            push!(tau_i,-Inf)
        end
    end
    push!(tau,tau_i)
end

# sort the k-th components
tau1 = [t[1] for t in tau]
tau2 = [t[2] for t in tau]

pi1 = sortperm(tau1)

pi2 = sortperm(tau2)

tau1_sorted = push!([],-Inf)
append!(tau1_sorted,tau1[pi1])
push!(tau1_sorted,Inf)

tau2_sorted = push!([],-Inf)
append!(tau2_sorted,tau2[pi2])
push!(tau2_sorted,Inf)

function ttau(beta)
    return (dot(H1[2],beta),dot(H2[2],beta))
end

ttau([5,6])
ttau([3,2])
ttau([5,6]) <= ttau([3,2])

delta = []
delta_A = []
S_A = []
for l1 in 1:4
    _A1 = Vector{Vector{QQFieldElem}}()
    b1 = Vector{QQFieldElem}()
    if tau1_sorted[l1 + 1] == -Inf*one(QQ)
        continue
    end
    if tau1_sorted[l1 + 1] != Inf*one(QQ)
        push!(b1, tau1_sorted[l1 + 1]-1)
        push!(_A1,H1[2])
    end
    if tau1_sorted[l1] != -Inf*one(QQ)
        push!(b1,-tau1_sorted[l1])
        push!(_A1,-H1[2])
    end
    if length(b1) < 1
        continue
    end
    A1 = reduce(vcat,transpose.(_A1))
    p1 = polyhedron(A1,b1)
    for l2 in 1:4 
        _A2 = Vector{Vector{QQFieldElem}}()
        b2 = Vector{QQFieldElem}()
        if tau2_sorted[l2 + 1] == -Inf*one(QQ)
            continue
        end
        if tau2_sorted[l2 + 1] != Inf*one(QQ)
            push!(b2,tau2_sorted[l2 + 1]-1)
            push!(_A2,H2[2])
        end
        if tau2_sorted[l2] != -Inf*one(QQ)
            push!(b2,-tau2_sorted[l2])
            push!(_A2,-H2[2])
        end
        if length(b2) < 1
            continue
        end
        A2 = reduce(vcat,transpose.(_A2))
        p2 = polyhedron(A2,b2)
        p_12 = intersect(p1,p2)
        push!(delta,p_12)

        I1 = Vector{Int}()
        I2 = Vector{Int}()
        if tau1_sorted[l1] == -Inf*one(QQ)
            I1 = findall(x -> x == -Inf*one(QQ),tau1)
        else
            I1 = findall(x ->  x == -Inf*one(QQ) || (x <= tau1_sorted[l1]),tau1)
        end

        if tau2_sorted[l2] == -Inf*one(QQ)
            I2 = findall(x -> x == -Inf*one(QQ),tau2)
        else
            I2 = findall(x -> x == -Inf*one(QQ) || (x <= tau2_sorted[l2]),tau2)
        end
        A = intersect(I1,I2)
        if !(A in delta_A)
            push!(delta_A,A)
        end
        if A in [x[1] for x in S_A]
            index = findfirst(x -> x[1] == A, S_A)
            push!(S_A[index][2],(p_12,(l1,l2)))
        else
            push!(S_A,(intersect(I1,I2),[(p_12,(l1,l2))]))
        end
    end
end

[x[2] for x in S_A]

_S = []
for x in S_A
    _Sx = [(x[1],y) for y in x[2]]
    append!(_S,_Sx)
end

true_values = []
for s1 in _S, s2 in _S
    if s1 == s2
        continue
    end
    if issubset(s2[1],s1[1]) && s1[2][2] >= s2[2][2]
        K = s2[2][1] + (-1)*s1[2][1]
        if dim(intersect(K,kQ.cone)) > -1
            push!(true_values,(s1,s2))
        end
    else 
        continue
    end
end

[dim(x) for x in delta]

########################################
# some more tries on applying \Gamma_I and computing local cohomology 

S_A = []
for s in S[1]
    # divide A into triple
    A_0 = [a for a in s[1] if a <= j]
    A_1 = [a for a in s[1] if j < a <= k]
    A_2 = [a for a in s[1] if k < a]
    A = [A_0,A_1,A_2]

    J_S0 = free_module(QQ,length(A_0))
    J_S1 = free_module(QQ,length(A_1))
    J_S2 = free_module(QQ,length(A_2))

    #compute the maps by deleting rows and columns in phi and psi
    #phi
    rows = []
    for i in A_0
        push!(rows, phi_01[i,:])
    end
    _phi_del = transpose(hcat(rows...))

    columns = []
    for i in A_1
        push!(columns, _phi_del[:,i-j])
    end
    if length(columns) > 0
        phi_del = matrix(QQ, hcat(columns...))
    else
        phi_del = matrix(QQ,zeros(QQ,length(A_0),length(A_1)))
    end

    #psi
    rows = []
    for i in A_1
        push!(rows, psi_12[i-j,:])
    end
    _psi_del = transpose(hcat(rows...))

    columns = []
    for i in A_2
        push!(columns, _psi_del[:,i-k])
    end
    if length(columns) > 0
        psi_del = matrix(QQ,hcat(columns...))
    else 
        psi_del = matrix(QQ,zeros(QQ,length(A_1),length(A_2)))
    end

    f0 = hom(J_S0,J_S1,phi_del)
    f1 = hom(J_S1,J_S2,psi_del)
    H_012 = quo(kernel(f1)[1],image(f0)[1])[1]
    push!(S_A,([A,s[2],s[3]], H_012))
end

H = [s[2] for s in S_A]
filter(x -> dim(x) > 0, H)

# apply \Gamma_I
I = p_0 # maximal ideal

_J_0 = []
rows_phi = []
for i in eachindex(J_0)
    if issubset(I,J_0[i].prime) 
        push!(_J_0,J_0[i])
        push!(rows_phi,phi_01[i,:])
    end
end
_phi = transpose(hcat(rows_phi...))

columns_phi = []
rows_psi = []
_J_1 = []
for i in eachindex(J_1)
    if issubset(I,J_1[i].prime) 
        push!(_J_1,J_1[i])
        push!(columns_phi,_phi[:,i] )
        push!(rows_psi,psi_12[i,:])
    end
end
phi = matrix(QQ,transpose(hcat(columns_phi...)))
_psi = transpose(hcat(rows_psi...))

columns_psi = []
_J_2 = []
for i in eachindex(J_2)
    if issubset(I,J_2[i].prime) 
        push!(_J_2,J_2[i])
        push!(columns_psi, _psi[:,i])
    end
end
psi = matrix(QQ,transpose(hcat(columns_psi...)))

_J = [_J_0...,_J_1...,_J_2...]
j,k = length(_J_0), length(_J_0) + length(_J_1)


