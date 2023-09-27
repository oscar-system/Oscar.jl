import Oscar.gens, AbstractAlgebra.direct_sum, Oscar.invariant_ring
#export ReductiveGroup, reductive_group, representation_matrix, group, reynolds_operator, group_ideal, canonical_representation, natural_representation
#export InvariantRing, invariant_ring, gens, hilbert_ideal, derksen_ideal
##########################
#Reductive Groups
##########################
mutable struct ReductiveGroup
    group::Tuple{Symbol, Int}
    direct_sum::Tuple{Bool, Vector{Int}}
    direct_product::Tuple{Bool,Int}
    vector_space_dimension::Int
    rep_mat::AbstractAlgebra.Generic.MatSpaceElem
    group_ideal::MPolyIdeal
    reynolds_operator::Function
    canonical_representation::AbstractAlgebra.Generic.MatSpaceElem
    
    function ReductiveGroup(sym::Symbol, m::Int) #have not decided the representation yet
        G = new()
        if sym != :SL
            error("Only implemented for SLm")
        end
        G.group = (sym,m)
        G.reynolds_operator = reynolds__
        R,_ = PolynomialRing(QQ,"z"=> 1:m^2)
        M = matrix(R,m,m,[0 for i in 1:m^2])
        for i in 1:m
            M[1:m,i] = gens(R)[m*(i-1)+1:i*m]
        end
        G.canonical_representation = M
        #base ring of M has to be the same as the representation matrix when that is created later.
        return G
    end
    
    function ReductiveGroup(sym::Symbol, m::Int, sym_deg::Int)
        G = new()
        if sym != :SL
            error("Only implemented for SLm")
        end
        G.group = (sym, m)
        G.rep_mat = rep_mat_(m,sym_deg)
        G.vector_space_dimension = ncols(G.rep_mat)
        R = base_ring(G.rep_mat)
        M = matrix(R,m,m,[0 for i in 1:m^2])
        for i in 1:m
            M[1:m,i] = gens(R)[m*(i-1)+1:i*m]
        end
        G.group_ideal = ideal([det(M) - 1])
        G.reynolds_operator = reynolds__
        G.canonical_representation = M
        G.direct_sum = (false,[])
        G.direct_product = (false,1)
        return G
    end
    
    #direct sum
    function ReductiveGroup(sym::Symbol, m::Int, v::Vector{Int})
        G = new()
        if sym != :SL
            error("Only implemented for SLm")
        end
        G.group = (sym, m)
        G.direct_sum = (true,v)
        G.direct_product = (false,1)
        G.rep_mat = rep_mat_(m,v,true)
        G.vector_space_dimension = nrows(G.rep_mat)
        R = base_ring(G.rep_mat)
        M = matrix(R,m,m,[0 for i in 1:m^2])
        for i in 1:m
            M[1:m,i] = gens(R)[m*(i-1)+1:i*m]
        end
        G.group_ideal = ideal([det(M) - 1])
        G.canonical_representation = M
        G.reynolds_operator = reynolds__
        return G
    end
    
    function ReductiveGroup(sym::Symbol, m::Int, r::Int, direct_product::Bool)
        if direct_product
        G = new()
        G.group = (sym,m)
        v = [1 for i in 1:r]
        G.rep_mat = rep_mat_(m,v,false)
            G.vector_space_dimension = ncols(G.rep_mat)
        G.direct_sum = (false,[])
        G.direct_product = (true,r)
        R = base_ring(G.rep_mat)
        G.reynolds_operator = reynolds__
        M = matrix(R,r*m, r*m, [0 for i in 1:r^2*m^2])
        x = 1
        for k in 1:r
            for i in ((k-1)*m)+1:k*m
                for j in ((k-1)*m)+1:k*m
                    M[i,j] = gens(R)[x]
                    x+=1
                end
            end
        end
        G.canonical_representation = M
            v = Vector{elem_type(R)}()
            x = 1
            for i in 1:r
                mat = matrix(R,m,m,[0 for i in 1:m^2])
                for i in 1:m
                    for j in 1:m
                        mat[i,j] = gens(R)[x]
                        x+=1
                    end
                end
                push!(v, det(mat) - 1)
            end
            G.group_ideal = ideal(v)
        return G
        else
            return ReductiveGroup(sym,m,r)
        end
    end
end

function Base.show(io::IO, G::ReductiveGroup)
    if !G.direct_sum[1] && !G.direct_product[1]
        println(io, "Reductive group ", G.group[1], G.group[2])
        println(io, "acting on vector space of dimension ", G.vector_space_dimension)
    elseif G.direct_sum[1]
        println(io, "Reductive group ", G.group[1], G.group[2])
        println(io, "acting on direct sum of vector spaces of dimensions ", G.direct_sum[2])
    end
end

function reductive_group(sym::Symbol, m::Int, n::Int)
    return ReductiveGroup(sym,m,n)
end

function reductive_group(sym::Symbol, m::Int, v::Vector{Int})
    return ReductiveGroup(sym,m,v)
end

function reductive_group(sym::Symbol, m::Int, r::Int, dir_prod::Bool)
    return ReductiveGroup(sym,m,r,dir_prod)
end

representation_matrix(G::ReductiveGroup) = G.rep_mat
group(G::ReductiveGroup) = G.group
reynolds_operator(G::ReductiveGroup) = G.reynolds_operator
group_ideal(G::ReductiveGroup) = G.group_ideal
canonical_representation(G::ReductiveGroup) = G.canonical_representation
natural_representation(G::ReductiveGroup) = G.canonical_representation

###############

#for the representation matrices

function rep_mat_(m::Int, sym_deg::Int)
    n = binomial(m + sym_deg - 1, m - 1)
    mixed_ring, t, z, a = PolynomialRing(QQ, "t"=> 1:m, "z"=> (1:m, 1:m), "a" => 1:n)
    group_mat = matrix(mixed_ring, m,m,[z[i,j] for i in 1:m, j in 1:m])
    vars = [t[i] for i in 1:m]
    new_vars = vars*group_mat
    degree_basiss = degree_basis(mixed_ring,m, sym_deg)
    sum = mixed_ring()
    for j in 1:length(degree_basiss)
        prod = mixed_ring(1)
        factors_ = factorise(degree_basiss[j])
        for i in 1:length((factors_)[1:m])
            prod = prod*new_vars[i]^(factors_[i][2])
        end
        prod = prod*a[j]
        sum += prod
    end
    mons = collect(monomials(sum))
    coeffs = collect(coefficients(sum))
    mat = matrix(mixed_ring, n,n,[0 for i in 1:n^2])
    for i in 1:n
        for j in 1:n
            for k in 1:length(mons)
                if divides(mons[k], a[j])[1] && divides(mons[k], degree_basiss[i])[1]
                    mat[i,j] += coeffs[k]*numerator((mons[k]//degree_basiss[i])//a[j])
                end
            end
        end
    end
    #we have to return mat in a different ring! 
    group_ring, Z = PolynomialRing(QQ, "Z"=>(1:m, 1:m))
    mapp = hom(mixed_ring, group_ring, vcat([0 for i in 1:m], gens(group_ring), [0 for i in 1:n]))
    Mat = matrix(group_ring, n, n, [mapp(mat[i,j]) for i in 1:n, j in 1:n])
    return Mat            
end

function degree_basis(R::MPolyRing,m::Int, t::Int)
    v = R(1)
    genss = gens(R)
    n = length(gens(R)[1:m])
    C = zero_matrix(Int,n,n)
    for i in 1:n
      C[i,i] = -1 #find a better way to write this TODO
    end
    d = [0 for i in 1:n]
    A = [1 for i in 1:n]
    b = [t]
    P = Polyhedron((C,d),(A,b))
    L = lattice_points(P)
    W = Vector{MPolyRingElem}(undef,0)
    for l in L
      for i in 1:n
        v = v*genss[i]^l[i]
      end
    push!(W,v)
    v = R(1)
    end
    return W
end


#############################
#DIRECT SUM


function rep_mat_(m::Int, v::Vector{Int}, dir_sum::Bool)
    if dir_sum
        n = num_of_as(m,v)
        mixed_ring, t, z, a = PolynomialRing(QQ, "t"=> 1:m, "z"=> (1:m, 1:m), "a" => 1:n)
        group_mat = matrix(mixed_ring, m,m,[z[i,j] for i in 1:m, j in 1:m])
        vars = [t[i] for i in 1:m]
        new_vars = vars*group_mat
        big_matrix = matrix(mixed_ring, 0,0,[])
        for sym_deg in v
            degree_basiss = degree_basis(mixed_ring,m, sym_deg)
            sum = mixed_ring()
            for j in 1:length(degree_basiss)
                prod = mixed_ring(1)
                factors_ = factorise(degree_basiss[j])
                for i in 1:length((factors_)[1:m])
                    prod = prod*new_vars[i]^(factors_[i][2])
                end
                prod = prod*a[j]
                sum += prod
            end
            mons = collect(monomials(sum))
            coeffs = collect(coefficients(sum))
            N = binomial(m + sym_deg - 1, m - 1)
            mat = matrix(mixed_ring, N,N,[0 for i in 1:N^2])
            for i in 1:N
                for j in 1:N
                    for k in 1:length(mons)
                        if divides(mons[k], a[j])[1] && divides(mons[k], degree_basiss[i])[1]
                            mat[i,j] += coeffs[k]*numerator((mons[k]//degree_basiss[i])//a[j])
                        end
                    end
                end
            end
            big_matrix = direct_sum(big_matrix,mat)
        end
        #we have to return mat in a different ring! 
        group_ring, Z = PolynomialRing(QQ, "Z"=>(1:m, 1:m))
        mapp = hom(mixed_ring, group_ring, vcat([0 for i in 1:m], gens(group_ring), [0 for i in 1:n]))
        Mat = matrix(group_ring, nrows(big_matrix), nrows(big_matrix), [mapp(big_matrix[i,j]) for i in 1:nrows(big_matrix), j in 1:nrows(big_matrix)])
        return Mat
    else
        l = length(v)
        #there are l*m^2 z's
        ringg,Z = PolynomialRing(QQ,"Z"=>1:l*m^2)
        prod_mat = matrix(ringg,1,1,[1])
        for i in 1:l
            mat = rep_mat_(m,v[i])
            R = base_ring(mat)
            mapp = hom(R,ringg, gens(ringg)[((i-1)*m^2)+1:i*m^2])
            mat_ = matrix(ringg, nrows(mat), ncols(mat), [0 for i in 1:nrows(mat)*ncols(mat)])
            for i in 1:nrows(mat)
                for j in 1:ncols(mat)
                    mat_[i,j] = mapp(mat[i,j])
                end
            end
            prod_mat = kronecker_product(prod_mat, mat_)
        end
        return prod_mat
    end
end

function direct_sum(M::AbstractAlgebra.Generic.MatSpaceElem, N::AbstractAlgebra.Generic.MatSpaceElem)
    #parent(M[1,1]) == parent(N[1,1]) || @error("not same ring")
    mr = nrows(M)
    mc = ncols(M)
    nr = nrows(N)
    nc = ncols(N)
    mat_ = matrix(base_ring(M), mr+nr, mc+nc, [0 for i in 1:(mr+nr)*(mc+nc)])
    mat_[1:mr, 1:mc] = M
    mat_[mr+1:mr+nr, mc+1:mc+nc] = N
    return mat_
end

function num_of_as(m::Int, v::Vector{Int})
    max = maximum(v)
    sum = binomial(m + max - 1, m - 1)
    return sum
end

##########################
#Invariant Rings of Reductive groups
##########################
mutable struct InvariantRing
    poly_ring::MPolyDecRing #graded
    group::ReductiveGroup
    generators::Vector{MPolyDecRingElem}
    HilbertIdeal::MPolyIdeal
    DerksenIdeal::MPolyIdeal
    cached_rep_mat::AbstractAlgebra.Generic.MatSpaceElem

    #Not finished
    function InvariantRing(G::ReductiveGroup, rep_mat::AbstractAlgebra.Generic.MatSpaceElem)
        G.group[1] != SL && return nothing
        z = new()
        R = base_ring(rep_mat)
        base_ring(R) == QQ || error("must be rational field")
        m = Int(sqrt(ngens(R)))
        n = ncols(rep_mat)
        
        M = matrix(R,1,m,gens(R)[1:m])
        for i in 1:m-1 
            M = vcat(M, matrix(R,1,m,gens(R)[(i)*m+1:(i+1)*m]))
        end
        det_ = det(M)
        
        G.rep_mat = rep_mat
        G.group_ideal = det(M) - 1
        z.poly_ring, _ = grade(PolynomialRing(QQ, "X" => 1:n)[1])
        z.group = G
        z.generators = inv_generators(G, z.poly_ring)
        return z
    end
    
    #This should be the only one that is used
    function InvariantRing(G::ReductiveGroup) #here G already contains information n and rep_mat
        z = new()
        n = G.vector_space_dimension
        z.poly_ring, __ = grade(PolynomialRing(QQ, "X" => 1:n)[1])
        z.group = G
        #right now these will be computed 2 or 3 times. TODO.
        z.DerksenIdeal, z.cached_rep_mat = proj_of_image_ideal(G)
        z.HilbertIdeal = ideal(generators(G, z.DerksenIdeal))
        z.generators = inv_generators(G, z.poly_ring, z.cached_rep_mat, z.DerksenIdeal)
        return z
    end
end

invariant_ring(G::ReductiveGroup) = InvariantRing(G)
gens(R::InvariantRing) = R.generators
hilbert_ideal(R::InvariantRing) = R.HilbertIdeal
derksen_ideal(R::InvariantRing) = R.DerksenIdeal

function Base.show(io::IO, R::InvariantRing) #TODO compact printing
    #if get(io, :supercompact, false)
        print(io, "Invariant Ring of", "\n")
        show(io, R.poly_ring)
        print(io, " under group action of ", R.group.group[1], R.group.group[2], "\n", "\n")
        print(io, "Generated by ", "\n")
    join(io, R.generators, "\n")
    #else
     #   print(io, "Invariant Ring under ")
      #  show(io, R.group[1], R.group[2])
    #end
end

#I think this should work for everyone except tensor products.
function image_ideal(G::ReductiveGroup)
    rep_mat = G.rep_mat
    R = base_ring(rep_mat)
    n = G.vector_space_dimension
    m = G.group[2]
    r = G.direct_product[2]
    mixed_ring_xy, x, y, zz = PolynomialRing(QQ, "x"=>1:n, "y"=>1:n, "zz"=>1:r*m^2)
    # naming the variables zz here and z in the group ring
    #for determinant - 
    #M1 = matrix(mixed_ring_xy,m,m,[zz[i,j] for i in 1:m for j in 1:m])
    ztozz = hom(R,mixed_ring_xy, gens(mixed_ring_xy)[(2*n)+1:(2*n)+(r*m^2)])
    genss_ = gens(G.group_ideal)
    genss = [ztozz(genss_[i]) for i in 1:length(genss_)]
    #rep_mat in the new ring
    new_rep_mat = matrix(mixed_ring_xy,n,n,[ztozz(rep_mat[i,j]) for i in 1:n, j in 1:n])
    new_vars = new_rep_mat*[x[i] for i in 1:n] #changed the order - didn't help. 
    ideal_vect = [y[i] - new_vars[i] for i in 1:n] #check here TODO
    ideal_vect = vcat(ideal_vect,genss)
    return ideal(mixed_ring_xy, ideal_vect), new_rep_mat
end

#tested
function factorise(x::MPolyRingElem)
    #x has to be a monomial. No check implemented. 
    R = parent(x)
    V = exponent_vector(x,1)
    Factorisation = Vector{Tuple{elem_type(R), Int}}()
    for i in 1:length(V)
        push!(Factorisation, (gens(R)[i], V[i]))
    end
    return Factorisation
end


function proj_of_image_ideal(G::ReductiveGroup)
    W = image_ideal(G)
    mixed_ring_xy = base_ring(W[2])
    n = G.vector_space_dimension
    m = G.group[2]
    r = G.direct_product[2]
    #use parallelised groebner bases here. This is the bottleneck!
    return eliminate(W[1], gens(mixed_ring_xy)[(2*n)+1:(2*n)+(r*m^2)]), W[2]
end


#evaluate at y = 0 
function generators(G::ReductiveGroup, X::MPolyIdeal)
    n = G.vector_space_dimension
    m = G.group[2]
    gbasis = gens(X) 
    length(gbasis) == 0 && return gbasis,new_rep_mat
    mixed_ring_xy = parent(gbasis[1])
    
    #to evaluate gbasis at y = 0
    ev_gbasis = Vector{Union{AbstractAlgebra.Generic.MPoly{nf_elem}, MPolyRingElem}}(undef,0)
    for elem in gbasis
        b = mixed_ring_xy()
        mons = collect(monomials(elem))
        coeffs = collect(coefficients(elem))
        for i in 1:length(mons)
            if (exponent_vector(mons[i],1)[n+1:2*n] == [0 for i in 1:n])
                b += coeffs[i]*mons[i]
            end
        end
        b != 0 && push!(ev_gbasis, b)
    end
    
    #find representatives mod (det_ -1)
    #Did not make a difference. 
    #map_z_to_zz = hom(parent(det_), mixed_ring_xy, gens(mixed_ring_xy)[2*n+1:(2*n)+m^2])
    #new_det_ = map_z_to_zz(det_)
    #R,projection_ = quo(mixed_ring_xy, ideal(mixed_ring_xy, new_det_ -1))
    #V = Vector{MPolyRingElem}(undef,0)
    #for elem in ev_gbasis
    #    push!(V,lift(projection_(elem)))
    #end
    #@show V
    
    
    #grading starts here. 
    mixed_ring_graded, (x,y,zz) = grade(mixed_ring_xy)
    mapp = hom(mixed_ring_xy, mixed_ring_graded, gens(mixed_ring_graded))
    ev_gbasis_new = [mapp(ev_gbasis[i]) for i in 1:length(ev_gbasis)]
    if length(ev_gbasis_new) == 0
        return [mixed_ring_graded()], new_rep_mat
    end
    return minimal_generating_set(ideal(ev_gbasis_new))
end

#now we have to perform reynolds operation. This will happen in mixed_ring_xy. 
#the elements returned will be in the polynomial ring K[X]
function inv_generators(G::ReductiveGroup, ringg::MPolyRing, M::AbstractAlgebra.Generic.MatSpaceElem, X::MPolyIdeal)
    genss = generators(G, X)
    if length(genss) == 0
        return Vector{elem_type(ringg)}()
    end
    mixed_ring_xy = parent(genss[1])
    R = base_ring(M)
    m = G.group[2]
    n = G.vector_space_dimension
    r = G.direct_product[2]
    mapp = hom(R, mixed_ring_xy, gens(mixed_ring_xy))
    new_rep_mat = matrix(mixed_ring_xy,n,n,[mapp(M[i,j]) for i in 1:n, j in 1:n])
    #we need det_
    det_ = det(G.canonical_representation)
    mapp_ = hom(parent(det_), mixed_ring_xy, gens(mixed_ring_xy)[(2*n)+1:(2*n)+(r*m^2)])
    new_det = mapp_(det_)
    if G.group[1] == :SL && !G.direct_product[1] #TODO other types of reductive groups
        new_gens_wrong_ring = [G.reynolds_operator(genss[i], new_rep_mat, new_det, m) for i in 1:length(genss)]
    else 
        return nothing
    end
    img_genss = vcat(gens(ringg), zeros(ringg, n+r*m^2))
    mixed_to_ring = hom(mixed_ring_xy, ringg, img_genss)
    new_gens = Vector{MPolyDecRingElem}(undef,0)
    for elemm in new_gens_wrong_ring
        if elemm != 0
        push!(new_gens, mixed_to_ring(elemm))
        end
    end
    if length(new_gens) == 0
        return [ringg()]
    end
    
    #remove ugly coefficients: 
    V= Vector{QQFieldElem}[]
    for elem in new_gens
        V = vcat(V,collect(coefficients(elem)))
    end
    maxx = maximum([abs(denominator(V[i])) for i in 1:length(V)])
    minn = minimum([abs(numerator(V[i])) for i in 1:length(V)])
    new_gens_ = Vector{MPolyDecRingElem}(undef,0)
    for elem in new_gens
        if denominator((maxx*elem)//minn) != 1
            error("den not 1")
        end
        push!(new_gens_, numerator((maxx*elem)//minn))
    end
    return new_gens_
end


function mu_star(new_rep_mat::AbstractAlgebra.Generic.MatSpaceElem)
    mixed_ring_xy = parent(new_rep_mat[1,1])
    n = ncols(new_rep_mat)
    vars = matrix(mixed_ring_xy,n,1,[gens(mixed_ring_xy)[i] for i in 1:n])
    new_vars = new_rep_mat*vars 
    D = Dict([])
    for i in 1:n
        push!(D, gens(mixed_ring_xy)[i]=>new_vars[i])
    end
    return D
end

function reynolds__(elem::MPolyDecRingElem, new_rep_mat::AbstractAlgebra.Generic.MatSpaceElem, new_det::MPolyDecRingElem, m::Int)
    D = mu_star(new_rep_mat)
    mixed_ring_xy = parent(elem)
    sum = mixed_ring_xy()
    #mu_star: 
    for monomial in monomials(elem)
        k = mixed_ring_xy(1)
        factors = factorise(monomial)
        for j in 1:length(factors)
            if factors[j][2] != 0
                neww = getindex(D, factors[j][1])
                k = k*(neww^(factors[j][2]))
            end
        end
        sum += k
    end
    t = needed_degree(sum, m)
    if !divides(t, m)[1]
        return parent(elem)()
    else
        p = divexact(t, m)
    end
    num = omegap(p, new_det, sum)
    #num = omegap(p, new_det, elem)
    den = omegap(p, new_det, (new_det)^p)
    if !(denominator(num//den)==1)
        error("denominator of reynolds not rational")
    end
    return numerator(num//den)
end

function needed_degree(elem::MPolyDecRingElem, m::Int)
    R = parent(elem)
    n = numerator((ngens(R) - m^2)//2)
    extra_ring, zzz = PolynomialRing(base_ring(R), "zzz"=>1:m^2)
    mapp = hom(R,extra_ring, vcat([1 for i in 1:2*n], gens(extra_ring)))
    return total_degree(mapp(elem))
end

#works
function omegap(p::Int, det_::MPolyDecRingElem, f::MPolyDecRingElem)
    parent(det_) == parent(f) || error("Omega process ring error")
    action_ring = parent(det_)
    detp = (det_)^p
    monos = collect(monomials(detp))
    coeffs = collect(coefficients(detp))
    h = action_ring()
    for i in 1:length(monos)
        exp_vect = exponent_vector(monos[i], 1)
        x = f
        for i in 1:length(exp_vect)
            for j in 1:exp_vect[i]
                x = derivative(x, gens(action_ring)[i])
                if x == 0
                    break
                end
            end
        end
        h += coeffs[i]*x
    end
    return h
end

function reduce_gens_(v::Vector{MPolyDecRingElem})
    new_gens_ = [v[1]]
    for i in 1:length(v)-1
        if !(v[i+1] in ideal(v[1:i]))
            push!(new_gens_, v[i+1])
        end
    end
    return new_gens_
end




