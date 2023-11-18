import Oscar.gens, AbstractAlgebra.direct_sum, Oscar.invariant_ring
#export ReductiveGroup, reductive_group, representation_matrix, group, mu_star, reynolds_operator, group_ideal, canonical_representation, natural_representation, vector_space_dimension, 
#export InvariantRing, invariant_ring, gens, hilbert_ideal, derksen_ideal
##########################
#Reductive Groups
##########################
mutable struct ReductiveGroup
    field::Field #characteristic zero. implement check? 
    group::Tuple{Symbol, Int}
    group_ideal::MPolyIdeal
    reynolds_operator::Function
    canonical_representation::AbstractAlgebra.Generic.MatSpaceElem
    
    # function ReductiveGroup(sym::Symbol, m::Int, fld::Field) 
    #     G = new()
    #     if sym != :SL
    #         error("Only implemented for SLm")
    #     end
    #     G.field = fld
    #     G.group = (sym,m)
    #     G.reynolds_operator = reynolds_slm
    #     R, M0 = polynomial_ring(fld, :z => (1:m,1:m))
    #     M = matrix(M0)
    #     G.canonical_representation = M
    #     G.group_ideal = ideal([det(M) - 1])
    #     #base ring of M has to be the same as the representation matrix when that is created later.
    #     return G
    # end

    function ReductiveGroup(sym::Symbol, m::Int, fld::Field) #have not decided the representation yet
        R, _ = polynomial_ring(fld, :z => (1:m,1:m))
        return ReductiveGroup(sym, m, R)
    end
    
    function ReductiveGroup(sym::Symbol, m::Int, ring::MPolyRing) #the ring input is the group ring
        G = new()
        if sym != :SL
            error("Only implemented for SLm")
        end
        G.field = base_ring(ring)
        @assert m^2  == ngens(ring)
        G.group = (sym,m)
        G.reynolds_operator = reynolds__
        M = matrix(ring, m, m, gens(R))
        G.canonical_representation = M
        G.group_ideal = ideal([det(M) - 1])
        #base ring of M has to be the same as the representation matrix when that is created later.
        return G
    end
end

function Base.show(io::IO, G::ReductiveGroup)
    print(io, "Reductive group ", G.group[1], G.group[2])
end

function reductive_group(sym::Symbol, m::Int, F::Field)
    return ReductiveGroup(sym,m,F)
end

function reductive_group(sym::Symbol, m::Int, R::MPolyRing)
    return ReductiveGroup(sym,m,R)
end

group(G::ReductiveGroup) = G.group
reynolds_operator(G::ReductiveGroup) = G.reynolds_operator
group_ideal(G::ReductiveGroup) = G.group_ideal
canonical_representation(G::ReductiveGroup) = G.canonical_representation
natural_representation(G::ReductiveGroup) = G.canonical_representation

#####################
#Representation
#####################

mutable struct RepresentationReductiveGroup
    group::ReductiveGroup
    rep_mat::AbstractAlgebra.Generic.MatSpaceElem
    sym_deg::Tuple{Bool, Int}
    reynolds_v::Function
    
    #representation of group G over symmetric degree d
    function RepresentationReductiveGroup(G::ReductiveGroup, d::Int)
        R = new()
        R.group = G
        R.rep_mat = rep_mat_(G, d)
        R.sym_deg = (true, d)
        R.reynolds_v = reynolds_v_slm
        return R
    end
    
    #does not check the given matrix.
    function RepresentationReductiveGroup(G::ReductiveGroup, M::AbstractAlgebra.Generic.MatSpaceElem)
        @req base_ring(M) == base_ring(G.group_ideal) "Group ideal and representation matrix must have same parent ring"
        R = new()
        R.group = G
        R.rep_mat = M
        R.sym_deg = (false, 0)
        R.reynolds_v = reynolds_v_slm
        return R
    end
end

function representation_reductive_group(G::ReductiveGroup, d::Int)
    return RepresentationReductiveGroup(G, d)
end

representation_matrix(R::RepresentationReductiveGroup) = R.rep_mat
reductive_group(R::RepresentationReductiveGroup) = R.group
function vector_space_dimension(R::RepresentationReductiveGroup)
    return ncols(R.rep_mat)
end

function Base.show(io::IO, R::RepresentationReductiveGroup)
    println(io, "Representation of ", R.group.group[1], R.group.group[2])
    if R.sym_deg[1]
        print(io, "over symmetric forms of degree ", R.sym_deg[2])
    else
        println(io, "with representation matrix")
        show(io, R.rep_mat)
    end
end

###############

#computes the representation matrices of SL_m acting over m-forms of symmetric degree sym_deg
function rep_mat_(G::ReductiveGroup, sym_deg::Int)
    G.group[1] == :SL || error("Only implemented for SLm")
    m = G.group[2]
    n = binomial(m + sym_deg - 1, m - 1)
    mixed_ring, t, z = PolynomialRing(G.field, "t"=> 1:m, "z"=> (1:m, 1:m))
    group_mat = matrix(mixed_ring, z)
    new_vars = group_mat*t

    b = reverse(degree_basis(mixed_ring,m, sym_deg))
    
    # transform the b elements
    new_vars_plus_rest = vcat(new_vars, gens(mixed_ring)[m+1:nvars(mixed_ring)])
    images_of_b = [evaluate(f, new_vars_plus_rest) for f in b]

    mat = zero_matrix(mixed_ring, n, n)
    for i in 1:n
        f = images_of_b[i]
        # express f as a linear combination of degree_basis
        for t in terms(f)
            for j in 1:n
                if divides(t, b[j])[1]
                    c = leading_coefficient(b[j])
                    mat[i,j] += t / b[j] * c
                    #mat[i,j] += t / b[j]  # FIXME: this is the natural thing but gives "wrong" (?) results
                end
            end
        end
    end
    #we have to return mat in a different ring! 
    group_ring = base_ring(G.group_ideal)
    mapp = hom(mixed_ring, group_ring, vcat([0 for i in 1:m], gens(group_ring)))
    Mat = matrix(group_ring, n, n, [mapp(mat[i,j]) for i in 1:n, j in 1:n])
    return Mat            
end

function degree_basis(R::MPolyRing, m::Int, t::Int)
    C = zero_matrix(Int, m, m)
    for i in 1:m
      C[i,i] = -1 #find a better way to write this TODO
    end
    d = zeros(n)
    A = ones(n)
    b = [t]
    P = Polyhedron((C, d), (A, b))
    L = lattice_points(P)
    W = [ R([multinomial(t,l)], [l]) for l in L ]
    return W
end

function multinomial(n::Int, v::Union{Vector{Int64},PointVector{ZZRingElem}})
    l = length(v)
    x = 1
    for i in 1:l
        x = x*factorial(v[i])
    end
    return Int(factorial(n)/x)
end

#############################
#DIRECT SUM

# #TODO
# function rep_mat_(m::Int, v::Vector{Int}, dir_sum::Bool)
#     if dir_sum
#         n = num_of_as(m,v)
#         mixed_ring, t, z, a = PolynomialRing(QQ, "t"=> 1:m, "z"=> (1:m, 1:m), "a" => 1:n)
#         group_mat = matrix(mixed_ring, m,m,[z[i,j] for i in 1:m, j in 1:m])
#         vars = [t[i] for i in 1:m]
#         new_vars = vars*group_mat
#         big_matrix = matrix(mixed_ring, 0,0,[])
#         for sym_deg in v
#             b = degree_basis(mixed_ring,m, sym_deg)
#             sum_ = mixed_ring()
#             for j in 1:length(b)
#                 prod_ = mixed_ring(1)
#                 factors_ = factorise(b[j])
#                 for i in 1:length((factors_)[1:m])
#                     prod_ = prod_*new_vars[i]^(factors_[i][2])
#                 end
#                 prod_ = prod_*a[j]
#                 sum_ += prod_
#             end
#             mons = collect(monomials(sum_))
#             coeffs = collect(coefficients(sum_))
#             N = binomial(m + sym_deg - 1, m - 1)
#             mat = matrix(mixed_ring, N,N,[0 for i in 1:N^2])
#             for i in 1:N
#                 for j in 1:N
#                     for k in 1:length(mons)
#                         if divides(mons[k], a[i])[1] && divides(mons[k], b[j])[1]
#                             mat[i,j] += coeffs[k]*numerator((mons[k]//b[j])//a[i])
#                         end
#                     end
#                 end
#             end
#             big_matrix = direct_sum(big_matrix,mat)
#         end
#         #we have to return mat in a different ring! 
#         group_ring, Z = PolynomialRing(QQ, "Z"=>(1:m, 1:m))
#         mapp = hom(mixed_ring, group_ring, vcat([0 for i in 1:m], gens(group_ring), [0 for i in 1:n]))
#         Mat = matrix(group_ring, nrows(big_matrix), nrows(big_matrix), [mapp(big_matrix[i,j]) for i in 1:nrows(big_matrix), j in 1:nrows(big_matrix)])
#         return Mat
#     else
#         l = length(v)
#         #there are l*m^2 z's
#         ringg,Z = PolynomialRing(QQ,"Z"=>1:l*m^2)
#         prod_mat = matrix(ringg,1,1,[1])
#         for i in 1:l
#             mat = rep_mat_(m,v[i])
#             R = base_ring(mat)
#             mapp = hom(R,ringg, gens(ringg)[((i-1)*m^2)+1:i*m^2])
#             mat_ = matrix(ringg, nrows(mat), ncols(mat), [0 for i in 1:nrows(mat)*ncols(mat)])
#             for i in 1:nrows(mat)
#                 for j in 1:ncols(mat)
#                     mat_[i,j] = mapp(mat[i,j])
#                 end
#             end
#             prod_mat = kronecker_product(prod_mat, mat_)
#         end
#         return prod_mat
#     end
# end

function direct_sum(M::AbstractAlgebra.Generic.MatSpaceElem, N::AbstractAlgebra.Generic.MatSpaceElem)
    #parent(M[1,1]) == parent(N[1,1]) || @error("not same ring")
    mr = nrows(M)
    mc = ncols(M)
    nr = nrows(N)
    nc = ncols(N)
    mat_ = zero_matrix(base_ring(M), mr+nr, mc+nc)
    mat_[1:mr, 1:mc] = M
    mat_[mr+1:mr+nr, mc+1:mc+nc] = N
    return mat_
end

function num_of_as(m::Int, v::Vector{Int})
    max = maximum(v)
    return binomial(m + max - 1, m - 1)
end

##########################
#Invariant Rings of Reductive groups
##########################
mutable struct InvariantRing
    field::Field
    poly_ring::MPolyDecRing #graded
    
    group::ReductiveGroup
    representation::RepresentationReductiveGroup
    
    fundamental::Vector{MPolyDecRingElem}
    
    reynolds_operator::Function
    
    NullConeIdeal::MPolyIdeal
    
    #This should be the only one that is used
    function InvariantRing(R::RepresentationReductiveGroup) #here G already contains information n and rep_mat
        z = new()
        n = ncols(R.rep_mat)
        z.representation = R
        z.group = R.group
        G = z.group
        z.field = G.field
        z.poly_ring, __ = graded_polynomial_ring(G.field, "X" => 1:n)
        z.reynolds_operator = reynolds_v_slm
        I, M = proj_of_image_ideal(G, R.rep_mat)
        z.NullConeIdeal = ideal(generators(G, I, R.rep_mat))
        z.fundamental = inv_generators(z.NullConeIdeal, G, z.poly_ring, M, I, z.reynolds_operator)
        return z
    end
    
    function InvariantRing(R::RepresentationReductiveGroup, ring::MPolyRing)
        n = ncols(R.rep_mat)
        n == ngens(ring) || error("The given polynomial ring is not compatible.")
        z = new()
        z.representation = R
        z.group = R.group
        G = R.group
        z.field = G.field
        if is_graded(ring)
            z.poly_ring = ring[1]
        else
            z.poly_ring = grade(ring)
        end
        z.reynolds_operator = reynolds_v_slm
        I, M = proj_of_image_ideal(G, R.rep_mat)
        z.NullConeIdeal = ideal(generators(G, I, R.rep_mat))
        z.fundamental = inv_generators(z.NullConeIdeal, G, z.poly_ring, M, I, z.reynolds_operator)
        return z
    end
end

invariant_ring(R::RepresentationReductiveGroup) = InvariantRing(R)
fundamental_invariants(R::InvariantRing) = R.fundamental
null_cone_ideal(R::InvariantRing) = R.NullConeIdeal

function Base.show(io::IO, R::InvariantRing) #TODO compact printing
    #if get(io, :supercompact, false)
        print(io, "Invariant Ring of", "\n")
        show(io, R.poly_ring)
        print(io, " under group action of ", R.group.group[1], R.group.group[2], "\n", "\n")
        print(io, "Generated by ", "\n")
    join(io, R.fundamental, "\n")
    #else
     #   print(io, "Invariant Ring under ")
      #  show(io, R.group[1], R.group[2])
    #end
end

#I think this should work for everyone except tensor products.
function image_ideal(G::ReductiveGroup, rep_mat::AbstractAlgebra.Generic.MatSpaceElem)
    R = base_ring(rep_mat)
    n = ncols(rep_mat)
    m = G.group[2]
    #r = G.direct_product[2]
    mixed_ring_xy, x, y, zz = PolynomialRing(G.field, "x"=>1:n, "y"=>1:n, "zz"=>1:m^2)
    # naming the variables zz here and z in the group ring
    #for determinant - 
    #M1 = matrix(mixed_ring_xy,m,m,[zz[i,j] for i in 1:m for j in 1:m])
    ztozz = hom(R,mixed_ring_xy, gens(mixed_ring_xy)[(2*n)+1:(2*n)+(m^2)])
    genss = [ztozz(f) for f in gens(G.group_ideal)]
    #rep_mat in the new ring
    new_rep_mat = matrix(mixed_ring_xy,n,n,[ztozz(rep_mat[i,j]) for i in 1:n, j in 1:n])
    new_vars = new_rep_mat*x
    ideal_vect = y - new_vars 
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


function proj_of_image_ideal(G::ReductiveGroup, rep_mat::AbstractAlgebra.Generic.MatSpaceElem)
    W = image_ideal(G, rep_mat)
    mixed_ring_xy = base_ring(W[2])
    n = ncols(rep_mat)
    m = G.group[2]
    #r = G.direct_product[2]
    #use parallelised groebner bases here. This is the bottleneck!
    return eliminate(W[1], gens(mixed_ring_xy)[(2*n)+1:(2*n)+(m^2)]), W[2]
end


#evaluate at y = 0 
function generators(G::ReductiveGroup, X::MPolyIdeal, rep_mat::AbstractAlgebra.Generic.MatSpaceElem)
    n = ncols(rep_mat)
    m = G.group[2]
    gbasis = gens(X) 
    length(gbasis) == 0 && return gbasis,new_rep_mat
    mixed_ring_xy = parent(gbasis[1])
    
    #to evaluate gbasis at y = 0
    ev_gbasis = elem_type(mixed_ring_xy)[]
    for elem in gbasis
        b = mixed_ring_xy()
        for t in terms(elem)
            if exponent_vector(t,1)[n+1:2*n] == [0 for i in 1:n]
                b += t
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
        return [mixed_ring_graded()]
    end
    return minimal_generating_set(ideal(ev_gbasis_new))
end

#now we have to perform reynolds operation. This will happen in mixed_ring_xy. 
#the elements returned will be in the polynomial ring K[X]
function inv_generators(I::MPolyIdeal, G::ReductiveGroup, ringg::MPolyRing, M::AbstractAlgebra.Generic.MatSpaceElem, X::MPolyIdeal, reynolds_function::Function)
    genss = gens(I)
    if length(genss) == 0
        return Vector{elem_type(ringg)}()
    end
    mixed_ring_xy = parent(genss[1])
    R = base_ring(M)
    m = G.group[2]
    n = ncols(M)
    mapp = hom(R, mixed_ring_xy, gens(mixed_ring_xy))
    new_rep_mat = matrix(mixed_ring_xy,n,n,[mapp(M[i,j]) for i in 1:n, j in 1:n])
    #we need det_
    det_ = det(G.canonical_representation)
    mapp_ = hom(parent(det_), mixed_ring_xy, gens(mixed_ring_xy)[(2*n)+1:(2*n)+(m^2)])
    new_det = mapp_(det_)
    if G.group[1] == :SL #TODO other types of reductive groups
        new_gens_wrong_ring = [reynolds_function(genss[i], new_rep_mat, new_det, m) for i in 1:length(genss)]
    else 
        return nothing
    end
    img_genss = vcat(gens(ringg), zeros(ringg, n+m^2))
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
    V= Vector{FieldElem}[]
    new_gens_ = Vector{MPolyDecRingElem}(undef,0)
    for elem in new_gens
        V = vcat(V,collect(coefficients(elem)))
        maxx = maximum([abs(denominator(V[i])) for i in 1:length(V)])
        minn = minimum([abs(numerator(V[i])) for i in 1:length(V)])
        if denominator((maxx*elem)//minn) != 1
            error("den not 1")
        end
        push!(new_gens_, numerator((maxx*elem)//minn))
        V= Vector{FieldElem}[]
    end
    return new_gens_
end


function mu_star(new_rep_mat::AbstractAlgebra.Generic.MatSpaceElem)
    mixed_ring_xy = parent(new_rep_mat[1,1])
    n = ncols(new_rep_mat)
    vars = matrix(mixed_ring_xy,1,n,[gens(mixed_ring_xy)[i] for i in 1:n])
    new_vars = vars*transpose(new_rep_mat)
    D = Dict([])
    for i in 1:n
        push!(D, gens(mixed_ring_xy)[i]=>new_vars[i])
    end
    return D
end

function reynolds_v_slm(elem::MPolyDecRingElem, new_rep_mat::AbstractAlgebra.Generic.MatSpaceElem, new_det::MPolyDecRingElem, m::Int)
    D = mu_star(new_rep_mat)
    mixed_ring_xy = parent(elem)
    sum_ = mixed_ring_xy()
    #mu_star: 
    mons = collect(monomials(elem))
    coeffs = collect(coefficients(elem))
    for i in 1:length(mons)
        k = mixed_ring_xy(1)
        factors = factorise(mons[i])
        for j in 1:length(factors)
            if factors[j][2] != 0
                neww = getindex(D, factors[j][1])
                k = k*(neww^(factors[j][2]))
            end
        end
        sum_ += coeffs[i]*k
    end
    t = needed_degree(sum_, m)
    if !divides(t, m)[1]
        return parent(elem)()
    else
        p = divexact(t, m)
    end
    #num = omegap_(p, new_det, sum_)
    #den = omegap_(p, new_det, (new_det)^p)
    #if !(denominator(num//den)==1)
    #    error("denominator of reynolds not rational")
    #end
    #return numerator(num//den)
    return reynolds_slm(sum_, new_det, p)
end

function reynolds_slm(elem::MPolyElem, det_::MPolyElem, p::Int)
    num = omegap_(p,det_, elem)
    den = omegap_(p,det_,det_^p)
    if !(denominator(num//den)==1)
        error("denominator of reynolds not rational")
    end
    return numerator(num//den)
end

function needed_degree(elem::MPolyDecRingElem, m::Int)
    R = parent(elem)
    n = numerator((ngens(R) - m^2)//2)
    extra_ring, _= PolynomialRing(base_ring(R), "z"=>1:m^2)
    mapp = hom(R,extra_ring, vcat([1 for i in 1:2*n], gens(extra_ring)))
    return total_degree(mapp(elem))
end

#works #Unused.
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
                @show derivative(x, gen(action_ring,i)) == derivative(x,i)
                x = derivative(x, gen(action_ring,i))
                if x == 0
                    break
                end
            end
        end
        h += coeffs[i]*x
    end
    return h
end

#alternate function for omegap. Works the same.
function omegap_(p::Int, det_::MPolyDecRingElem, f::MPolyDecRingElem)
    parent(det_) == parent(f) || error("Omega process ring error")
    action_ring = parent(det_)
    monos = collect(monomials(det_))
    coeffs = collect(coefficients(det_))
    for i in 1:p
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
        f = h
    end
    return f
end

#####################callable reynold's operator

function reynolds_operator(elem::MPolyElem, X::RepresentationReductiveGroup)
    X.group.group[1] == :SL || error("Only implemented for SLm")
    vector_ring = parent(elem)
    G = X.group
    n = Int(ngens(vector_ring))
    n == ncols(X.rep_mat) || error("group not compatible with element")
    m = G.group[2]
    R, _ = grade(PolynomialRing(G.field,"x"=>1:n, "y"=>1:n, "z"=>(1:m, 1:m))[1])
    map1 = hom(vector_ring, R, gens(R)[1:n])
    new_elem = map1(elem)
    group_ring = base_ring(X.rep_mat)
    map = hom(group_ring, R, gens(R)[2*n+1:2*n+m^2])
    new_rep_mat = matrix(R,n,n,[map(X.rep_mat[i,j]) for i in 1:n, j in 1:n])
    new_det = map(det(G.canonical_representation))
    f = X.reynolds_v(new_elem, new_rep_mat, new_det, m)
    reverse_map = hom(R, vector_ring, vcat(gens(vector_ring), [0 for i in 1:n+m^2]))
    return reverse_map(f)
end

#TODO
function mu_star(elem_::MPolyElem, R::RepresentationReductiveGroup)
    G = R.group
    n = ncols(R.rep_mat)
    ngens(parent(elem_)) == n || error("group not compatible with element")
    m = G.group[2]
    mixed_ring_xz, x, z = PolynomialRing(G.field, "x"=>1:n, "z"=>(1:m,1:m))
    group_ring = base_ring(R.rep_mat)
    map1 = hom(group_ring, mixed_ring_xz, gens(mixed_ring_xz)[n+1:n+m^2])
    new_rep_mat = matrix(mixed_ring_xz,n,n,[map1(R.rep_mat[i,j]) for i in 1:n, j in 1:n])
    D = mu_star(new_rep_mat)
    vector_ring = parent(elem_)
    map2 = hom(vector_ring, mixed_ring_xz, gens(mixed_ring_xz)[1:n])
    elem = map2(elem_)
    sum_ = mixed_ring_xz()
    #mu_star: 
    mons = collect(monomials(elem))
    coeffs = collect(coefficients(elem))
    for i in 1:length(mons)
        k = mixed_ring_xz(1)
        factors = factorise(mons[i])
        for j in 1:length(factors)
            if factors[j][2] != 0
                neww = getindex(D, factors[j][1])
                k = k*(neww^(factors[j][2]))
            end
        end
        sum_ += coeffs[i]*k
    end
    return sum_
end
