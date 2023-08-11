mutable struct InvariantRing
    field::Union{NumField, QQField}
    poly_ring::MPolyRing #CHECK
    group::Tuple{Symbol,Int64}
    group_equations::Union{Vector{QQMPolyRingElem},Vector{AbstractAlgebra.Generic.MPoly{nf_elem}}} 
    group_rep:: T where T <: AbstractAlgebra.Generic.MatSpaceElem
    generators::Vector{MPolyRingElem}

    function InvariantRing(sym::Symbol, rep_mat::T where T <:AbstractAlgebra.Generic.MatSpaceElem)
        #sym != SL && return nothing
        z = new()
        R = parent(rep_mat[1,1]) 
        z.field = base_ring(R)
        (typeof(z.field) <: NumField || typeof(z.field) == QQField) || @error("Field should be rational or number field")
        m = Int(sqrt(ngens(R)))
        n = ncols(rep_mat)
        z.poly_ring, _ = PolynomialRing(z.field, "X" => 1:n)
        z.group = (sym, m)
        M = matrix(R,1,m,gens(R)[1:m])
        for i in 1:m-1 
            M = vcat(M, matrix(R,1,m,gens(R)[(i)*m+1:(i+1)*m]))
        end
        det_ = det(M)
        z.group_equations = [det_ - 1]
        z.group_rep = rep_mat
        z.generators = inv_generators(rep_mat, z.poly_ring, det_)
        return z
    end

    #INCOMPLETE. TODO
    function InvariantRing(sym::Symbol, m::Int64, field::Field, sym_deg::Int64)
        #sym != SL && return nothing
        #sym_rep_mat = TODO
        z = new()
        z.field = field
        #need to check if this is a NumField
        z.group_rep = rep_mat_(m, field, sym_deg)
        n = ncols(z.group_rep)
        z.poly_ring, __ = PolynomialRing(z.field, "X" => 1:n)
        z.group = (sym, m)
        group_poly_ring, Z = PolynomialRing(field, "Z"=>(1:m,1:m))
        M = matrix(group_poly_ring,m,m,[Z[i,j] for i in 1:m, j in 1:m])
        det_ = det(M)
        z.group_equations = [det_ - 1]
        z.generators = inv_generators(z.group_rep, z.poly_ring, det_)
        return z
    end
end

#tested
function image_ideal(rep_mat::T where T <:AbstractAlgebra.Generic.MatSpaceElem)
    R = parent(rep_mat[1,1])
    n = ncols(rep_mat)
    m = Int(sqrt(ngens(R)))
    K = base_ring(R)
    mixed_ring_xy, x, y, zz = PolynomialRing(K, "x"=>1:n, "y"=>1:n, "zz"=>(1:m,1:m))
    # naming the variables zz here and z in the group ring
    #for determinant - 
    M1 = matrix(mixed_ring_xy,m,m,[zz[i,j] for i in 1:m for j in 1:m])
    ztozz = hom(R,mixed_ring_xy, gens(mixed_ring_xy)[(2*n)+1:(2*n)+(m^2)])
    #rep_mat in the new ring
    new_rep_map = matrix(mixed_ring_xy,n,n,[image_in_map(rep_mat[i,j], ztozz) for i in 1:n, j in 1:n])
    new_vars = new_rep_map*[x[i] for i in 1:n]
    ideal_vect = [y[i] - new_vars[i] for i in 1:n]
    Base.push!(ideal_vect,det(M1) - 1)
    return (ideal(mixed_ring_xy, ideal_vect), new_rep_map)
end

#tested
function image_in_map(X::Union{MPolyRingElem,AbstractAlgebra.Generic.MPoly}, f::Oscar.MPolyAnyMap)
    #x in f.domain || error
    answer = f.codomain()
    mons = collect(monomials(X))
    coeffs = collect(coefficients(X))
    for i in 1:length(mons)
        Factorisation = factorise(mons[i])
        y = f.codomain(1)
        for i in 1:length(Factorisation)
            if Factorisation[i][2] != 0
                for j in 1:Factorisation[i][2]
                    y = y*f.img_gens[i]
                end
            end
        end
        answer += coeffs[i]y
    end
    return answer
end

#tested
function factorise(x::MPolyRingElem)
    #x has to be a monomial. No check implemented. 
    R = parent(x)
    V = exponent_vector(x,1)
    Factorisation = Vector{Tuple{MPolyRingElem,Int64}}(undef,0)
    for i in 1:length(V)
        Base.push!(Factorisation, (gens(R)[i], V[i]))
    end
    return Factorisation
end


function proj_of_image_ideal(rep_mat::T where T <:AbstractAlgebra.Generic.MatSpaceElem)
    W = image_ideal(rep_mat)
    mixed_ring_xy = base_ring(W[1])
    n = ncols(rep_mat)
    m = Int(sqrt(ngens(mixed_ring_xy) - 2*n))
    #use parallelised groebner bases here. This is the bottleneck!
    return (groebner_basis(eliminate(W[1], gens(mixed_ring_xy)[(2*n)+1:(2*n)+(m^2)])), W[2])
    #return (eliminate(W[1], gens(mixed_ring_xy)[(2*n)+1:(2*n)+(m^2)]))
end


#evaluate at y = 0 
function generators(rep_mat::T where T <:AbstractAlgebra.Generic.MatSpaceElem)
    (G, new_rep_mat) = proj_of_image_ideal(rep_mat)
    gbasis = collect(G) #is this needed?
    length(gbasis) == 0 && return gbasis,new_rep_mat
    mixed_ring_xy = parent(gbasis[1])
    #to evaluate gbasis at y = 0
    ev_gbasis = Vector{Union{AbstractAlgebra.Generic.MPoly{nf_elem}, MPolyRingElem}}(undef,0)
    n = ncols(rep_mat)
    for elem in gbasis
        b = mixed_ring_xy()
        for monomial in monomials(elem)
            for j in 1:n 
                if divides(monomial, gens(mixed_ring_xy)[n+j])[1]
                    @goto lab
                end
            end
            b += monomial 
            @label lab
        end
        b != 0 && Base.push!(ev_gbasis, b)
    end
    return ev_gbasis, new_rep_mat
end

#now we have to perform reynolds operation. This will happen in mixed_ring_xy. 
#the elements returned will be in the polynomial ring K[X]
function inv_generators(rep_mat::T where T <:AbstractAlgebra.Generic.MatSpaceElem, ringg::MPolyRing, det_::MPolyRingElem)
    genss, new_rep_mat = generators(rep_mat)
    if length(genss) == 0
        return Vector{MPolyRingElem}(undef,0)
    end
    mixed_ring_xy = parent(genss[1])
    m = Int(sqrt(ngens(parent(det_))))
    n = ncols(rep_mat)
    mapp_ = hom(parent(rep_mat[1,1]), mixed_ring_xy, gens(mixed_ring_xy)[(2*n)+1:(2*n)+(m^2)])
    new_det = image_in_map(det_, mapp_)
    #new_gens_wrong_ring = Vector{Union{MPolyRingElem,AbstractAlgebra.Generic.MPoly}}(undef,0) #what do i do here?
    #for gen in genss
     #   Base.push!(new_gens_wrong_ring, reynolds__(gen, m, new_rep_mat, new_det))
    #end
    new_gens_wrong_ring = [reynolds__(genss[i], new_rep_mat, new_det) for i in 1:length(genss)]
    img_genss = vcat(gens(ringg), zeros(ringg, n+m^2))
    mixed_to_ring = hom(mixed_ring_xy, ringg, img_genss)
    new_gens = Vector{MPolyRingElem}(undef,0)
    for elemm in new_gens_wrong_ring
        Base.push!(new_gens, image_in_map(elemm, mixed_to_ring))
    end
    return new_gens
end

function image_in_action_ring(X::MPolyRingElem, map::MPolyAnyMap)
    #X in f.domain || error
    n = ngens(f.domain)
    answer = parent(X)()
    mons = collect(monomials(X))
    coeffs = collect(coefficients(X))
    for i in 1:length(mons)
        Factorisation = factorise(mons[i])
        y = f.codomain(1)
        for i in 1:length(Factorisation)
            if Factorisation[i][2] != 0
                for j in 1:Factorisation[i][2]
                    y = y*f.img_gens[i+n]
                end
            end
        end
        answer += coeffs[i]*y
    end
    return answer
end

function mu_star(new_rep_mat::T where T <:AbstractAlgebra.Generic.MatSpaceElem)
    mixed_ring_xy = parent(new_rep_mat[1,1])
    n = ncols(new_rep_mat)
    vars = matrix(mixed_ring_xy,n,1,[gens(mixed_ring_xy)[i] for i in 1:n])
    new_vars = new_rep_mat*vars 
    D = Dict([])
    for i in 1:n
        Base.push!(D, gens(mixed_ring_xy)[i]=>new_vars[i])
    end
    return D
end

function reynolds__(elem::MPolyRingElem, new_rep_mat::T where T<:AbstractAlgebra.Generic.MatSpaceElem, new_det::MPolyRingElem)
    n = ncols(new_rep_mat)
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
    if !divides(total_degree(sum), n)[1]
        return parent(elem)()
    else
        p = divexact(total_degree(sum), n)
    end
    num = omegap(p, new_det, sum)
    den = omegap(p, new_det, (new_det)^p)
    if denominator(num//den) != 1
        @error("denominatior of reynolds not rational")
    end
    return numerator(num//den)
end

#works
function omegap(p::Int64, det_::MPolyRingElem, f::MPolyRingElem)
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

###############

#for the second constructor function of InvariantRing

function rep_mat_(m::Int64, field::Field, sym_deg::Int64)
    n = binomial(m + sym_deg - 1, m - 1)
    mixed_ring, t, z, a = PolynomialRing(field, "t"=> 1:m, "z"=> (1:m, 1:m), "a" => 1:n)
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
    group_ring, Z = PolynomialRing(field, "Z"=>(1:m, 1:m))
    mapp = hom(mixed_ring, group_ring, vcat([0 for i in 1:m], gens(group_ring), [0 for i in 1:n]))
    Mat = matrix(group_ring, n, n, [image_in_map(mat[i,j], mapp) for i in 1:n, j in 1:n])
    return Mat            
end

function degree_basis(R::MPolyRing,m::Int64, t::Int64)
    v = R(1)
    genss = gens(R)
    n = length(gens(R)[1:m])
    C = zero_matrix(Int64,n,n)
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



function Base.show(io::IO, R::InvariantRing)
    print(io, "Invariant Rinf of", "\n")
    print(io, R.poly_ring, "\n")
    print(io, "under group action of ", R.group[1], R.group[2], "\n")
    print(io, "represented as ", "\n")
    print(io, R.group_rep, "\n")
    print(io, "Generated by ", "\n")
    print(io, R.generators)
end