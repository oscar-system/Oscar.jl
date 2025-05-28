################################################################################
#
#  Tropical varieties
#  ==================
#  concrete subtype of TropicalVarietySupertype in variety_supertype.jl
#
################################################################################

@attributes mutable struct TropicalVariety{minOrMax,isEmbedded} <: TropicalVarietySupertype{minOrMax,isEmbedded}
    polyhedralComplex::PolyhedralComplex
    multiplicities::Vector{ZZRingElem}

    # tropical varieties need to be embedded
    function TropicalVariety{minOrMax,true}(Sigma::PolyhedralComplex,multiplicities::Vector{ZZRingElem}) where {minOrMax<:Union{typeof(min),typeof(max)}}
        return new{minOrMax,true}(Sigma,multiplicities)
    end
end



################################################################################
#
#  Printing
#
################################################################################

function Base.show(io::IO, tv::TropicalVariety{typeof(min), true})
    print(io, "Min tropical variety")
end
function Base.show(io::IO, tv::TropicalVariety{typeof(max), true})
    print(io, "Max tropical variety")
end



################################################################################
#
#  Constructors
#
################################################################################

@doc raw"""
    tropical_variety(Sigma::PolyhedralComplex, mult, minOrMax::Union{typeof(min),typeof(max)}=min)

Return the `TropicalVariety` whose polyhedral complex is `Sigma` with multiplicities `mult` and convention `minOrMax`. Here, `mult` is optional can be specified as a `Vector{ZZRingElem}` which represents a list of multiplicities on the maximal polyhedra in the order of `maximal_polyhedra(Sigma)`.  If `mult` is unspecified, then all multiplicities are set to one.

# Examples
```jldoctest
julia> Sigma = polyhedral_complex(incidence_matrix([[1],[2]]), [[0],[1]])
Polyhedral complex in ambient dimension 1

julia> tropical_variety(Sigma)
Min tropical variety

julia> mult = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
2-element Vector{ZZRingElem}:
 1
 1

julia> tropical_variety(Sigma,mult,min)
Min tropical variety

julia> mult = ZZ.([1,2])
2-element Vector{ZZRingElem}:
 1
 2

julia> tropical_variety(Sigma,mult,max)
Max tropical variety

```
"""
function tropical_variety(Sigma::PolyhedralComplex, mult::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    return TropicalVariety{typeof(minOrMax),true}(Sigma,mult)
end
function tropical_variety(Sigma::PolyhedralComplex, minOrMax::Union{typeof(min),typeof(max)}=min)
    mult = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
    return tropical_variety(Sigma,mult,minOrMax)
end

function tropical_variety(TropH::TropicalHypersurface{minOrMax,true}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    TropV = tropical_variety(polyhedral_complex(TropH),multiplicities(TropH), convention(TropH))

    if has_attribute(TropH,:algebraic_polynomial)
        set_attribute!(TropV,:algebraic_ideal,ideal([get_attribute(TropH,:algebraic_polynomial)]))
    end
    if has_attribute(TropH,:tropical_semiring_map)
        set_attribute!(TropV,:tropical_semiring_map,get_attribute(TropH,:tropical_semiring_map))
    end

    return TropV
end


function tropical_variety(TropL::TropicalLinearSpace{minOrMax,true}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    TropV = tropical_variety(polyhedral_complex(TropL),multiplicities(TropL), convention(TropL))

    if has_attribute(TropL,:algebraic_ideal)
        set_attribute!(TropV,:algebraic_ideal,get_attribute(TropL,:algebraic_ideal))
    end
    if has_attribute(TropL,:tropical_semiring_map)
        set_attribute!(TropV,:tropical_semiring_map,get_attribute(TropL,:tropical_semiring_map))
    end

    return TropV
end


function tropical_variety(Sigma::Vector{<:Polyhedron}, multiplicities::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    return tropical_variety(polyhedral_complex(Sigma; non_redundant=true),multiplicities,minOrMax)
end




################################################################################
#
#  Properties
#  ----------
#  none, see variety_supertype.jl for common properties of all tropical variety types
#
################################################################################



################################################################################
#
#  Tropical varieties of polynomial ideals
#  ----------
#  References for computing tropical varieties via Groebner complex traversal:
#    T. Bogart, A. Jensen, D. Speyer, B. Sturmfels, R. Thomas: Computing tropical varieties
#    T. Markwig, Y. Ren: Computing tropical points over fields with valuation
#
################################################################################

@doc raw"""
    tropical_variety(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false, skip_saturation::Bool=false, skip_primary_decomposition::Bool=false)

Return the tropicalization of `I` with respect to `nu` as a `Vector{TropicalVariety}`.
If `nu==nothing`, will compute with respect to the trivial valuation and min convention.
If `weighted_polyhedral_complex_only==true`, will not cache any additional information.
If `skip_saturation==true`, will not saturate `I` with respect to the product of all variables.
If `skip_primary_decomposition==true`, will not decompose `I`.

!!! warning
    `tropical_variety` is currently under development and only works for ideals that primary decompose into
    principal, linear, and binomial ideals.

# Examples
```jldoctest
julia> R,(x,y) = QQ[:x, :y];

julia> I = ideal([(x^2+y)*(x+y^2)*(x+y)]);

julia> tropical_variety(I)
3-element Vector{TropicalVariety}:
 Min tropical variety
 Min tropical variety
 Min tropical variety

```
"""
function tropical_variety(I::Union{MPolyIdeal,MPolyRingElem}, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false, skip_saturation::Bool=false, skip_primary_decomposition::Bool=false)
    ###
    # Step 0.a: convert I to ideal if poly,
    #   initialize nu as the trivial valuation if not specified by user
    ###
    if I isa MPolyRingElem
        I = ideal(parent(I),[I])
    end
    if isnothing(nu)
        nu = tropical_semiring_map(coefficient_ring(I))
    end

    ###
    # Step 0.b: Saturate `I` and assert that `I` is not the whole ring
    ###
    if !skip_saturation
        R = base_ring(I)
        I = saturation(I,ideal([prod(gens(R))]))
    end
    @req !isone(I) "ideal contains a monomial, tropical varieties in OSCAR cannot be empty"

    ###
    # Step 0.c: Compute primary decomposition and tropical varieties of all primary factors
    ###
    toTropicalize = [I]
    if !skip_primary_decomposition
        toTropicalize = [ P for (P,_) in primary_decomposition(I) ]
    end

    tropicalVarieties = TropicalVariety[]
    for P in toTropicalize
        # compute a reduced GB to test whether `P` is principal, binomial, or linear
        GB = groebner_basis(P,complete_reduction=true)

        if length(GB)==1
            # P is principal
            TropV = tropical_variety_principal(P,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
        elseif max(length.(GB)...)==2
            # P is binomial
            TropV = tropical_variety_binomial(P,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
        elseif max(total_degree.(GB)...)==1
            # P linear
            TropV = tropical_variety_linear(P,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)
        else
            # P general
            error("general tropical varieties currently unsupported")
        end
        push!(tropicalVarieties,TropV)
    end

    return tropicalVarieties
end


################################################################################
#
#  Principal ideals
#
################################################################################

function tropical_variety_principal(I::MPolyIdeal,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)
    ###
    # Construct TropicalVariety from TropicalHypersurface
    ###
    g = first(gens(I))
    TropV = tropical_variety(tropical_hypersurface(g,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only))
    if !weighted_polyhedral_complex_only
        set_attribute!(TropV,:algebraic_ideal,I)
        set_attribute!(TropV,:tropical_semiring_map,nu)
    end
    return TropV
end



################################################################################
#
#  Binomial ideals
#
################################################################################

function tropical_variety_binomial(I::MPolyIdeal,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)
    ###
    # Construct matrix of exponent vector differences
    # and vector of coefficient valuation differences
    ###
    G = gens(I)
    A = matrix(ZZ,[first(collect(expv))-last(collect(expv)) for expv in exponents.(G)])
    b = [ QQ(nu(last(collect(coeff))/first(collect(coeff)))) for coeff in coefficients.(G)]

    ###
    # Compute tropical variety multiplicity
    ###
    snfAdiag = elementary_divisors(A)
    weight = abs(prod([m for m in snfAdiag if !iszero(m)]))

    ###
    # Constructing tropical variety set-theoretically
    ###
    A = QQMatrix(A)
    L = transpose(kernel(A, side = :right))
    can_solve, V = can_solve_with_solution(transpose(A),matrix(QQ,[b]),side=:left)
    @req can_solve "tropical variety cannot be empty"
    SigmaV = polyhedral_complex(IncidenceMatrix([[1]]), V, nothing, L)

    ###
    # Assemble tropical variety
    ###
    TropV = tropical_variety(SigmaV,[weight],convention(nu))
    if !weighted_polyhedral_complex_only
        set_attribute!(TropV,:algebraic_ideal,I)
        set_attribute!(TropV,:tropical_semiring_map,nu)
    end
    return TropV
end



################################################################################
#
#  Linear ideals
#
################################################################################

function tropical_variety_linear(I::MPolyIdeal,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)
    ###
    # Compute reduced Groebner basis (usually already cached),
    # and check whether the linear polynomials have a constant term
    ###
    R = base_ring(I)
    G = groebner_basis(I,complete_reduction=true)
    if min(total_degree.(Iterators.flatten(collect.(terms.(G))))...)==1
        # input homogneeous, construct TropicalVariety via TropicalLinearSpace
        TropV = tropical_variety(tropical_linear_space(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only))
        if !weighted_polyhedral_complex_only
            set_attribute!(TropV,:algebraic_ideal,I)
            set_attribute!(TropV,:tropical_semiring_map,nu)
        end
        return TropV
    else
        # input inhomogeneous, homogenise first
        Ih,_,_ = homogenize_pre_tropicalization(I)
        TropLh = tropical_linear_space(Ih,nu,weighted_polyhedral_complex_only=true)
        Sigma = dehomogenize_post_tropicalization(polyhedral_complex(TropLh))

        multiplicities = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
        TropV = tropical_variety(Sigma,multiplicities)
        if !weighted_polyhedral_complex_only
            set_attribute!(TropV,:algebraic_ideal,I)
            set_attribute!(TropV,:tropical_semiring_map,nu)
        end
        return TropV
    end
end



################################################################################
#
#  Zero-dimensional ideals
#
################################################################################

@doc raw"""
    Oscar.tropical_variety_zerodimensional(I::MPolyIdeal,nu::TropicalSemiringMap{QQField,ZZRingElem,<:Union{typeof(min),typeof(max)}}; precision::Int=64)

Internal function for computing zero-dimensional tropical varieties over p-adic numbers via
finite precision Eigenvalue computation.  Assumes without test that `I` is zero-dimensional.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> R,(x1,x2,x3) = polynomial_ring(QQ,3);

julia> I = ideal([28*x3^2 - 1*x3 - 1,
                  2*x2 - x3,
                  2*x1 - x2]);

julia> nu = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> TropI = Oscar.tropical_variety_zero_dimensional(I,nu)
Min tropical variety

julia> vertices(TropI)
2-element SubObjectIterator{PointVector{QQFieldElem}}:
 [-2, -1, 0]
 [-4, -3, -2]

julia> nu = tropical_semiring_map(QQ,3,max)
Map into Max tropical semiring encoding the 3-adic valuation on Rational field

julia> TropI = Oscar.tropical_variety_zero_dimensional(I,nu)
Max tropical variety

julia> vertices(TropI)
1-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0, 0]

```
"""
function tropical_variety_zero_dimensional(I::MPolyIdeal,nu::TropicalSemiringMap{QQField,ZZRingElem,<:Union{typeof(min),typeof(max)}}; precision::Int=64)
    # Construct the representation matrices of the multiplications by xi in K[x]/I
    _,x = number_field(I)
    mx = representation_matrix.(x)

    # Compute their simultaneous diagonalization numerically
    Qp = padic_field(uniformizer(nu), precision=precision)
    TropVDict = Oscar.simultaneous_diagonalization(map_entries.(Ref(Qp), mx))

    # Construct their tropical variety as a polyhedral complex consisting only of vertices
    # and a list of multiplicities
    TropVPoints = convention(nu)==min ? collect(values(TropVDict)) : -collect(values(TropVDict))
    TropVPointsUnique = unique(TropVPoints)
    Sigma = polyhedral_complex(IncidenceMatrix([[i] for i in 1:length(TropVPointsUnique)]), TropVPointsUnique)
    TropVMults = [ZZ(length(findall(isequal(p),TropVPoints))) for p in TropVPointsUnique]
    TropV = tropical_variety(Sigma,TropVMults,convention(nu))
    set_attribute!(TropV,:algebraic_points,collect(keys(TropVDict)))
    return TropV
end


###
# Code by Claus:
###
function slope_eigenspace(M::MatElem{T}) where T <: Hecke.NonArchLocalFieldElem
    f = charpoly(M)
    lf = Hecke.slope_factorization(f)
    # @req all(==(1), values(lf))

    se = Dict{typeof(f), typeof(M)}()
    k = base_ring(M)
    zk = maximal_order(k)

    for f = keys(lf)
        se[f] = kernel(f(M), side = :right) #hopefully, this is in rref
    end
    @assert sum(ncols(x) for x = values(se)) == nrows(M)
    return se
end

function _intersect(M::MatElem{T}, N::MatElem{T}) where T <: Hecke.FieldElem
    k = base_ring(M)
    I = [M N]
    PR = maximum(precision, I)
    pr = minimum(precision, I)
    if pr != PR
        for i = eachindex(I)
            I[i] = setprecision(I[i], pr)
        end
    end

    v = kernel(I, side = :right) #precision issues...
    l = M*v[1:ncols(M), 1:ncols(v)]
    return transpose(rref(transpose(l))[2])
end

function valuation_of_roots(f::PolyRingElem{<:Hecke.NonArchLocalFieldElem})
    iszero(f) && error("polynomial must not be zero")
    return (valuation(constant_coefficient(f)) - valuation(leading_coefficient(f)))//degree(f)
end

function simultaneous_diagonalization(v::Vector{<:MatElem{T}}) where T <: Hecke.NonArchLocalFieldElem

    k = base_ring(v[1])
    @assert all(x->base_ring(x) == k, v)
    n = nrows(v[1])
    @assert all(x->ncols(x) == nrows(x) == n, v)

    vv = map(slope_eigenspace, v)

    d = Dict(v => [valuation_of_roots(k)] for (k,v) = vv[1])
    @assert sum(ncols(x) for x = keys(d)) == n
    for i=2:length(vv)
        dd = typeof(d)()
        for (mat, pol_vec) = d
            for (p, m) = vv[i]
                j = _intersect(mat, m)
                if ncols(j) > 0
                    dd[j] = push!(copy(pol_vec), valuation_of_roots(p))
                end
            end
        end
        d = dd
        @assert sum(ncols(x) for x = keys(d)) == n
    end

    return d
end


# #=======
# tropical variety of an ideal
# todo: proper documentation
# Example:
# import Random
# K,s = RationalFunctionField(QQ,"s");
# Kx,(x1,x2,x3,x4) = polynomial_ring(K,4);
# val = TropicalSemiringMap(K,s);
# I = ideal([x1-s*x2+(s+1)*x3,3*x2-s^2*x3+(s^2+1)*x4]);
# Random.seed!(3847598273423);
# TropI = tropical_variety(I,val)
# =======#
# function tropical_variety(I::MPolyIdeal, val::TropicalSemiringMap, convention::Union{typeof(min),typeof(max)}=min)

#   ###
#   # Part 0: Preprocessing
#   #   Check whether valuation is on the coefficient ring of input polynomials,
#   #   homogenize input ideal if not homogeneous
#   ###
#   if coefficient_ring(base_ring(I))!=val.valued_field
#     error("input valuation not on coefficient ring of input ideal")
#   end
#   was_input_homogeneous = true
#   for g in groebner_basis(I,complete_reduction=true) # todo: replace GB computation with interreduction
#     if !sloppy_is_homogeneous(g)
#       was_input_homogeneous = false
#       I = homogenize(I)
#       break
#     end
#   end



#   ###
#   # Part 1: Tropical starting polyhedra (todo: avoid recomputation)
#   #   - compute and recompute starting points until they lie in the relative interior of maximal cells
#   #   - initialize working lists
#   ###

#   # Note: a working list entry consists of a triple (w,C,G) where
#   #   * G is a (tropical) Groebner basis
#   #   * C is the Groebner polyhedron
#   #   * w is the sum of vertices and rays of C
#   # In particular:
#   #   * w is a weight vector with respect to which G is a Groebner basis,
#   #   * w is compatible with coordinate permutations if symmetries exist,
#   #   * instead of comparing C or G it suffices to compare w.
#   working_list_todo = [] # list of groebner polyhedra with potentially unknown neighbors
#   working_list_done = [] # list of groebner polyhedra with known neighbors
#   facet_points_done = [] # list of facet points whose tropical links were computed and traversed

#   compute_starting_points = true
#   while compute_starting_points
#     print("computing random starting points... ")
#     starting_points = tropical_points(I,val)
#     println("done")

#     working_list_todo = []
#     working_list_done = []
#     facet_points_done = []
#     compute_starting_points = false
#     for starting_point in starting_points
#       print("computing groebner_basis for starting point ",starting_point,"... ")
#       G = groebner_basis(I,val,starting_point)
#       println("done")
#       C = groebner_polyhedron(G,val,starting_point)
#       w = anchor_point(C)

#       # if C is lower-dimensional, recompute all starting points
#       if dim(C)!=dim(I)
#         println("starting point on lower-dimensional cell, recomputing...")
#         compute_starting_points = true
#         break
#       end

#       # if (w,C,G) already is in working_list_todo, skip
#       i = searchsortedfirst(working_list_todo,(w,C,G),by=first)
#       if i<=length(working_list_todo) && working_list_todo[i][1]==w
#         continue
#       end
#       # otherwise, add (w,C,G) to todo list
#       insert!(working_list_todo, i, (w,C,G))
#     end
#   end



#   ###
#   # Part 2: Tropical traversal
#   ###
#   while !isempty(working_list_todo)
#     print("#working_list_todo: ",length(working_list_todo),"  ")
#     println("#working_list_done: ",length(working_list_done))

#     # pick a groebner polyhedron from todo list, add it to the done list, and compute its facet points
#     (w,C,G) = popfirst!(working_list_todo)
#     i = searchsortedfirst(working_list_done,(w,C,G),by=first)
#     insert!(working_list_done, i, (w,C,G))

#     points_to_traverse = facet_points(C)
#     for point_to_traverse in points_to_traverse
#       # if point was traversed before, skip
#       i = searchsortedfirst(facet_points_done,point_to_traverse)
#       if i<=length(facet_points_done) && facet_points_done[i]==point_to_traverse
#         continue
#       end
#       # otherwise add point_to_traverse to facet_points_done
#       insert!(facet_points_done, i, point_to_traverse)

#       directions_to_traverse = tropical_link(ideal(G),val,point_to_traverse) # todo, this output can be wrong
#       for direction_to_traverse in directions_to_traverse
#         # compute neighbor
#         print("computing groebner_basis for ",point_to_traverse,direction_to_traverse,"... ")
#         G_neighbor = groebner_flip(G,val,w,point_to_traverse,direction_to_traverse)
#         println("done")
#         C_neighbor = groebner_polyhedron(G_neighbor,val,point_to_traverse,perturbation=direction_to_traverse)
#         w_neighbor = anchor_point(C_neighbor)

#         # if neighbor is already in done list, skip
#         i = searchsortedfirst(working_list_done,
#                               (w_neighbor,C_neighbor,G_neighbor),
#                               by=first)
#         if i<=length(working_list_done) && working_list_done[i][1]==w_neighbor
#           continue
#         end
#         # if neighbor is already in todo list, skip
#         i = searchsortedfirst(working_list_todo,
#                               (w_neighbor,C_neighbor,G_neighbor),
#                               by=first)
#         if i<=length(working_list_todo) && working_list_todo[i][1]==w_neighbor
#           continue
#         end
#         # otherwise, add neighbor to todo list
#         insert!(working_list_todo, i, (w_neighbor,C_neighbor,G_neighbor))
#       end
#     end
#   end



#   ###
#   # Part 3: Postprocessing
#   ###
#   # 3.0: dehomogenize data if input was homogenized
#   if !was_input_homogeneous
#     n = length(gens(base_ring(I)))

#     # 3.0.1: dehomogenize Groebner polyhedra
#     zeroth_unit_vector_as_row_vector = zeros(Int,1,n)
#     zeroth_unit_vector_as_row_vector[1,1] = 1
#     dehomogenising_hyperplane = polyhedron((zeros(Int,0,n),zeros(Int,0)),
#                                            (zeroth_unit_vector_as_row_vector,[1]))
#     for wCG in working_list_done
#       wCG[2] = intersect(wCG[2],dehomogenising_hyperplane)
#     end
#     # 3.0.2: check that initial ideals are distinct (todo)
#   end

#   # 3.1: construct PolyhedralComplex
#   # 3.1.1: construct incidence_matrix, vertices_and_rays, and far_vertices
#   incidence_matrix = Vector{Vector{Int}}()
#   verts_rays = Vector{Vector{Polymake.Rational}}()
#   far_vertices = Vector{Int}()
#   for (w,C,G) in working_list_done
#     incidence_vector = Vector{Int}()
#     for vert in vertices(C)
#       i = findfirst(isequal(vert),verts_rays)
#       if i === nothing
#         # if vert does not occur in verts_rays
#         # add it to verts_rays
#         push!(verts_rays,vert)
#         push!(incidence_vector,length(verts_rays))
#       else
#         push!(incidence_vector,i)
#       end
#     end
#     for ray in rays(C)
#       i = findfirst(isequal(ray),verts_rays)
#       if i === nothing || !(i in far_vertices)
#         # if ray does not occur in verts_rays or if it occurs but not as a ray,
#         # add it to verts_rays
#         push!(verts_rays,ray)
#         push!(far_vertices,length(verts_rays))
#         push!(incidence_vector,length(verts_rays))
#       else
#         push!(incidence_vector,i)
#       end
#     end
#     push!(incidence_matrix,incidence_vector)

#   end
#   verts_rays_matrix = permutedims(reduce(hcat, verts_rays)) # convert Vector{Vector} to Matrix

#   # 3.1.2: construct lineality space
#   (w,C,G) = first(working_list_done)
#   lineality_space_gens = matrix(QQ,lineality_space(C))



#   # 3.2: Construct lists for weight_vectors, initial_ideals and multiplicities
#   weight_vectors = [w for (w,C,G) in working_list_done]
#   initial_ideals = [ideal(initial(G,val,w)) for (w,C,G) in working_list_done]
#   multiplicities = [multiplicity(inI) for inI in initial_ideals]



#   PC = polyhedral_complex(IncidenceMatrix(incidence_matrix),
#   verts_rays_matrix,
#   far_vertices,
#   lineality_space_gens)

#   verts_rays_perm  = collect(vertices_and_rays(PC))
#   verts_rays_perm = Vector{Int64}.(verts_rays_perm)
#   permutation = [findfirst(isequal(l), verts_rays_perm) for l in verts_rays]
# #  inc = [findall(incidence_matrix[i, :]) for i in 1:size(incidence_matrix, 1)]
#   new_incidence = [permutation[incidence] for incidence in incidence_matrix]
#   mults = Dict(new_incidence[i] => multiplicities[i] for i in 1:length(working_list_done))

#   TropI = TropicalVariety{typeof(max),true}(PC, mults)


#   set_attribute!(TropI,:weight_vectors,weight_vectors)
#   set_attribute!(TropI,:initial_ideals,initial_ideals)

#   return TropI
# end



# #=======
# Example:
# P = cube(4)
# anchor_point(P)
# facet_points(P)
# =======#
# function anchor_point(P::Polyhedron)
#   # compute the sum of vertices and rays in homogenized coordinates
#   pt = convert(Vector{QQFieldElem},sum([vertices(P)...,rays(P)...]))
#   pushfirst!(pt,n_vertices(P))

#   # project to orthogonal complement of lineality space if necessary
#   if lineality_dim(P)>0
#     pt = Polymake.Matrix{Polymake.Rational}(vcat(transpose(pt)))
#     Polymake.common.project_to_orthogonal_complement(pt, P.pm_polytope.LINEALITY_SPACE)
#     pt = convert(Matrix{QQFieldElem}, pt)[1,:]
#   end

#   # rescale until first entry is 1 and remove it
#   pt = [pt[i]//pt[1] for i in 2:length(pt)]
#   return pt
# end

# function facet_points(P::Polyhedron)
#   points = []
#   for facet in faces(P,dim(P)-1)
#     if length(vertices(facet))>0 # skipping facets at infinity
#       push!(points,anchor_point(facet))
#     end
#   end
#   return points
# end
