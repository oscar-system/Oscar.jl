###
# Tropical varieties in Oscar
# ===========================
###



###
# 0. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention, affecting initial ideals
# EMB = true or false:
#   embedded or abstract tropical variety
#   embedded tropical variety = weighted polyhedral complex in euclidean space
#   abstract tropical variety = weighted hypergraph with enumerated vertices
###

@attributes mutable struct TropicalVariety{M,EMB} <: TropicalVarietySupertype{M,EMB}
    polyhedralComplex::PolyhedralComplex
    multiplicities::Dict{Vector{Int}, Int}
    function TropicalVariety{M,EMB}(Sigma::PolyhedralComplex) where {M,EMB}
        return new{M,EMB}(Sigma)
    end
end
export TropicalVariety


function pm_object(T::TropicalVariety)
    if has_attribute(T,:polymake_bigobject)
        return get_attribute(T,:polymake_bigobject)
    end
    error("pm_object(T::TropicalVariety): Has no polymake bigobject.")
end



###
# 1. Printing
# -----------
###

function Base.show(io::IO, tv::TropicalVariety{M, EMB}) where {M, EMB}
    if EMB
        print(io, "A $(repr(M)) tropical variety of dimension $(dim(tv)) embedded in $(ambient_dim(tv))-dimensional Euclidian space")
    else
        print(io, "An abstract $(repr(M)) tropical variety of dimension $(dim(tv))")
    end
end


###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalVariety()

Construct the embedded tropical variety of a polynomial ideal over a (possibly trivially) valued field

# Examples
"""
# todo: Dartmouth
# function TropicalVariety()
#
#     return #...
# end



@doc Markdown.doc"""
    TropicalVariety{M,EMB}(Sigma::PolyhedralComplex)

Construct the abstract tropical variety from a polyhedral complex

# Examples
```jldoctest
julia> IM = IncidenceMatrix([[1,2],[1,3],[1,4]]);

julia> VR = [0 0; 1 0; 0 1; -1 -1];

julia> far_vertices = [2,3,4];

julia> Sigma = PolyhedralComplex(IM, VR, far_vertices);

julia> tropicalLine = TropicalVariety{min,true}(Sigma)
"""


###
# 3. Basic properties
# -------------------
###



###
# 4. Tropical varieties of polynomial ideals
# ------------------------------------------
# References for computing tropical varieties via traversal:
#   T. Bogart, A. Jensen, D. Speyer, B. Sturmfels, R. Thomas: Computing tropical varieties
#   T. Markwig, Y. Ren: Computing tropical points over fields with valuation
# Reference for the general structure of tropical varieties:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
###

function homogenize(I::MPolyIdeal)
  G = groebner_basis(I,complete_reduction=true)

  Kx = base_ring(I)
  K = coefficient_ring(Kx)
  x = symbols(Kx)
  Kxh,_ = PolynomialRing(K,vcat([:xh],x))

  Gh = Vector{elem_type(Kx)}(undef,length(G))
  for (i,g) in enumerate(G)
    gh = MPolyBuildCtx(Kxh)
    d = max([sum(expv) for expv in exponent_vectors(g)]...) # degree of g
    for (c,alpha) in zip(coefficients(g),exponent_vectors(g))
      pushfirst!(alpha,d-sum(alpha)) # homogenize exponent vector
      push_term!(gh,c,alpha)
    end
    Gh[i] = finish(gh)
  end

  return ideal(Gh)
end
export homogenize


#=======
tropical variety of an ideal
todo: proper documentation
Example:
import Random
K,s = RationalFunctionField(QQ,"s");
Kx,(x1,x2,x3,x4) = PolynomialRing(K,4);
val = ValuationMap(K,s);
I = ideal([x1-s*x2+(s+1)*x3,3*x2-s^2*x3+(s^2+1)*x4]);
Random.seed!(3847598273423);
TropI = tropical_variety(I,val)
=======#
function tropical_variety(I::MPolyIdeal, val::ValuationMap, convention::Union{typeof(min),typeof(max)}=min)

  ###
  # Part 0: Preprocessing
  #   Check whether valuation is on the coefficient ring of input polynomials,
  #   homogenize input ideal if not homogeneous
  ###
  if coefficient_ring(base_ring(I))!=val.valued_field
    error("input valuation not on coefficient ring of input ideal")
  end
  was_input_homogeneous = true
  for g in groebner_basis(I,complete_reduction=true) # todo: replace GB computation with interreduction
    if !sloppy_is_homogeneous(g)
      was_input_homogeneous = false
      I = homogenize(I)
      break
    end
  end



  ###
  # Part 1: Tropical starting polyhedra (todo: avoid recomputation)
  #   - compute and recompute starting points until they lie in the relative interior of maximal cells
  #   - initialize working lists
  ###

  # Note: a working list entry consists of a triple (w,C,G) where
  #   * G is a (tropical) Groebner basis
  #   * C is the Groebner polyhedron
  #   * w is the sum of vertices and rays of C
  # In particular:
  #   * w is a weight vector with respect to which G is a Groebner basis,
  #   * w is compatible with coordinate permutations if symmetries exist,
  #   * instead of comparing C or G it suffices to compare w.
  working_list_todo = [] # list of groebner polyhedra with potentially unknown neighbours
  working_list_done = [] # list of groebner polyhedra with known neighbours
  facet_points_done = [] # list of facet points whose tropical links were computed and traversed

  compute_starting_points = true
  while compute_starting_points
    print("computing random starting points... ")
    starting_points = tropical_points(I,val)
    println("done")

    working_list_todo = []
    working_list_done = []
    facet_points_done = []
    compute_starting_points = false
    for starting_point in starting_points
      print("computing groebner_basis for starting point ",starting_point,"... ")
      G = groebner_basis(I,val,starting_point)
      println("done")
      C = groebner_polyhedron(G,val,starting_point)
      w = anchor_point(C)

      # if C is lower-dimensional, recompute all starting points
      if dim(C)!=dim(I)
        println("starting point on lower-dimensional cell, recomputing...")
        compute_starting_points = true
        break
      end

      # if (w,C,G) already is in working_list_todo, skip
      i = searchsortedfirst(working_list_todo,(w,C,G),by=x->x[1])
      if i<=length(working_list_todo) && working_list_todo[i][1]==w
        continue
      end
      # otherwise, add (w,C,G) to todo list
      insert!(working_list_todo, i, (w,C,G))
    end
  end



  ###
  # Part 2: Tropical traversal
  ###
  while !isempty(working_list_todo)
    print("#working_list_todo: ",length(working_list_todo),"  ")
    println("#working_list_done: ",length(working_list_done))

    # pick a groebner polyhedron from todo list, add it to the done list, and compute its facet points
    (w,C,G) = popfirst!(working_list_todo)
    i = searchsortedfirst(working_list_done,(w,C,G),by=x->x[1])
    insert!(working_list_done, i, (w,C,G))

    points_to_traverse = facet_points(C)
    for point_to_traverse in points_to_traverse
      # if point was traversed before, skip
      i = searchsortedfirst(facet_points_done,point_to_traverse)
      if i<=length(facet_points_done) && facet_points_done[i]==point_to_traverse
        continue
      end
      # otherwise add point_to_traverse to facet_points_done
      insert!(facet_points_done, i, point_to_traverse)

      directions_to_traverse = tropical_link(ideal(G),val,point_to_traverse) # todo, this output can be wrong
      for direction_to_traverse in directions_to_traverse
        # compute neighbour
        print("computing groebner_basis for ",point_to_traverse,direction_to_traverse,"... ")
        G_neighbour = groebner_flip(G,val,w,point_to_traverse,direction_to_traverse)
        println("done")
        C_neighbour = groebner_polyhedron(G_neighbour,val,point_to_traverse,pertubation=direction_to_traverse)
        w_neighbour = anchor_point(C_neighbour)

        # if neighbour is already in done list, skip
        i = searchsortedfirst(working_list_done,
                              (w_neighbour,C_neighbour,G_neighbour),
                              by=x->x[1])
        if i<=length(working_list_done) && working_list_done[i][1]==w_neighbour
          continue
        end
        # if neighbour is already in todo list, skip
        i = searchsortedfirst(working_list_todo,
                              (w_neighbour,C_neighbour,G_neighbour),
                              by=x->x[1])
        if i<=length(working_list_todo) && working_list_todo[i][1]==w_neighbour
          continue
        end
        # otherwise, add neighbour to todo list
        insert!(working_list_todo, i, (w_neighbour,C_neighbour,G_neighbour))
      end
    end
  end



  ###
  # Part 3: Postprocessing
  ###
  # 3.0: dehomogenize data if input was homogenized
  if !was_input_homogeneous
    n = length(gens(base_ring(I)))

    # 3.0.1: dehomogenize Groebner polyhedra
    zeroth_unit_vector_as_row_vector = zeros(Int,1,n)
    zeroth_unit_vector_as_row_vector[1,1] = 1
    dehomogenising_hyperplane = Polyhedron((zeros(Int,0,n),zeros(Int,0)),
                                           (zeroth_unit_vector_as_row_vector,[1]))
    for wCG in working_list_done
      wCG[2] = intersect(wCG[2],dehomogenising_hyperplane)
    end
    # 3.0.2: check that initial ideals are distinct (todo)
  end

  # 3.1: construct PolyhedralComplex
  # 3.1.1: construct incidence_matrix, vertices_and_rays, and far_vertices
  incidence_matrix = Vector{Vector{Int}}()
  vertices_and_rays = Vector{Vector{Polymake.Rational}}()
  far_vertices = Vector{Int}()
  for (w,C,G) in working_list_done
    incidence_vector = Vector{Int}()
    for vert in vertices(C)
      i = findfirst(isequal(vert),vertices_and_rays)
      if i == nothing
        # if vert does not occur in vertices_and_rays
        # add it to vertices_and_rays
        push!(vertices_and_rays,vert)
        push!(incidence_vector,length(vertices_and_rays))
      else
        push!(incidence_vector,i)
      end
    end
    for ray in rays(C)
      i = findfirst(isequal(ray),vertices_and_rays)
      if i == nothing || !(i in far_vertices)
        # if ray does not occur in vertices_and_rays or if it occurs but not as a ray,
        # add it to vertices_and_rays
        push!(vertices_and_rays,ray)
        push!(far_vertices,length(vertices_and_rays))
        push!(incidence_vector,length(vertices_and_rays))
      else
        push!(incidence_vector,i)
      end
    end
    push!(incidence_matrix,incidence_vector)

  end
  vertices_and_rays = permutedims(reduce(hcat, vertices_and_rays)) # convert Vector{Vector} to Matrix

  # 3.1.2: construct lineality space
  (w,C,G) = first(working_list_done)
  lineality_space_gens = matrix(QQ,lineality_space(C))



  # 3.2: Construct lists for weight_vectors, initial_ideals and multiplicities
  weight_vectors = [w for (w,C,G) in working_list_done]
  initial_ideals = [ideal(initial(G,val,w)) for (w,C,G) in working_list_done]
  multiplicities = [multiplicity(inI) for inI in initial_ideals]

  mults = Dict(incidence_matrix[i] => multiplicities[i] for i in 1:length(working_list_done))
  TropI = TropicalVariety{typeof(max),true}(PolyhedralComplex(IncidenceMatrix(incidence_matrix),
                                                              vertices_and_rays,
                                                              far_vertices,
                                                              lineality_space_gens),
					    mults)





  set_attribute!(TropI,:weight_vectors,weight_vectors)
  set_attribute!(TropI,:initial_ideals,initial_ideals)

  return TropI
end
export tropical_variety



#=======
Example:
P = cube(4)
anchor_point(P)
facet_points(P)
=======#
function anchor_point(P::Polyhedron)
  # compute the sum of vertices and rays in homogenized coordinates
  pt = convert(Vector{fmpq},sum([vertices(P)...,rays(P)...]))
  pushfirst!(pt,nvertices(P))

  # project to orthogonal complement of lineality space if necessary
  if lineality_dim(P)>0
    pt = Polymake.Matrix{Polymake.Rational}(vcat(transpose(pt)))
    Polymake.common.project_to_orthogonal_complement(pt, P.pm_polytope.LINEALITY_SPACE)
    pt = convert(Matrix{fmpq}, pt)[1,:]
  end

  # rescale until first entry is 1 and remove it
  pt = [pt[i]//pt[1] for i in 2:length(pt)]
  return pt
end
export anchor_point

function facet_points(P::Polyhedron)
  points = []
  for facet in faces(P,dim(P)-1)
    if length(vertices(facet))>0 # skipping facets at infinity
      push!(points,anchor_point(facet))
    end
  end
  return points
end
export facet_points
