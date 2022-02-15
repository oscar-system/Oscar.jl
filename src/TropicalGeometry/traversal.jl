###
# Traversing on tropical varieties
# ================================
#
# References for computing tropical varieties via traversal:
#   T. Bogart, A. Jensen, D. Speyer, B. Sturmfels, R. Thomas: Computing tropical varieties
#   T. Markwig, Y. Ren: Computing tropical points over fields with valuation
# Reference for the general structure of tropical varieties:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
###


#=======
tropical variety of an ideal
todo: proper documentation
Example:
import Random
K,s = RationalFunctionField(QQ,"s");
Kx,(x1,x2,x3,x4) = PolynomialRing(K,4);
val = ValuationMap(K,s);
I = ideal([x1-s*x2+(s+1)*x3,3*x2-s^2*x3+(s^2+1)*x4]);
Random.seed!(13371337);
tropical_variety(I,val)
=======#
function tropical_variety(I::MPolyIdeal, val::ValuationMap)

  ###
  # Part 1: Computing starting points and initializing working lists
  ###
  print("computing starting points... ")
  starting_points = tropical_points(I,val)
  println("done")

  # Note: in this function, a Groebner polyhedra means a triple (w,C,G) where
  #   * G is a (tropical) Groebner basis
  #   * C is the actual Groebner polyhedron
  #   * w is the sum of vertices and rays of C
  # In particular, w is a weight vector with respect to which G is a Groebner basis,
  #   w is compatible with coordinate permutations if symmetries exist,
  #   and instead of comparing C it suffices to compare w.
  maximal_cells_todo = [] # list of groebner polyhedra with potentially unknown neighbours
  maximal_cells_done = [] # list of groebner polyhedra with known neighbours
  facets_done = []        # list of facet points whose tropical links were computed and traversed
  for starting_point in starting_points
    G = groebner_basis(I,val,starting_point)
    C = groebner_polyhedron(I,val,starting_point,skip_groebner_basis_computation=true)
    w = anchor_point(C)
    push!(maximal_cells_todo, (w,C,G))
  end

  ###
  # Part 2: Traversing the tropical variety
  ###
  while !isempty(maximal_cells_todo)
    println("length(maximal_cells_todo): ",length(maximal_cells_todo))
    println("length(maximal_cells_done): ",length(maximal_cells_done))

    # pick a groebner polyhedron from todo list, add it to the done list, and compute its facet points
    (w,C,G) = popfirst!(maximal_cells_todo)
    push!(maximal_cells_done, (w,C,G))
    points_to_traverse = facet_points(C)
    println("points_to_traverse:")
    println(points_to_traverse)

    for point_to_traverse in points_to_traverse
      # if point was traversed before, skip
      i = searchsortedfirst(facets_done,point_to_traverse)
      if i<=length(facets_done) && facets_done[i]==point_to_traverse
        continue
      end
      # otherwise add point_to_traverse to facets_done
      insert!(facets_done, i, point_to_traverse)

      # compute tropical link
      directions_to_traverse = tropical_link(ideal(G),val,point_to_traverse)
      println("directions_to_traverse:")
      println(directions_to_traverse)

      for direction_to_traverse in directions_to_traverse
        # compute neighbour
        println("groebner_flip")
        G_neighbor = groebner_flip(G,val,w,point_to_traverse,direction_to_traverse)
        println("groebner_polyhedron")
        C_neighbor = groebner_polyhedron(G_neighbor,val,point_to_traverse,pertubation=direction_to_traverse)
        w_neighbor = anchor_point(C_neighbor)

        # if neighbour is in done list, skip
        i = searchsortedfirst(maximal_cells_done,
                              (w_neighbor,C_neighbor,G_neighbor),
                              by=x->x[1])
        if i<=length(maximal_cells_done) && maximal_cells_done[i][1]==w_neighbor
          continue
        end

        # if neighbour is in todo list, skip
        i = searchsortedfirst(maximal_cells_todo,
                              (w_neighbor,C_neighbor,G_neighbor),
                              by=x->x[1])
        if i<=length(maximal_cells_todo) && maximal_cells_todo[i][1]==w_neighbor
          continue
        end
        # otherwise, add data to todo list
        insert!(maximal_cells_todo, i, (w_neighbor,C_neighbor,G_neighbor))
      end
    end
  end

  ###
  # Part 3: Preparing data to return
  ###
  weights_and_groebner_bases = []
  incidence_matrix = Vector{Vector{Int}}()
  vertices_and_rays = Vector{Vector{Polymake.Rational}}()
  far_vertices = Vector{Int}()
  for (w,C,G) in maximal_cells_done
    push!(weights_and_groebner_bases,(w,G))

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

  println("incidence_matrix:")
  println(incidence_matrix)
  println("vertices_and_rays:")
  println(vertices_and_rays)
  Trop_I = PolyhedralComplex(IncidenceMatrix(incidence_matrix),
                             vertices_and_rays,
                             far_vertices)

  return Trop_I, weights_and_groebner_bases
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
