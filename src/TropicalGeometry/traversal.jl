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
  print("computing starting points... ")
  starting_points = tropical_points(I,val)
  println("done")

  groebner_polyhedra_todo = []
  groebner_polyhedra_done = []
  links_done = []
  for starting_point in starting_points
    G = groebner_basis(I,val,starting_point)
    C = groebner_polyhedron(I,val,starting_point,skip_groebner_basis_computation=true)
    w = anchor_point(C)
    push!(groebner_polyhedra_todo, (w,C,G))
  end

  while !isempty(groebner_polyhedra_todo)
    println("length(groebner_polyhedra_todo): ",length(groebner_polyhedra_todo))
    println("length(groebner_polyhedra_done): ",length(groebner_polyhedra_done))

    (w,C,G) = popfirst!(groebner_polyhedra_todo)
    points_to_traverse = facet_points(C)
    println("points_to_traverse:")
    println(points_to_traverse)
    for point_to_traverse in points_to_traverse

      # skip facet_point if it lies in links_done
      i = searchsortedfirst(links_done,point_to_traverse)
      if i<=length(links_done) && links_done[i]==point_to_traverse # todo: is facet_point the sum of the vertices?
        continue
      end
      # add facet_point to links_done
      insert!(links_done, i, point_to_traverse)

      # compute tropical link
      directions_to_traverse = tropical_link(ideal(G),val,point_to_traverse)
      println("directions_to_traverse:")
      println(directions_to_traverse)
      for direction_to_traverse in directions_to_traverse
        # compute neighboring data
        println("groebner_flip")
        G_neighbor = groebner_flip(G,val,w,point_to_traverse,direction_to_traverse)
        println("groebner_polyhedron")
        C_neighbor = groebner_polyhedron(G_neighbor,val,point_to_traverse,pertubation=direction_to_traverse) # todo: this computes an unnecessary GB
        println("summing vertices")
        w_neighbor = anchor_point(C)

        # if neighboring polyhedra is in done list, skip
        i = searchsortedfirst(groebner_polyhedra_done,
                              (w_neighbor,C_neighbor,G_neighbor),
                              by=x->x[1])
        if i<=length(groebner_polyhedra_done) && groebner_polyhedra_done[i][1]==w_neighbor
          continue
        end

        # if neighboring polyhedra is in todo list, skip
        i = searchsortedfirst(groebner_polyhedra_todo,
                              (w_neighbor,C_neighbor,G_neighbor),
                              by=x->x[1])
        if i<=length(groebner_polyhedra_todo) && groebner_polyhedra_todo[i][1]==w_neighbor
          continue
        end
        # otherwise, add data to todo list
        insert!(groebner_polyhedra_todo, i, (w_neighbor,C_neighbor,G_neighbor))
      end
    end
  end
end
export tropical_variety



#=======
Example:
P = cube(4)
anchor_point(P)
facet_points(P)
=======#
function anchor_point(P::Polyhedron) # todo: simplify this function
  # compute the sum of vertices and rays in homogenized coordinates
  pt = convert(Vector{fmpq},sum([vertices(P)...,rays(P)...]))
  pushfirst!(pt,nvertices(P))

  # project to orthogonal complement of lineality space if necessary
  if lineality_dim(P)>0
    pt = Polymake.Matrix{Polymake.Rational}(vcat(transpose(pt)))
    Polymake.common.project_to_orthogonal_complement(pt, nf.pm_fan.LINEALITY_SPACE)
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
