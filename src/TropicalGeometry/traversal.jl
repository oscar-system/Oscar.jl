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
Random.seed!(13371337)
K,s = RationalFunctionField(QQ,"s")
Kx,(x1,x2,x3,x4) = PolynomialRing(K,4)
val = ValuationMap(K,s)
I = ideal([x1-s*x2+(s+1)*x3,3*x2-s^2*x3+(s^2+1)*x4])
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
    w = convert(Vector{fmpq},sum(vertices(C)))
    push!(groebner_polyhedra_todo, (w,C,G))
  end

  while !isempty(groebner_polyhedra_todo)
    println("length(groebner_polyhedra_todo): ",length(groebner_polyhedra_todo))
    println("length(groebner_polyhedra_done): ",length(groebner_polyhedra_done))

    (w,C,G) = popfirst!(groebner_polyhedra_todo)

    facet_points_to_traverse = facet_points(C)
    for facet_point in facet_points_to_traverse

      # skip facet_point if it lies in links_done
      index = searchsortedfirst(links_done,facet_point)
      if links_done[index]==facet_point # todo: is facet_point the sum of the vertices?
        continue
      end
      # add facet_point to links_done
      insert!(links_done, index, facet_point)

      # compute tropical link
      facet_link = tropical_link(ideal(G),val,facet_point)
      for u in facet_link
        # compute neighboring data
        G_neighbor = groebner_flip(ideal(G),val,facet_point,u)
        C_neighbor = groebner_polyhedron(I,val,starting_point) # todo: this computes an unnecessary GB
        w_neighbor = convert(Vector{fmpq},sum(vertices(C)))

        # if neighboring polyhedra is in done list, skip
        index = searchsortedfirst(groebner_polyhedra_done,w_neighbor,by=x->x[1]t)
        if groebner_polyhedra_done[index][1]==w_neighbor
          continue
        end

        # if neighboring polyhedra is in todo list, skip
        index = searchsortedfirst(groebner_polyhedra_done,w_neighbor,by=x->x[1]t)
        if groebner_polyhedra_done[index][1]==w_neighbor
          continue
        end
        # otherwise, add data to todo list
        insert!(groebner_polyhedra_todo, index, (w_neighbor,C_neighbor,G_neighbor))
      end
    end
  end
end
export tropical_variety



#=======
tropical variety of an ideal
todo: proper documentation
Example:
P = cube(4)
facet_points(P)
=======#
function facet_points(P::Polyhedron)
  points = []
  for facet in faces(P,dim(P)-1)
    push!(points,convert(Vector{fmpq},relative_interior_point(facet)))
  end
  return points
end
export facet_points
