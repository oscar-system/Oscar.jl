###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################

@doc Markdown.doc"""
   orbit_polytope(V::AbstractVecOrMat, G::PermGroup)

Construct the convex hull of the orbit of the point(s) in $V$ under the action of $G$.
"""
function orbit_polytope(V::AbstractMatrix, G::PermGroup)
   if size(V)[2] != degree(G)
      throw(ArgumentError("Dimension of points and group degree need to be the same."))
   end
   generators = PermGroup_to_polymake_array(G)
   pmGroup = Polymake.group.PermutationAction(GENERATORS=generators)
   pmPolytope = Polymake.polytope.orbit_polytope(homogenize(V,1), pmGroup)
   return Polyhedron(pmPolytope)
end
function orbit_polytope(V::AbstractVector, G::PermGroup)
   return orbit_polytope(Matrix(reshape(V,(1,length(V)))), G)
end

@doc Markdown.doc"""
   cube(d [, u, l])

Construct the $[-1,1]$-cube in dimension $d$. If $u$ and $l$ are given, the $[l,u]$-cube in dimension $d$ is returned.
""" cube(d) = Polyhedron(Polymake.polytope.cube(d))
cube(d, u, l) = Polyhedron(Polymake.polytope.cube(d, u, l))



"""
    newton_polytope(poly)

Compute the Newton polytope of the given polynomial `poly`.
"""
function newton_polytope(f)
    exponents = reduce(hcat, Oscar.exponent_vectors(f))'
    convex_hull(exponents)
end




"""
   intersect(P::Polyhedron, Q::Polyhedron)

   Intersect two polyhedra.
"""
function intersect(P::Polyhedron, Q::Polyhedron)
   return Polyhedron(Polymake.polytope.intersection(pm_polytope(P), pm_polytope(Q)))
end


"""
   minkowski_sum(P::Polyhedron, Q::Polyhedron)

   Minkowski sum of two polyhedra.
"""
function minkowski_sum(P::Polyhedron, Q::Polyhedron; algorithm::Symbol=:standard)
   if algorithm == :standard
      return Polyhedron(Polymake.polytope.minkowski_sum(pm_polytope(P), pm_polytope(Q)))
   elseif algorithm == :fukuda
      return Polyhedron(Polymake.polytope.minkowski_sum_fukuda(pm_polytope(P), pm_polytope(Q)))
   else
      throw(ArgumentError("Unknown minkowski sum `algorithm` argument :" * string(algorithm)))
   end
end




#TODO: documentation  + extend to different fields.

"""
   +(P::Polyhedron, Q::Polyhedron)

   Minkowski sum of two polyhedra.
"""
+(P::Polyhedron, Q::Polyhedron) = minkowski_sum(P,Q)


#TODO: extend to different fields

"""
   *(k::Int, Q::Polyhedron)

   Returns the scaled polyhedron `kQ`.
"""
*(k::Int, P::Polyhedron) = Polyhedron(Polymake.polytope.scale(pm_polytope(P),k))


"""
   *(P::Polyhedron, k::Int)

   Returns the scaled polyhedron `kP`.
"""
*(P::Polyhedron,k::Int) = k*P


"""
   +(P::Polyhedron, v::AbstractVector)

   Returns the translation `P+v` of `P` by the vector `v`.
"""
function +(P::Polyhedron,v::AbstractVector)
    if ambient_dim(P) != length(v)
        throw(ArgumentError("Translation vector not correct dimension"))
    else
        return Polyhedron(Polymake.polytope.translate(pm_polytope(P),Polymake.Vector{Polymake.Rational}(v)))
    end
end


"""
   +(v::AbstractVector,P::Polyhedron)

   Returns the translation `v+P` of `P` by the vector `v`.
"""
+(v::AbstractVector,P::Polyhedron) = P+v

@doc Markdown.doc"""

   simplex(d[,n])

Construct the simplex which is the convex hull of the standard basis vectors
along with the origin in R^$d$, optionally scaled by $n$.
"""
simplex(d::Int64,n) = Polyhedron(Polymake.polytope.simplex(d,n))
simplex(d::Int64) = Polyhedron(Polymake.polytope.simplex(d))


@doc Markdown.doc"""

   cross(d[,n])

Construct a $d$-dimensional cross polytope around origin with vertices located at $\pm e_i$ for each unit vector $e_i$ of $R^d$.
If $n$ is not given, construct the unit cross polytope around origin.
"""
cross(d::Int64,n) = Polyhedron(Polymake.polytope.cross(d,n))
cross(d::Int64) = Polyhedron(Polymake.polytope.cross(d))

@doc Markdown.doc"""

   archimedean_solid(s)

Construct an Archimedean solid with the name given by String `s` from the list
below.  The polytopes are realized with floating point numbers and thus not
exact; Vertex-facet-incidences are correct in all cases.

# Arguments
- `s::String`: the name of the desired Archimedean solid

    Possible values:

      "truncated_tetrahedron" : Truncated tetrahedron.
          Regular polytope with four triangular and four hexagonal facets.
      "cuboctahedron" : Cuboctahedron.
          Regular polytope with eight triangular and six square facets.
      "truncated_cube" : Truncated cube.
          Regular polytope with eight triangular and six octagonal facets.
      "truncated_octahedron" : Truncated Octahedron.
          Regular polytope with six square and eight hexagonal facets.
      "rhombicuboctahedron" : Rhombicuboctahedron.
          Regular polytope with eight triangular and 18 square facets.
      "truncated_cuboctahedron" : Truncated Cuboctahedron.
          Regular polytope with 12 square, eight hexagonal and six octagonal
          facets.
      "snub_cube" : Snub Cube.
          Regular polytope with 32 triangular and six square facets.
          The vertices are realized as floating point numbers.
          This is a chiral polytope.
      "icosidodecahedron" : Icosidodecahedon.
          Regular polytope with 20 triangular and 12 pentagonal facets.
      "truncated_dodecahedron" : Truncated Dodecahedron.
          Regular polytope with 20 triangular and 12 decagonal facets.
      "truncated_icosahedron" : Truncated Icosahedron.
          Regular polytope with 12 pentagonal and 20 hexagonal facets.
      "rhombicosidodecahedron" : Rhombicosidodecahedron.
          Regular polytope with 20 triangular, 30 square and 12 pentagonal
          facets.
      "truncated_icosidodecahedron" : Truncated Icosidodecahedron.
          Regular polytope with 30 square, 20 hexagonal and 12 decagonal
          facets.
      "snub_dodecahedron" : Snub Dodecahedron.
          Regular polytope with 80 triangular and 12 pentagonal facets. The
          vertices are realized as floating point numbers. This is a chiral
          polytope.
"""
archimedean_solid(s::String) = Polyhedron(Polymake.polytope.archimedean_solid(s))


@doc Markdown.doc"""

    upper_bound_theorem(d::Int, n::Int)

Returns a polyhedron which contains the combinatioral data shared by all
simplicial d-polytopes with n vertices with the maximal number of facets as
given by McMullen's Upper-Bound-Theorem. Essentially, one can read the
``H_VECTOR`` and ``F_VECTOR`` of a polytope that attains the McMullen's
upperbounds.

Arguments:
- `d::Int`: the dimension
- `n::Int`: the number of vertices
"""
upper_bound_theorem(d::Int,n::Int) = Polyhedron(Polymake.polytope.upper_bound_theorem(d,n))
