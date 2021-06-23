###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################

@doc Markdown.doc"""
   orbit_polytope(V, G)

Construct the convex hull of the orbit of one or several points under the action of a permutation group.

# Arguments
- `V::AbstractVecOrMat`: Initial point(s).
- `P::PermGroup`: A permutation group.

# Examples
This will construct the $3$-dimensional permutahedron:
```julia-repl
julia> V = [1 2 3];

julia> G = symmetric_group(3);

julia> P = orbit_polytope(V, G)
A polyhedron in ambient dimension 3

julia> collect(vertices(P))
6-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 2 3
 pm::Vector<pm::Rational>
1 3 2
 pm::Vector<pm::Rational>
2 1 3
 pm::Vector<pm::Rational>
2 3 1
 pm::Vector<pm::Rational>
3 1 2
 pm::Vector<pm::Rational>
3 2 1
```
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
   cube(d [, l, u])

Construct the $[-1,1]$-cube in dimension $d$. If $l$ and $u$ are given, the $[l,u]$-cube in dimension $d$ is returned.

# Arguments
- `d::Int`: Dimension of the cube.
- `l::Rational`: Lower bound for each coordinate.
- `u::Rational`: Upper bound for each coordinate.

# Examples
In this example the 5-dimensional unit cube is constructed to ask for one of its properties:
```julia-repl
julia> C = cube(5,0,1);

julia> normalized_volume(C)
120
```
"""
cube(d) = Polyhedron(Polymake.polytope.cube(d))
cube(d, l, u) = Polyhedron(Polymake.polytope.cube(d, u, l))



"""
    newton_polytope(poly)

Compute the Newton polytope of the given polynomial `poly`.

# Arguments
- `poly::Polynomial`: A multivariate polynomial.

# Examples
```julia-repl
julia> S, (x, y) = PolynomialRing(ZZ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Integer Ring, fmpz_mpoly[x, y])

julia> f = x^3*y + 3x*y^2 + 1
x^3*y + 3*x*y^2 + 1

julia> NP = newton_polytope(f)
A polyhedron in ambient dimension 2

julia> collect(vertices(NP))
3-element Array{Polymake.Vector{Polymake.Rational},1}:
 pm::Vector<pm::Rational>
3 1
 pm::Vector<pm::Rational>
1 2
 pm::Vector<pm::Rational>
0 0
```
"""
function newton_polytope(f)
    exponents = reduce(hcat, Oscar.exponent_vectors(f))'
    convex_hull(exponents)
end




@doc Markdown.doc"""
   intersect(P, Q)

Intersect two polyhedra.

# Arguments
- `P::Polyhedron`: First polyhedron.
- `Q::Polyhedron`: Second polyhedron.

# Examples
The positive orthant of the plane is the intersection of the two halfspaces with $x>0$ and $y>0$ respectively.
```julia-repl
julia> UH1 = convex_hull([0 0],[1 0],[0 1]);

julia> UH2 = convex_hull([0 0],[0 1],[1 0]);

julia> PO = intersect(UH1, UH2)
A polyhedron in ambient dimension 2

julia> collect(rays(PO))
2-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0
 pm::Vector<pm::Rational>
0 1
```
"""
function intersect(P::Polyhedron, Q::Polyhedron)
   return Polyhedron(Polymake.polytope.intersection(pm_polytope(P), pm_polytope(Q)))
end


"""
   minkowski_sum(P::Polyhedron, Q::Polyhedron)

Minkowski sum of two polyhedra.

# Arguments
- `P::Polyhedron`: First polyhedron.
- `Q::Polyhedron`: Second polyhedron.

# Examples
The Minkowski sum of a square and the 2-dimensional cross-polytope is an octagon:
```julia-repl
julia> P = cube(2);

julia> Q = cross(2);

julia> M = minkowski_sum(P, Q)
A polyhedron in ambient dimension 2

julia> nvertices(M)
8
```
"""
function minkowski_sum(P::Polyhedron, Q::Polyhedron; algorithm::Symbol=:standard)
   if algorithm == :standard
      return Polyhedron(Polymake.polytope.minkowski_sum(pm_polytope(P), pm_polytope(Q)))
   elseif algorithm == :fukuda
      return Polyhedron(Polymake.polytope.minkowski_sum_fukuda(pm_polytope(P), pm_polytope(Q)))
   else
      throw(ArgumentError("Unknown minkowski sum `algorithm` argument:" * string(algorithm)))
   end
end




#TODO: documentation  + extend to different fields.

"""
   +(P::Polyhedron, Q::Polyhedron)

Minkowski sum of two polyhedra.

# Arguments
- `P::Polyhedron`: First polyhedron.
- `Q::Polyhedron`: Second polyhedron.

# Examples
The Minkowski sum of a square and the 2-dimensional cross-polytope is an octagon:
```julia-repl
julia> P = cube(2);

julia> Q = cross(2);

julia> M = minkowski_sum(P, Q)
A polyhedron in ambient dimension 2

julia> nvertices(M)
8
```
"""
+(P::Polyhedron, Q::Polyhedron) = minkowski_sum(P,Q)


#TODO: extend to different fields

@doc Markdown.doc"""
   *(k::Int, Q::Polyhedron)

Return the scaled polyhedron `kQ`.

# Arguments
- `k::Int`: Scaling factor.
- `Q::Polyhedron`: A polyhedron.

# Examples
Scaling an $n$-dimensional bounded polyhedron by the factor $k$ results in the volume being scaled by $k^n$.
This example confirms the statement for the 6-dimensional cube and $k = 2$.
```julia-repl
julia> C = cube(6);

julia> SC = 2*C
A polyhedron in ambient dimension 6

julia> volume(SC)//volume(C)
64
```
"""
*(k::Int, P::Polyhedron) = Polyhedron(Polymake.polytope.scale(pm_polytope(P),k))


@doc Markdown.doc"""
   *(P::Polyhedron, k::Int)

Return the scaled polyhedron `kP`.

# Arguments
- `k::Int`: Scaling factor.
- `Q::Polyhedron`: A polyhedron.

# Examples
Scaling an $n$-dimensional bounded polyhedron by the factor $k$ results in the volume being scaled by $k^n$.
This example confirms the statement for the 6-dimensional cube and $k = 2$.
```julia-repl
julia> C = cube(6);

julia> SC = C*2
A polyhedron in ambient dimension 6

julia> volume(SC)//volume(C)
64
```
"""
*(P::Polyhedron,k::Int) = k*P


@doc Markdown.doc"""
   +(P::Polyhedron, v::AbstractVector)

Return the translation `P+v` of `P` by the vector `v`.

# Arguments
- `P::Polyhedron`: A polyhedron.
- `v::AbstractVector`: A vector of the same dimension as the ambient space of `P`.

# Examples
We construct a polyhedron from its $V$-description. Shifting it by the right vector reveals that its inner geometry
corresponds to that of the 3-simplex.
```julia-repl
julia> P = convex_hull([100 200 300; 101 200 300; 100 201 300; 100 200 301]);

julia> v = [-100, -200, -300];

julia> S = P + v
A polyhedron in ambient dimension 3

julia> collect(vertices(S))
4-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
0 0 0
 pm::Vector<pm::Rational>
1 0 0
 pm::Vector<pm::Rational>
0 1 0
 pm::Vector<pm::Rational>
0 0 1
```
"""
function +(P::Polyhedron,v::AbstractVector)
    if ambient_dim(P) != length(v)
        throw(ArgumentError("Translation vector not correct dimension"))
    else
        return Polyhedron(Polymake.polytope.translate(pm_polytope(P),Polymake.Vector{Polymake.Rational}(v)))
    end
end


@doc Markdown.doc"""
   +(v::AbstractVector,P::Polyhedron)

Return the translation `v+P` of `P` by the vector `v`.

# Arguments
- `P::Polyhedron`: A polyhedron.
- `v::AbstractVector`: A vector of the same dimension as the ambient space of `P`.

# Examples
We construct a polyhedron from its $V$-description. Shifting it by the right vector reveals that its inner geometry
corresponds to that of the 3-simplex.
```julia-repl
julia> P = convex_hull([100 200 300; 101 200 300; 100 201 300; 100 200 301]);

julia> v = [-100, -200, -300];

julia> S = v + P
A polyhedron in ambient dimension 3

julia> collect(vertices(S))
4-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
0 0 0
 pm::Vector<pm::Rational>
1 0 0
 pm::Vector<pm::Rational>
0 1 0
 pm::Vector<pm::Rational>
0 0 1
```
"""
+(v::AbstractVector,P::Polyhedron) = P+v

@doc Markdown.doc"""

   simplex(d[,n])

Construct the simplex which is the convex hull of the standard basis vectors
along with the origin in $\mathbb{R}^d$, optionally scaled by $n$.

# Arguments
- `d::Int`: Dimension of the simplex (and its ambient space).
- `n::Scalar`: Scaling factor.

# Examples
Here we take a look at the facets of the 7-simplex and a scaled 7-simplex:
```julia-repl
julia> s = simplex(7)
A polyhedron in ambient dimension 7

julia> collect(facets(s))
8-element Vector{Tuple{Polymake.Vector{Polymake.Rational}, Polymake.Rational}}:
 (pm::Vector<pm::Rational>
-1 0 0 0 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 -1 0 0 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 -1 0 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 -1 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 0 -1 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 0 0 -1 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 0 0 0 -1, 0)
 (pm::Vector<pm::Rational>
1 1 1 1 1 1 1, 1)

julia> t = simplex(7, 5)
A polyhedron in ambient dimension 7

julia> collect(facets(t))
8-element Vector{Tuple{Polymake.Vector{Polymake.Rational}, Polymake.Rational}}:
 (pm::Vector<pm::Rational>
-1 0 0 0 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 -1 0 0 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 -1 0 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 -1 0 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 0 -1 0 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 0 0 -1 0, 0)
 (pm::Vector<pm::Rational>
0 0 0 0 0 0 -1, 0)
 (pm::Vector<pm::Rational>
1 1 1 1 1 1 1, 5)
```
"""
simplex(d::Int64,n) = Polyhedron(Polymake.polytope.simplex(d,n))
simplex(d::Int64) = Polyhedron(Polymake.polytope.simplex(d))


@doc Markdown.doc"""

   cross(d[,n])

Construct a $d$-dimensional cross polytope around origin with vertices located at $\pm e_i$ for each unit vector $e_i$ of $R^d$, scaled by $n$.

# Arguments
- `d::Int`: Dimension of the cross polytope (and its ambient space).
- `n::Scalar`: Scaling factor.

# Examples
Here we print the facets of a non-scaled and a scaled 3-dimensional cross polytope:
```julia-repl
julia> C = cross(3)
A polyhedron in ambient dimension 3

julia> collect(facets(C))
8-element Vector{Tuple{Polymake.Vector{Polymake.Rational}, Polymake.Rational}}:
 (pm::Vector<pm::Rational>
1 1 1, 1)
 (pm::Vector<pm::Rational>
-1 1 1, 1)
 (pm::Vector<pm::Rational>
1 -1 1, 1)
 (pm::Vector<pm::Rational>
-1 -1 1, 1)
 (pm::Vector<pm::Rational>
1 1 -1, 1)
 (pm::Vector<pm::Rational>
-1 1 -1, 1)
 (pm::Vector<pm::Rational>
1 -1 -1, 1)
 (pm::Vector<pm::Rational>
-1 -1 -1, 1)

julia> D = cross(3, 2)
A polyhedron in ambient dimension 3

julia> collect(facets(D))
8-element Vector{Tuple{Polymake.Vector{Polymake.Rational}, Polymake.Rational}}:
 (pm::Vector<pm::Rational>
1 1 1, 2)
 (pm::Vector<pm::Rational>
-1 1 1, 2)
 (pm::Vector<pm::Rational>
1 -1 1, 2)
 (pm::Vector<pm::Rational>
-1 -1 1, 2)
 (pm::Vector<pm::Rational>
1 1 -1, 2)
 (pm::Vector<pm::Rational>
-1 1 -1, 2)
 (pm::Vector<pm::Rational>
1 -1 -1, 2)
 (pm::Vector<pm::Rational>
-1 -1 -1, 2)
```
"""
cross(d::Int64,n) = Polyhedron(Polymake.polytope.cross(d,n))
cross(d::Int64) = Polyhedron(Polymake.polytope.cross(d))

@doc Markdown.doc"""

   archimedean_solid(s)

Construct an Archimedean solid with the name given by String `s` from the list below.
The polytopes are realized with floating point numbers and thus not exact; Vertex-facet-incidences are correct in all cases.

# Arguments
- `s::String`: The name of the desired Archimedean solid.

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
          Regular polytope with 12 square, eight hexagonal and six octagonal facets.
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
          Regular polytope with 20 triangular, 30 square and 12 pentagonal facets.
      "truncated_icosidodecahedron" : Truncated Icosidodecahedron.
          Regular polytope with 30 square, 20 hexagonal and 12 decagonal facets.
      "snub_dodecahedron" : Snub Dodecahedron.
          Regular polytope with 80 triangular and 12 pentagonal facets.
          The vertices are realized as floating point numbers.
          This is a chiral polytope.

# Examples
```julia-repl
julia> T = archimedean_solid("cuboctahedron")
A polyhedron in ambient dimension 3

julia> sum([nvertices(F) for F in faces(T, 2)] .== 3)
8

julia> sum([nvertices(F) for F in faces(T, 2)] .== 4)
6

julia> nfacets(T)
14
```
"""
archimedean_solid(s::String) = Polyhedron(Polymake.polytope.archimedean_solid(s))
