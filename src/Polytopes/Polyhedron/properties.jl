###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

#TODO: take into account lineality space

"""
    faces(as::Type{T} = Polyhedron, P::Polyhedron, face_dim::Int)

Return the faces of `P` of dimension `face_dim` as an iterator over the type of
object given by `as`.

Optional arguments for `as` include
* `Polyhedron`/`Polyhedra`.

# Examples
A `Vector` containing the six sides of the 3-dimensional cube can be obtained
via the following input:
```jldoctest
julia> F = faces(Polyhedron, cube(3), 2)
Oscar.PolyhedronFacePolyhedronIterator(A polyhedron in ambient dimension 3, 2)

julia> collect(F)
6-element Vector{Any}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
```
"""
function faces(as::Type{T}, P::Polyhedron, face_dim::Int) where T<:Union{Polyhedron, Polyhedra}
    if (face_dim < 0)
        return nothing
    end
    pfaces = Polymake.polytope.faces_of_dim(pm_polytope(P),face_dim-size(lineality_space(P),1))
    nfaces = length(pfaces)
    rfaces = Vector{Int64}()
    sizehint!(rfaces, nfaces)
    for index in 1:nfaces
        isfar = true
        for v in pfaces[index]
            #Checking that the first coordinate is zero is equivalent
            #  to being a vertex of the far face
            if !iszero(pm_polytope(P).VERTICES[1+v[1],1])
                isfar = false
                break
            end
        end
        if !isfar
            push!(rfaces, index)
        end
    end
    return PolyhedronOrConeIterator{AsTypeIdentities(as)}(P.pm_polytope.VERTICES,pfaces[rfaces], P.pm_polytope.LINEALITY_SPACE)
end
"""
    faces(P::Polyhedron, face_dim::Int)

Return the faces of `P` of dimension `face_dim` as an iterator over `Polyhedron` objects.

# Examples
A `Vector` containing the six sides of the 3-dimensional cube can be obtained via the following input:
```jldoctest
julia> F = faces(cube(3),2)
Oscar.PolyhedronFacePolyhedronIterator(A polyhedron in ambient dimension 3, 2)

julia> collect(F)
6-element Vector{Any}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
```
"""
faces(P::Polyhedron, face_dim::Int) = faces(Polyhedron, P, face_dim)

"""
    vertices(as, P)

Return an iterator over the vertices of `P` in the format defined by `as`.

Optional arguments for `as` include
* `Points`.

See also: [`vertices_as_point_matrix`](@ref).

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and
a square:
```jldoctest
julia> P = simplex(2) + cube(2);

julia> collect(vertices(Points, P))
5-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
-1 -1
 pm::Vector<pm::Rational>
2 -1
 pm::Vector<pm::Rational>
2 1
 pm::Vector<pm::Rational>
-1 2
 pm::Vector<pm::Rational>
1 2
```
"""
function vertices(as::Type{T}, P::Polyhedron) where T<:Points
    vertices = pm_polytope(P).VERTICES
    rvertices = Vector{Int64}()
    sizehint!(rvertices, size(vertices, 1))
    for index in 1:size(vertices, 1)
        if !iszero(vertices[index, 1])
            push!(rvertices, index)
        end
    end
    return PointIterator{AsTypeIdentities(as), Polymake.Rational}(vertices[rvertices, 2:end])
end

"""
    vertices(P::Polyhedron)

Return an iterator over the vertices of a polyhedron `P` as points.

See also `vertices_as_point_matrix`.

# Examples
The following code computes the vertices of the Minkowski sum of a triangle and
a square:
```jldoctest
julia> P = simplex(2) + cube(2);

julia> collect(vertices(P))
5-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
-1 -1
 pm::Vector<pm::Rational>
2 -1
 pm::Vector<pm::Rational>
2 1
 pm::Vector<pm::Rational>
-1 2
 pm::Vector<pm::Rational>
1 2
```
"""
vertices(P::Polyhedron) = vertices(Points, P)


"""
    nrays(P::Polyhedron)

Return the number of rays of `P`.

# Examples
Reflecting the input, the upper half-plane indeed has one ray.
```jldoctest
julia> UH = convex_hull([0 0],[0 1],[1 0]);

julia> nrays(UH)
1
```
"""
nrays(P::Polyhedron) = length(pm_polytope(P).FAR_FACE)

@doc Markdown.doc"""
    nvertices(P::Polyhedron)

Return the number of vertices of `P`.

# Examples
The 3-cube's number of vertices can be obtained with this input:
```jldoctest
julia> C = cube(3);

julia> nvertices(C)
8
```
"""
nvertices(P::Polyhedron) = pm_polytope(P).N_VERTICES - nrays(P)


@doc Markdown.doc"""
    rays(as::Type{T} = Points, P::Polyhedron)

Return a minimal set of generators of the cone of unbounded directions of `P`
(i.e. its rays) in the format defined by `as`.

Optional arguments for `as` include
* `Points`.

See also `rays_as_point_matrix`.

# Examples
We can verify that the positive orthant of the plane is generated by the two
rays in positive unit direction:
```jldoctest
julia> PO = convex_hull([0 0], [1 0; 0 1]);

julia> collect(rays(Points, PO))
2-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0
 pm::Vector<pm::Rational>
0 1
```
"""
function rays(as::Type{T}, P::Polyhedron) where T<:Union{Ray, Rays}
    vertices = pm_polytope(P).VERTICES
    rrays = Vector{Int64}()
    sizehint!(rrays, size(vertices, 1))
    for index in 1:size(vertices, 1)
        if iszero(vertices[index, 1])
            push!(rrays, index)
        end
    end
    return PointIterator{AsTypeIdentities(as), Polymake.Rational}(vertices[rrays, 2:end])
end

@doc Markdown.doc"""
    rays(P::Polyhedron)

Return minimal set of generators of the cone of unbounded directions of `P`
(i.e. its rays) as points.

See also `rays_as_point_matrix`.

# Examples
We can verify that the positive orthant of the plane is generated by the two
rays in positive unit direction:
```jldoctest
julia> PO = convex_hull([0 0], [1 0; 0 1]);

julia> collect(rays(PO))
2-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0
 pm::Vector<pm::Rational>
0 1
```
"""
rays(P::Polyhedron) = rays(Ray,P)

"""
    nfacets(P::Polyhedron)

Return the number of facets of `P`.

# Examples
The number of facets of the 5-dimensional cross polytope can be retrieved via
the following line:
```jldoctest
julia> nfacets(cross(5))
32
```
"""
nfacets(P::Polyhedron) = pm_polytope(P).N_FACETS

@doc Markdown.doc"""
    facets(as::Type{T} = Halfspaces, P::Polyhedron)

Return the facets of `P` in the format defined by `as`.

The allowed values for `as` are
* `Halfspaces`,
* `Polyhedron`/`Polyhedra`.

See also `facets_as_halfspace_matrix_pair`.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```jldoctest
julia> C = cube(3);

julia> collect(facets(Polyhedron, C))
6-element Vector{Polyhedron}:
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3
 A polyhedron in ambient dimension 3

julia> collect(facets(Halfspaces, C))
6-element Vector{Tuple{Polymake.Vector{Polymake.Rational}, Polymake.Rational}}:
 (pm::Vector<pm::Rational>
-1 0 0, 1)
 (pm::Vector<pm::Rational>
1 0 0, 1)
 (pm::Vector<pm::Rational>
0 -1 0, 1)
 (pm::Vector<pm::Rational>
0 1 0, 1)
 (pm::Vector<pm::Rational>
0 0 -1, 1)
 (pm::Vector<pm::Rational>
0 0 1, 1)
```
"""
facets(as::Type{T}, P::Polyhedron) where T<:Union{Halfspace, Halfspaces, Pair, Polyhedron, Polyhedra} = HalfSpaceIterator{AsTypeIdentities(as)}(decompose_hdata(pm_polytope(P).FACETS)...)

@doc Markdown.doc"""
    facets(P::Polyhedron)

Return the facets of `P` as halfspaces.

See also `facets_as_halfspace_matrix_pair`.

# Examples
We can retrieve the six facets of the 3-dimensional cube this way:
```jldoctest
julia> C = cube(3);

julia> collect(facets(C))
6-element Vector{Tuple{Polymake.Vector{Polymake.Rational}, Polymake.Rational}}:
 (pm::Vector<pm::Rational>
-1 0 0, 1)
 (pm::Vector<pm::Rational>
1 0 0, 1)
 (pm::Vector<pm::Rational>
0 -1 0, 1)
 (pm::Vector<pm::Rational>
0 1 0, 1)
 (pm::Vector<pm::Rational>
0 0 -1, 1)
 (pm::Vector<pm::Rational>
0 0 1, 1)

```
"""
facets(P::Polyhedron) = facets(Halfspace, P)

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################
"""
    volume(P::Polyhedron)

Return the (Euclidean) volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> volume(C)
4
```
"""
volume(P::Polyhedron) = (pm_polytope(P)).VOLUME

@doc Markdown.doc"""
    normalized_volume(P::Polyhedron)

Return the (normalized) volume of `P`.

# Examples
```jldoctest
julia> C = cube(2);

julia> normalized_volume(C)
8
```
"""
normalized_volume(P::Polyhedron) = factorial(dim(P))*(pm_polytope(P)).VOLUME

"""
    dim(P::Polyhedron)

Return the dimension of `P`.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V);

julia> dim(P)
2
```
"""
dim(P::Polyhedron) = Polymake.polytope.dim(pm_polytope(P))


@doc Markdown.doc"""
    lattice_points(P::Polyhedron)

Return the integer points contained in the bounded polyhedron `P`.

# Examples
```jldoctest
julia> S = 2 * simplex(2);

julia> lattice_points(S)
6-element Vector{Polymake.VectorAllocated{Polymake.Integer}}:
 pm::Vector<pm::Integer>
0 0
 pm::Vector<pm::Integer>
0 1
 pm::Vector<pm::Integer>
0 2
 pm::Vector<pm::Integer>
1 0
 pm::Vector<pm::Integer>
1 1
 pm::Vector<pm::Integer>
2 0
```
"""
function lattice_points(P::Polyhedron)
    if pm_polytope(P).BOUNDED
        lat_pts = pm_polytope(P).LATTICE_POINTS_GENERATORS[1]
        return PointIterator{Points, Polymake.Integer}(lat_pts[:, 2:end])
    else
        throw(ArgumentError("Polyhedron not bounded"))
    end
end

@doc Markdown.doc"""
    ambient_dim(P::Polyhedron)

Return the ambient dimension of `P`.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V);

julia> ambient_dim(P)
3
```
"""
ambient_dim(P::Polyhedron) = Polymake.polytope.ambient_dim(pm_polytope(P))

"""
    codim(P::Polyhedron)

Return the codimension of `P`.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> P = convex_hull(V);

julia> codim(P)
1
```
"""
codim(P::Polyhedron) = ambient_dim(P)-dim(P)

###############################################################################
## Points properties
###############################################################################


# Previously: This implementation is not correct. Ask Taylor.
# Taylor: lineality space generators always look like [0, v] so
#  v is a natural output.
@doc Markdown.doc"""
    lineality_space(P::Polyhedron)

Return a matrix whose row span is the lineality space of `P`.

# Examples
Despite not being reflected in this construction of the upper half-plane,
its lineality in $x$-direction is recognized:
```jldoctest
julia> UH = convex_hull([0 0],[0 1; 1 0; -1 0]);

julia> lineality_space(UH)
pm::Matrix<pm::Rational>
1 0
```
"""
lineality_space(P::Polyhedron) = PointIterator{Ray, Polymake.Rational}(dehomogenize(P.pm_polytope.LINEALITY_SPACE))


@doc Markdown.doc"""
    recession_cone(P::Polyhedron)

Return the recession cone of `P`.

# Examples
```jldoctest
julia> P = Polyhedron([1 -2; -1 1; -1 0; 0 -1],[2,1,1,1]);

julia> collect(vertices(P))
3-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
0 -1
 pm::Vector<pm::Rational>
-1 0
 pm::Vector<pm::Rational>
-1 -1

julia> recession_cone(P)
A polyhedral cone in ambient dimension 2

julia> collect(rays(recession_cone(P)))
2-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 1/2
 pm::Vector<pm::Rational>
1 1
```
"""
recession_cone(P::Polyhedron) = Cone(Polymake.polytope.recession_cone(pm_polytope(P)))

###############################################################################
## Boolean properties
###############################################################################
"""
    isfeasible(P::Polyhedron)

Check whether `P` is feasible, i.e. non-empty.

# Examples
```jldoctest
julia> P = Polyhedron([1 -1; -1 1; -1 0; 0 -1],[-1,-1,1,1]);

julia> isfeasible(P)
false
```
"""
isfeasible(P::Polyhedron) = pm_polytope(P).FEASIBLE

"""
    contains(P::Polyhedron, v::AbstractVector)

Check whether `P` contains `v`.

# Examples
The positive orthant only contains vectors with non-negative entries:
```jldoctest
julia> PO = Polyhedron([-1 0; 0 -1], [0, 0]);

julia> contains(PO, [1, 2])
true

julia> contains(PO, [1, -2])
false
```
"""
contains(P::Polyhedron, v::AbstractVector) = Polymake.polytope.contains(pm_polytope(P), [1; v])

"""
    issmooth(P::Polyhedron)

Check whether `P` is smooth.

# Examples
A cube is always smooth.
```jldoctest
julia> C = cube(8);

julia> issmooth(C)
true
```
"""
issmooth(P::Polyhedron) = pm_polytope(P).SMOOTH


"""
    isnormal(P::Polyhedron)

Check whether `P` is normal.

# Examples
The 3-cube is normal.
```jldoctest
julia> C = cube(3)
A polyhedron in ambient dimension 3

julia> isnormal(C)
true
```
But this pyramid is not:
```jldoctest
julia> P = convex_hull([0 0 0; 0 1 1; 1 1 0; 1 0 1]);

julia> isnormal(P)
false
```
"""
isnormal(P::Polyhedron) = pm_polytope(P).NORMAL


"""
    isbounded(P::Polyhedron)

Check whether `P` is bounded.

# Examples
```jldoctest
julia> P = Polyhedron([1 -3; -1 1; -1 0; 0 -1],[1,1,1,1]);

julia> isbounded(P)
false
```
"""
isbounded(P::Polyhedron) = pm_polytope(P).BOUNDED


"""
    isfulldimensional(P::Polyhedron)

Check whether `P` is full-dimensional.

# Examples
```jldoctest
julia> V = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

julia> isfulldimensional(convex_hull(V))
false
```
"""
isfulldimensional(P::Polyhedron) = pm_polytope(P).FULL_DIM

@doc Markdown.doc"""
    f_vector(P::Polyhedron)

Compute the vector $(f₁,f₂,...,f_{(dim(P)-1))$` where $f_i$ is the number of
faces of $P$ of dimension $i$.

# Examples
Here we compute the f-vector of the 5-cube:
```jldoctest
julia> f_vector(cube(5))
5-element Vector{Int64}:
 32
 80
 80
 40
 10
```
"""
function f_vector(P::Polyhedron)
    f_vec=[length(collect(faces(P,i))) for i in 0:dim(P)-1]
    return f_vec
end


@doc Markdown.doc"""
    support_function(P::Polyhedron; convention::Symbol = :max)

Produce a function $h(ω) = max\{dot(x,ω)\ |\ x \in P\}$. $max$ may be changed
to $min$ by setting `convention = :min`.

# Examples
```jldoctest
julia> P = cube(3) + simplex(3);

julia> φ = support_function(P);

julia> φ([1,2,3])
9

julia> ψ = support_function(P, convention = :min);

julia> ψ([1,2,3])
-6
```
"""
function support_function(P::Polyhedron; convention = :max)
    function h(ω::AbstractVector)
        lp=LinearProgram(P,ω; convention = convention)
        return solve_lp(lp)[1]
    end
    return h
end

@doc Markdown.doc"""
    print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

Trivial inequalities are counted but omitted. They are included if `trivial` is
set to `true`.

# Examples
The 3-cube is given by $-1 ≦ x_i ≦ 0 ∀ i ∈ \{1, 2, 3\}$.
```jldoctest
julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5])
1: -x₁ + 4*x₃ + 5*x₄ ≦ 0
2: 4*x₁ + 4*x₂ + 4*x₃ + 3*x₄ ≦ 1
3: x₁ ≦ 2
5: 0 ≦ -4
6: 9*x₁ + 9*x₂ + 9*x₃ + 9*x₄ ≦ 5

julia> print_constraints([-1 0 4 5; 4 4 4 3; 1 0 0 0; 0 0 0 0; 0 0 0 0; 9 9 9 9], [0, 1, 2, 3, -4, 5]; trivial = true)
1: -x₁ + 4*x₃ + 5*x₄ ≦ 0
2: 4*x₁ + 4*x₂ + 4*x₃ + 3*x₄ ≦ 1
3: x₁ ≦ 2
4: 0 ≦ 3
5: 0 ≦ -4
6: 9*x₁ + 9*x₂ + 9*x₃ + 9*x₄ ≦ 5
```
"""
function print_constraints(A::AnyVecOrMat, b::AbstractVector; trivial::Bool = false)
    for i in 1:length(b)
        terms = Vector{String}(undef, size(A)[2])
        first = true
        for j in 1:size(A)[2]
            if iszero(A[i, j])
                terms[j] = ""
            else
                if isone(abs(A[i, j]))
                    terms[j] = first ? string(isone(A[i, j]) ? "x" : "-x" , ['₀'+ d for d in digits(j)]...) :
                        string(isone(A[i, j]) ? " + x" : " - x", ['₀'+ d for d in digits(j)]...)
                else
                    terms[j] = first ? string(A[i, j], "*x", ['₀'+ d for d in digits(j)]...) :
                        string(A[i, j] < 0 ? " - " : " + ", abs(A[i, j]), "*x", ['₀'+ d for d in digits(j)]...)
                end
                first = false
            end
        end
        if first
            if b[i] >= 0 && !trivial
                continue
            end
            terms[1] = "0"
        end
        println(string(i, ": ", terms..., " ≦ ", b[i]))
    end
end

@doc Markdown.doc"""
    print_constraints(P::Polyhedron; trivial::Bool = false)

Pretty print the constraints given by $P(A,b) = \{ x |  Ax ≤ b \}$.

Trivial inequalities are counted but omitted. They are included if `trivial` is
set to `true`.

# Examples
The 3-cube is given by $-1 ≦ x_i ≦ 0 ∀ i ∈ \{1, 2, 3\}$.
```jldoctest
julia> print_constraints(cube(3))
1: -x₁ ≦ 1
2: x₁ ≦ 1
3: -x₂ ≦ 1
4: x₂ ≦ 1
5: -x₃ ≦ 1
6: x₃ ≦ 1
```
"""
print_constraints(P::Polyhedron; trivial::Bool = false) = print_constraints(facets_as_halfspace_matrix_pair(P)...; trivial = trivial)
