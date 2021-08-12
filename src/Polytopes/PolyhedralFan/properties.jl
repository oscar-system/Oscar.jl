###############################################################################
###############################################################################
### Iterators
###############################################################################
###############################################################################

"""
    rays(PF::PolyhedralFan)

Return the rays of `PF`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C);

julia> collect(rays(NF))
6-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0 0
 pm::Vector<pm::Rational>
-1 0 0
 pm::Vector<pm::Rational>
0 1 0
 pm::Vector<pm::Rational>
0 -1 0
 pm::Vector<pm::Rational>
0 0 1
 pm::Vector<pm::Rational>
0 0 -1
```
"""
rays(PF::PolyhedralFan) = PointIterator{Polymake.Vector, Int64}(pm_fan(PF).RAYS)

"""
    maximal_cones(PF::PolyhedralFan)

Return the maximal cones of `PF`.

# Examples
Here we ask for the the number of rays for each maximal cone of the face fan of
the 3-cube and use that `maximal_cones` returns an iterator.
```jldoctest
julia> PF = face_fan(cube(3));

julia> for c in maximal_cones(PF)
       println(nrays(c))
       end
4
4
4
4
4
4
```
"""

#TODO: should the documentation mention maximal_cones_as_incidence_matrix?
#      similarly for cone ray iterators and facet iterators?
@doc Markdown.doc"""
    maximal_cones(PF::PolyhedralFan, as = :cones)

Return an iterator over the maximal cones of `PF`.
"""
function maximal_cones(PF::PolyhedralFan)
   MaximalConeIterator(PF)
end

struct MaximalConeIterator
    PF::PolyhedralFan
end

function Base.iterate(iter::MaximalConeIterator, index = 1)
    if index > nmaximal_cones(iter.PF)
        return nothing
    end
    current_cone = Cone(Polymake.fan.cone(pm_fan(iter.PF), index - 1))
    return (current_cone, index + 1)
end
Base.length(iter::MaximalConeIterator) = nmaximal_cones(iter.PF)

# function cones(as::Type{T}, PF::PolyhedralFan, cone_dim::Int) where {T}
#     rtype = AsTypeIdentitiesC(as)
#     if (cone_dim < 0)
#         return nothing
#     end
#     rcones = Polymake.fan.cones_of_dim(PF.pm_fan,cone_dim-length(lineality_space(PF)))
#     return PolyhedronOrConeIterator{rtype}(PF.pm_fan.RAYS,rcones, PF.pm_fan.LINEALITY_SPACE)
# end

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################

###############################################################################
## Scalar properties
###############################################################################

"""
    dim(PF::PolyhedralFan)

Return the dimension of `PF`.

# Examples
This fan in the plane contains a 2-dimensional cone and is thus 2-dimensional
itself.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));

julia> dim(PF)
2
```
"""
dim(PF::PolyhedralFan) = pm_fan(PF).FAN_DIM

"""
    nmaximal_cones(PF::PolyhedralFan)

Return the number of maximal cones of `PF`.

# Examples
The cones given in this construction are non-redundant. Thus there are two
maximal cones.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));

julia> nmaximal_cones(PF)
2
```
"""
nmaximal_cones(PF::PolyhedralFan) = pm_fan(PF).N_MAXIMAL_CONES

"""
    ambient_dim(PF::PolyhedralFan)

Return the ambient dimension `PF`, which is the dimension of the embedding
space.

This is equal to the dimension of the fan if and only if the fan is
full-dimensional.

# Examples
The normal fan of the 4-cube is embedded in the same ambient space.
```jldoctest
julia> ambient_dim(normal_fan(cube(4)))
4
```
"""
ambient_dim(PF::PolyhedralFan) = pm_fan(PF).FAN_AMBIENT_DIM

"""
    nrays(PF::PolyhedralFan)

Return the number of rays of `PF`.

# Examples
The 3-cube has 8 vertices. Accordingly, its face fan has 8 rays.
```jldoctest
julia> nrays(face_fan(cube(3)))
8
```
"""
nrays(PF::PolyhedralFan) = pm_fan(PF).N_RAYS


###############################################################################
## Points properties
###############################################################################

@doc Markdown.doc"""
    lineality_space(PF::PolyhedralFan)

Return a non-redundant matrix whose rows are generators of the lineality space
of `PF`.

# Examples
This fan consists of two cones, one containing all the points with $y ≤ 0$ and
one containing all the points with $y ≥ 0$. The fan's lineality is the common
lineality of these two cones, i.e. in $x$-direction.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 0; 0 -1], IncidenceMatrix([[1, 2, 3], [3, 4, 1]]))
A polyhedral fan in ambient dimension 2

julia> lineality_space(PF)
pm::Matrix<pm::Rational>
1 0

```
"""
lineality_space(PF::PolyhedralFan) = PointIterator{Polymake.Vector, Polymake.Rational}(pm_fan(PF).LINEALITY_SPACE)


"""
    rays_as_point_matrix(PF::PolyhedralFan)

Return the rays of `PF` as rows of a matrix.

# Examples
The rays of a face fan (centered at the origin) point into the direction of the
polytope's vertices; here we see this for the 3-cube.
```jldoctest
julia> collect(rays(face_fan(cube(3))))
8-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
-1 -1 -1
 pm::Vector<pm::Rational>
1 -1 -1
 pm::Vector<pm::Rational>
-1 1 -1
 pm::Vector<pm::Rational>
1 1 -1
 pm::Vector<pm::Rational>
-1 -1 1
 pm::Vector<pm::Rational>
1 -1 1
 pm::Vector<pm::Rational>
-1 1 1
 pm::Vector<pm::Rational>
1 1 1

```
"""
rays_as_point_matrix(PF::PolyhedralFan) = pm_fan(PF).RAYS


@doc Markdown.doc"""
    maximal_cones_as_incidence_matrix(PF::PolyhedralFan)

Return the maximal cones of `PF` as an incidence matrix.

The rows of the output correspond to the maximal cones and the columns
correspond to the rays.

# Examples
From this output we can read that the normal fan of the 2-dimensional cube
consists of four cones, generated by $\{ \pm e_1, \pm e_2 \}$, where $e_1$ and
$e_2$ are the two unit vectors. These cones correspond to the four quadrants.
```jldoctest
julia> PF = normal_fan(cube(2));

julia> collect(rays(PF))
4-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0
 pm::Vector<pm::Rational>
-1 0
 pm::Vector<pm::Rational>
0 1
 pm::Vector<pm::Rational>
0 -1

julia> maximal_cones_as_incidence_matrix(PF)
4×4 Matrix{Bool}:
 1  0  1  0
 0  1  1  0
 1  0  0  1
 0  1  0  1
```
"""
function maximal_cones_as_incidence_matrix(PF::PolyhedralFan)
   IncidenceMatrix(pm_fan(PF).MAXIMAL_CONES)
end

###############################################################################
## Boolean properties
###############################################################################
"""
    issmooth(PF::PolyhedralFan)

Determine whether `PF` is smooth.

# Examples
Even though the cones of this fan cover the positive orthant together, one of
these und thus the whole fan is not smooth.
```jldoctest
julia> PF = PolyhedralFan([0 1; 2 1; 1 0], IncidenceMatrix([[1, 2], [2, 3]]));

julia> issmooth(PF)
false
```
"""
issmooth(PF::PolyhedralFan) = pm_fan(PF).SMOOTH_FAN

"""
    isregular(PF::PolyhedralFan)

Determine whether `PF` is regular, i.e. the normal fan of a polytope.

# Examples
This fan is not complete and thus not regular.
```jldoctest
julia> PF = PolyhedralFan([1 0; 0 1; -1 -1], IncidenceMatrix([[1, 2], [3]]));

julia> isregular(PF)
false
```
"""
isregular(PF::PolyhedralFan) = pm_fan(PF).REGULAR

"""
    iscomplete(PF::PolyhedralFan)

Determine whether `PF` is complete, i.e. its support, the set-theoretic union
of its cones, covers the whole space.

# Examples
Normal fans of polytopes are complete.
```jldoctest
julia> iscomplete(normal_fan(cube(3)))
true
```
"""
iscomplete(PF::PolyhedralFan) = pm_fan(PF).COMPLETE
