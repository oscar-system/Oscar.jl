@doc Markdown.doc"""
    polyhedron(td::ToricDivisor)

Construct the polyhedron $P_D$ of a torus invariant divisor $D:=td$ as in 4.3.2
of [CLS11](@cite). The lattice points of this polyhedron correspond to the
global sections of the divisor.

# Examples
The polyhedron of the divisor with all coefficients equal to zero is a point,
if the ambient variety is complete. Changing the coefficients corresponds to
moving hyperplanes. One direction moves the hyperplane away from the origin,
the other moves it across. In the latter case there are no global sections
anymore and the polyhedron becomes empty.
```
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td0 = ToricDivisor(H, [0,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isfeasible(polyhedron(td0))
true

julia> dim(polyhedron(td0))
0

julia> td1 = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isfeasible(polyhedron(td1))
true

julia> td2 = ToricDivisor(H, [-1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isfeasible(polyhedron(td2))
false
```
"""
function polyhedron(td::ToricDivisor)
    if !has_attribute(td, :polyhedron)
        set_attribute!(td, :polyhedron, Polyhedron(pm_tdivisor(td).SECTION_POLYTOPE))
    end
    return get_attribute(td, :polyhedron)
end
export polyhedron


@doc Markdown.doc"""
    coefficients(td::ToricDivisor)

Identify the coefficients of a toric divisor in the group of torus invariant Weil divisors.

# Examples
```
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> D = ToricDivisor(H, [1,2,3,4])
A torus invariant divisor on a normal toric variety

julia> coefficients(D)
4-element Vector{Int64}:
 1
 2
 3
 4
```
"""
function coefficients(td::ToricDivisor)
    if !has_attribute(td, :coefficients)
        set_attribute!(td, :coefficients, Vector{Int}(Polymake.common.primitive(pm_tdivisor(td).COEFFICIENTS)))
    end
    return get_attribute(td, :coefficients)
end
export coefficients
