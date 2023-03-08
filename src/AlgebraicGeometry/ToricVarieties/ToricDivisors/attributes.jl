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
```jldoctest
julia> F4 = hirzebruch_surface(4)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> td0 = toric_divisor(F4, [0,0,0,0])
Torus-invariant, non-prime divisor on a normal toric variety

julia> is_feasible(polyhedron(td0))
true

julia> dim(polyhedron(td0))
0

julia> td1 = toric_divisor(F4, [1,0,0,0])
Torus-invariant, prime divisor on a normal toric variety

julia> is_feasible(polyhedron(td1))
true

julia> td2 = toric_divisor(F4, [-1,0,0,0])
Torus-invariant, non-prime divisor on a normal toric variety

julia> is_feasible(polyhedron(td2))
false
```
"""
@attr Polyhedron polyhedron(td::ToricDivisor) = Polyhedron(pm_tdivisor(td).SECTION_POLYTOPE)


@doc Markdown.doc"""
    coefficients(td::ToricDivisor)

Identify the coefficients of a toric divisor in the group of torus invariant Weil divisors.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> D = toric_divisor(F4, [1, 2, 3, 4])
Torus-invariant, non-prime divisor on a normal toric variety

julia> coefficients(D)
4-element Vector{ZZRingElem}:
 1
 2
 3
 4
```
"""
function coefficients(td::ToricDivisor)
    return td.coeffs
end


@doc Markdown.doc"""
    toric_variety(td::ToricDivisor)

Return the toric variety of a torus-invariant Weil divisor.

# Examples
```jldoctest
julia> F4 = hirzebruch_surface(4);

julia> D = toric_divisor(F4, [1, 2, 3, 4]);

julia> toric_variety(D)
Normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor
```
"""
function toric_variety(td::ToricDivisor)
    return td.toric_variety
end
