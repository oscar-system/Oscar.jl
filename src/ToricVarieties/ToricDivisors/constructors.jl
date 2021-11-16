######################
# 1: The Julia type for ToricDivisors
######################

struct ToricDivisor
           polymake_divisor::Polymake.BigObject
end
export ToricDivisor


function pm_tdivisor(td::ToricDivisor)
    return td.polymake_divisor
end

######################
# 2: Generic constructors
######################

@doc Markdown.doc"""
    ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{Int})

Construct the torus invariant divisor on the normal toric variety `v` as linear combination of the torus invariant prime divisors of `v`. The coefficients of thi linear combination are passed as list of integers as first argument.

# Examples
```jldoctest
julia> show(ToricDivisor(toric_projective_space(2), [1,1,2]))
A torus invariant divisor on a normal toric variety
```
"""
function ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{Int})
    if length(coeffs) != pm_object(v).N_RAYS
        throw(ArgumentError("Number of coefficients needs to match number of prime divisors!"))
    end
    ptd = Polymake.fulton.TDivisor(COEFFICIENTS=coeffs)
    pmntv = pm_object(v)
    Polymake.add(pmntv, "DIVISOR", ptd)
    return ToricDivisor(ptd)
end
export ToricDivisor


######################
# 3: Special constructors
######################

@doc Markdown.doc"""
    DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{Int})

Construct the torus invariant divisor associated to a character of the normal toric variety `v`.

# Examples
```jldoctest
julia> show(DivisorOfCharacter(toric_projective_space(2), [1,2]))
A torus invariant divisor on a normal toric variety
```
"""
function DivisorOfCharacter(v::AbstractNormalToricVariety, character::Vector{Int})
    if length(character) != rank(character_lattice(v))
        throw(ArgumentError("Character must consist of $(rank(character_lattice(v))) integers!"))
    end
    f = map_from_character_to_principal_divisors(v)
    char = sum(character .* gens(domain(f)))
    coeffs = [Int(x) for x in transpose(f(char).coeff)][:,1]
    return ToricDivisor(v, coeffs)
end
export DivisorOfCharacter


######################
# 5: Attributes
######################

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
    pmtd = pm_tdivisor(td)
    return Polyhedron(pmtd.SECTION_POLYTOPE)
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
    return Vector{Int}(Polymake.common.primitive(pm_tdivisor(td).COEFFICIENTS))
end
export coefficients


###############################################################################
###############################################################################
### 5: Display
###############################################################################
###############################################################################
function Base.show(io::IO, td::ToricDivisor)
    print(io, "A torus invariant divisor on a normal toric variety")
end

