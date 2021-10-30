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
    if length(coeffs) != pm_ntv(v).N_RAYS
        throw(ArgumentError("Number of coefficients needs to match number of prime divisors!"))
    end
    ptd = Polymake.fulton.TDivisor(COEFFICIENTS=coeffs)
    pmntv = pm_ntv(v)
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
# 4: Properties
######################

@doc Markdown.doc"""
    iscartier(d::ToricDivisor)

Checks if the divisor `d` is Cartier.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> iscartier(td)
true
```
"""
iscartier(td::ToricDivisor) = pm_tdivisor(td).CARTIER::Bool
export iscartier


@doc Markdown.doc"""
    isprincipal(td::ToricDivisor) 

Determine whether the toric divisor `td` is principal.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isprincipal(td)
false
```
"""
isprincipal(td::ToricDivisor) = pm_tdivisor(td).PRINCIPAL::Bool
export isprincipal


@doc Markdown.doc"""
    isbasepoint_free(td::ToricDivisor) 

Determine whether the toric divisor `td` is basepoint free.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isbasepoint_free(td)
true
```
"""
isbasepoint_free(td::ToricDivisor) = pm_tdivisor(td).BASEPOINT_FREE::Bool
export isbasepoint_free


@doc Markdown.doc"""
    iseffective(td::ToricDivisor) 

Determine whether the toric divisor `td` is effective.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> iseffective(td)
true
```
"""
iseffective(td::ToricDivisor) = pm_tdivisor(td).EFFECTIVE::Bool
export iseffective


@doc Markdown.doc"""
    isintegral(td::ToricDivisor) 

Determine whether the toric divisor `td` is integral.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isintegral(td)
true
```
"""
isintegral(td::ToricDivisor) = pm_tdivisor(td).INTEGRAL::Bool
export isintegral


@doc Markdown.doc"""
    isample(td::ToricDivisor) 

Determine whether the toric divisor `td` is ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isample(td)
false
```
"""
isample(td::ToricDivisor) = pm_tdivisor(td).AMPLE::Bool
export isample


@doc Markdown.doc"""
    isvery_ample(td::ToricDivisor) 

Determine whether the toric divisor `td` is very ample.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isvery_ample(td)
false
```
"""
isvery_ample(td::ToricDivisor) = pm_tdivisor(td).VERY_AMPLE::Bool
export isvery_ample


@doc Markdown.doc"""
    isnef(td::ToricDivisor) 

Determine whether the toric divisor `td` is nef.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isnef(td)
true
```
"""
isnef(td::ToricDivisor) = pm_tdivisor(td).NEF::Bool
export isnef


@doc Markdown.doc"""
    isq_cartier(td::ToricDivisor) 

Determine whether the toric divisor `td` is Q-Cartier.
# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isq_cartier(td)
true
```
"""
isq_cartier(td::ToricDivisor) = pm_tdivisor(td).Q_CARTIER::Bool
export isq_cartier


@doc Markdown.doc"""
    isprime_divisor(td::ToricDivisor) 

Determine whether the toric divisor `td` is a prime divisor.

# Examples
```jldoctest
julia> H = hirzebruch_surface(4)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> td = ToricDivisor(H, [1,0,0,0])
A torus invariant divisor on a normal toric variety

julia> isprime_divisor(td)
true
```
"""
function isprime_divisor(td::ToricDivisor)
    pm_divisor = td.polymake_divisor
    coeffs = [Int(c) for c in Polymake.common.primitive(pm_divisor.COEFFICIENTS)]
    if sum(coeffs) != 1
        return false
    end
    return all(y -> (y == 1 || y == 0), coeffs)
end
export isprime_divisor


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
### 6: Display
###############################################################################
###############################################################################
function Base.show(io::IO, td::ToricDivisor)
    print(io, "A torus invariant divisor on a normal toric variety")
end

