##############################################
# 1: The Julia type for toric algebraic cycles
##############################################

@attributes mutable struct RationalEquivalenceClass
    v::AbstractNormalToricVariety
    p::MPolyQuoElem
    RationalEquivalenceClass(v::AbstractNormalToricVariety,p::MPolyQuoElem) = new(v,p)
end
export RationalEquivalenceClass


####################################################
# 2: Generic constructors
####################################################

@doc Markdown.doc"""
    RationalEquivalenceClass(v::AbstractNormalToricVariety, coefficients::Vector{T}) where {T <: IntegerUnion}

Construct the rational equivalence class of algebraic cycles corresponding to a linear combination of cones.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(P2, [1,2,3])
A torus-invariant, non-prime divisor on a normal toric variety

julia> RationalEquivalenceClass(d)
A rational equivalence class on a normal toric variety represented by 6V(x3)
```
"""
function RationalEquivalenceClass(v::AbstractNormalToricVariety, coefficients::Vector{T}) where {T <: IntegerUnion}
    if !(is_complete(v) && is_simplicial(v))
        throw(ArgumentError("Currently, the Chow ring is only supported for toric varieties that are both complete and simplicial"))
    end
    if length(coefficients) != length(cones(v))
        throw(ArgumentError("The number of coefficients must match the number of all cones (but the trivial one) in the fan of the toric variety"))
    end
    mons = gens_of_rational_equivalence_classes(v)
    return RationalEquivalenceClass(v,sum(coefficients[i]*mons[i] for i in 1:length(coefficients)))
end


####################################################
# 3: Special constructors
####################################################

@doc Markdown.doc"""
    RationalEquivalenceClass(d::ToricDivisor)

Construct the rational equivalence class of algebraic cycles corresponding to the toric divisor `d`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> d = ToricDivisor(P2, [1,2,3])
A torus-invariant, non-prime divisor on a normal toric variety

julia> RationalEquivalenceClass(d)
A rational equivalence class on a normal toric variety represented by 6V(x3)
```
"""
function RationalEquivalenceClass(d::ToricDivisor)
    v = toric_variety(d)
    if istrivial(d)
        return RationalEquivalenceClass(v,zero(chow_ring(v)))
    end
    coeffs = coefficients(d)
    indets = gens(chow_ring(v))
    return RationalEquivalenceClass(v, sum(coeffs[k]*indets[k] for k in 1:length(coeffs)))
end


@doc Markdown.doc"""
    RationalEquivalenceClass(c::ToricDivisorClass)

Construct the algebraic cycle corresponding to the toric divisor class `c`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> tdc = ToricDivisorClass(P2, [2])
A divisor class on a normal toric variety

julia> RationalEquivalenceClass(tdc)
A rational equivalence class on a normal toric variety represented by 2V(x3)
```
"""
RationalEquivalenceClass(c::ToricDivisorClass) = RationalEquivalenceClass(toric_divisor(c))


@doc Markdown.doc"""
    RationalEquivalenceClass(l::ToricLineBundle)

Construct the toric algebraic cycle corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> l = ToricLineBundle(P2, [2])
A toric line bundle on a normal toric variety

julia> polynomial(RationalEquivalenceClass(l))
2*x3
```
"""
RationalEquivalenceClass(l::ToricLineBundle) = RationalEquivalenceClass(toric_divisor(l))


@doc Markdown.doc"""
    RationalEquivalenceClass(cc::CohomologyClass)

Construct the toric algebraic cycle corresponding to the cohomology class `cc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> (x1,x2,x3)=gens(cohomology_ring(P2))
3-element Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}:
 x1
 x2
 x3

julia> cc = CohomologyClass(P2,x1+x2)
A cohomology class on a normal toric variety given by x1 + x2

julia> RationalEquivalenceClass(cc)
A rational equivalence class on a normal toric variety represented by 2V(x3)
```
"""
RationalEquivalenceClass(cc::CohomologyClass) = RationalEquivalenceClass(toric_variety(cc), polynomial(cc))


@doc Markdown.doc"""
    RationalEquivalenceClass(cc::ClosedSubvarietyOfToricVariety)

Construct the rational equivalence class of algebraic
cycles of a closed subvariety of a normal toric variety.

# Examples
```jldoctest
julia> ntv = NormalToricVariety(Oscar.normal_fan(Oscar.cube(2)))
A normal toric variety

julia> set_coordinate_names(ntv,["x1","x2","y1","y2"]);

julia> (x1, x2, y1, y2) = gens(cox_ring(ntv))
4-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x1
 x2
 y1
 y2

julia> sv = ClosedSubvarietyOfToricVariety(ntv, [x1^2+x1*x2+x2^2,y2])
A closed subvariety of a normal toric variety

julia> RationalEquivalenceClass(sv)
A rational equivalence class on a normal toric variety represented by 2V(x2,y2)
```
"""
function RationalEquivalenceClass(sv::ClosedSubvarietyOfToricVariety)
    v = toric_variety(sv)
    indets = gens(chow_ring(v))
    mons = [[m for m in monomials(p)][1] for p in gens(defining_ideal(sv))]
    expos = [matrix(ZZ,[k for k in exponent_vectors(mons[k])]) for k in 1:length(mons)]
    coeffs = [1 for i in 1:length(mons)]
    new_mons = MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}[]
    for k in 1:length(mons)
      mon = 1
      for j in 1:ncols(expos[k])
        if expos[k][1,j] != 0
          coeffs[k] = coeffs[k] * expos[k][1,j]
          mon = mon * indets[j]
        end
      end
      push!(new_mons,mon)
    end
    classes = [coeffs[k]*RationalEquivalenceClass(v,new_mons[k]) for k in 1:length(mons)]
    return prod(classes)
end

####################################################
# 4: Addition, subtraction and scalar multiplication
####################################################

function Base.:+(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    if toric_variety(ac1) !== toric_variety(ac2)
        throw(ArgumentError("The rational equivalence classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) + polynomial(ac2))
end


function Base.:-(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    if toric_variety(ac1) !== toric_variety(ac2)
        throw(ArgumentError("The rational equivalence classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) - polynomial(ac2))
end


Base.:*(c::fmpq, ac::RationalEquivalenceClass) = RationalEquivalenceClass(toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac))
Base.:*(c::Rational{Int64}, ac::RationalEquivalenceClass) = RationalEquivalenceClass(toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac))
Base.:*(c::T, ac::RationalEquivalenceClass) where {T <: IntegerUnion} = RationalEquivalenceClass(toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac))


####################################################
# 5: Intersection product
####################################################

function Base.:*(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    if toric_variety(ac1) !== toric_variety(ac2)
        throw(ArgumentError("The rational equivalence classes must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) * polynomial(ac2))
end


Base.:^(ac::RationalEquivalenceClass, p::T) where {T <: IntegerUnion} = RationalEquivalenceClass(toric_variety(ac), polynomial(ac)^p)


function Base.:*(ac::RationalEquivalenceClass, sv::ClosedSubvarietyOfToricVariety)
    if toric_variety(ac) !== toric_variety(sv)
        throw(ArgumentError("The rational equivalence class and the closed subvariety must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ac * RationalEquivalenceClass(sv)
end


function Base.:*(sv::ClosedSubvarietyOfToricVariety, ac::RationalEquivalenceClass)
    if toric_variety(ac) !== toric_variety(sv)
        throw(ArgumentError("The rational equivalence class and the closed subvariety must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ac * RationalEquivalenceClass(sv)
end


function Base.:*(sv1::ClosedSubvarietyOfToricVariety, sv2::ClosedSubvarietyOfToricVariety)
    if toric_variety(sv1) !== toric_variety(sv2)
        throw(ArgumentError("The closed subvarieties must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return RationalEquivalenceClass(sv1) * RationalEquivalenceClass(sv2)
end


####################################################
# 6: Equality
####################################################

function Base.:(==)(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    return toric_variety(ac1) === toric_variety(ac2) && iszero(polynomial(ac1-ac2))
end


####################################################
# 7: Display
####################################################

function Base.show(io::IO, ac::RationalEquivalenceClass)
    if is_trivial(ac)
        join(io, "A trivial rational equivalence class on a normal toric variety")
    else
      # otherwise, extract properties to represent the rational equivalence class
      r = representative(ac)
      coeffs = [c for c in coefficients(r)]
      expos = [matrix(ZZ,[k for k in exponent_vectors(m)]) for m in monomials(r)]
      indets = gens(chow_ring(toric_variety(ac)))

      # form string to be printed
      properties_string = String[]
      for i in 1:length(coeffs)
          m = String[]
          for j in 1:ncols(expos[i])
              for k in 1:expos[i][1,j]
                  push!(m,string(indets[j]))
              end
          end
          tmp = join(m,",")
          if i == 1 && coeffs[i] == 1
              push!(properties_string, "A rational equivalence classon a normal toric variety represented by V($tmp)")
          elseif i == 1 && coeffs[i] != 1
              push!(properties_string, "A rational equivalence class on a normal toric variety represented by $(coeffs[i])V($tmp)")
          elseif i > 1 && coeffs[i] == 1
              push!(properties_string, "+V($tmp)")
          elseif i > 1 && coeffs[i] > 0 && coeffs[i] != 1
              push!(properties_string, "+$(coeffs[i])V($tmp)")
          elseif i > 1 && coeffs[i] < 0
              push!(properties_string, "$(coeffs[i])V($tmp)")
          end
      end
      # print information
      join(io, properties_string)
    end
end
