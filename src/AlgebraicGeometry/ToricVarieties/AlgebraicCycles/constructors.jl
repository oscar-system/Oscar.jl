##############################################
# 1: The Julia type for toric algebraic cycles
##############################################

@attributes mutable struct RationalEquivalenceClass
    v::NormalToricVarietyType
    p::MPolyQuoRingElem
    RationalEquivalenceClass(v::NormalToricVarietyType, p::MPolyQuoRingElem) = new(v, p)
end


####################################################
# 2: Generic constructors
####################################################

@doc raw"""
    rational_equivalence_class(v::NormalToricVarietyType, p::MPolyQuoRingElem)

Construct the rational equivalence class of algebraic cycles corresponding to a linear combination of cones.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> chow_ring(P2)
Quotient
  of multivariate polynomial ring in 3 variables x1, x2, x3
    over rational field
  by ideal (x1 - x3, x2 - x3, x1*x2*x3)

julia> (x1, x2, x3) = gens(chow_ring(P2))
3-element Vector{MPolyQuoRingElem{QQMPolyRingElem}}:
 x1
 x2
 x3

julia> rational_equivalence_class(P2, x1)
Rational equivalence classon a normal toric variety represented by V(x3)
```
"""
function rational_equivalence_class(v::NormalToricVarietyType, p::MPolyQuoRingElem)
    @req (is_simplicial(v) && is_complete(v)) "Currently, algebraic cycles are only supported for toric varieties that are simplicial and complete"
    @req parent(p) == chow_ring(v) "The polynomial must reside in the Chow ring of the toric variety"
    return RationalEquivalenceClass(v, p)
end


@doc raw"""
    rational_equivalence_class(v::NormalToricVarietyType, coefficients::Vector{T}) where {T <: IntegerUnion}

Construct the rational equivalence class of algebraic cycles corresponding to a linear combination of cones.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> rational_equivalence_class(P2, [6, 5, 4, 3, 2, 1])
Rational equivalence class on a normal toric variety represented by 15V(x1,x3)+6V(x3)
```
"""
function rational_equivalence_class(v::NormalToricVarietyType, coefficients::Vector{T}) where {T <: IntegerUnion}
    @req (is_simplicial(v) && is_complete(v)) "Currently, algebraic cycles are only supported for toric varieties that are simplicial and complete"
    @req length(coefficients) == n_cones(v) "The number of coefficients must match the number of all cones (but the trivial one) in the fan of the toric variety"
    mons = gens_of_rational_equivalence_classes(v)
    return RationalEquivalenceClass(v, sum(coefficients[i]*mons[i] for i in 1:length(coefficients)))
end


####################################################
# 3: Special constructors
####################################################

@doc raw"""
    rational_equivalence_class(d::ToricDivisor)

Construct the rational equivalence class of algebraic cycles corresponding to the toric divisor `d`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> d = toric_divisor(P2, [1, 2, 3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> rational_equivalence_class(d)
Rational equivalence class on a normal toric variety represented by 6V(x3)
```
"""
function rational_equivalence_class(d::ToricDivisor)
    v = toric_variety(d)
    if is_trivial(d)
        return RationalEquivalenceClass(v, zero(chow_ring(v)))
    end
    coeffs = coefficients(d)
    indets = gens(chow_ring(v))
    return RationalEquivalenceClass(v, sum(coeffs[k]*indets[k] for k in 1:length(coeffs)))
end


@doc raw"""
    rational_equivalence_class(c::ToricDivisorClass)

Construct the algebraic cycle corresponding to the toric divisor class `c`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> tdc = toric_divisor_class(P2, [2])
Divisor class on a normal toric variety

julia> rational_equivalence_class(tdc)
Rational equivalence class on a normal toric variety represented by 2V(x3)
```
"""
rational_equivalence_class(c::ToricDivisorClass) = rational_equivalence_class(toric_divisor(c))


@doc raw"""
    RationalEquivalenceClass(l::ToricLineBundle)

Construct the toric algebraic cycle corresponding to the toric line bundle `l`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(P2, [2])
Toric line bundle on a normal toric variety

julia> polynomial(rational_equivalence_class(l))
2*x3
```
"""
rational_equivalence_class(l::ToricLineBundle) = rational_equivalence_class(toric_divisor(l))


@doc raw"""
    rational_equivalence_class(cc::CohomologyClass)

Construct the toric algebraic cycle corresponding to the cohomology class `cc`.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> (x1, x2, x3) = gens(cohomology_ring(P2))
3-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 x1
 x2
 x3

julia> cc = cohomology_class(P2, x1+x2)
Cohomology class on a normal toric variety given by x1 + x2

julia> rational_equivalence_class(cc)
Rational equivalence class on a normal toric variety represented by 2V(x3)
```
"""
rational_equivalence_class(cc::CohomologyClass) = RationalEquivalenceClass(toric_variety(cc), polynomial(chow_ring(toric_variety(cc)), cc))


@doc raw"""
    rational_equivalence_class(sv::ClosedSubvarietyOfToricVariety)

Construct the rational equivalence class of algebraic
cycles of a closed subvariety of a normal toric variety.

# Examples
```jldoctest
julia> ntv = normal_toric_variety(Oscar.normal_fan(Oscar.cube(2)))
Normal toric variety

julia> set_coordinate_names(ntv, ["x1", "x2", "y1", "y2"]);

julia> (x1, x2, y1, y2) = gens(cox_ring(ntv))
4-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 y1
 y2

julia> sv = closed_subvariety_of_toric_variety(ntv, [x1^2+x1*x2+x2^2, y2])
Closed subvariety of a normal toric variety

julia> rational_equivalence_class(sv)
Rational equivalence class on a normal toric variety represented by 2V(x2,y2)
```
"""
function rational_equivalence_class(sv::ClosedSubvarietyOfToricVariety)
    v = toric_variety(sv)
    indets = gens(chow_ring(v))
    mons = [[m for m in monomials(p)][1] for p in gens(defining_ideal(sv))]
    expos = [matrix(ZZ, [k for k in AbstractAlgebra.exponent_vectors(mons[k])]) for k in 1:length(mons)]
    coeffs = [1 for i in 1:length(mons)]
    new_mons = MPolyQuoRingElem{QQMPolyRingElem}[]
    for k in 1:length(mons)
      mon = 1
      for j in 1:ncols(expos[k])
        if expos[k][1, j] != 0
          coeffs[k] = coeffs[k] * expos[k][1, j]
          mon = mon * indets[j]
        end
      end
      push!(new_mons, mon)
    end
    classes = [coeffs[k]*RationalEquivalenceClass(v, new_mons[k]) for k in 1:length(mons)]
    return prod(classes)
end


####################################################
# 4: Addition, subtraction and scalar multiplication
####################################################

function Base.:+(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    @req toric_variety(ac1) === toric_variety(ac2) "The rational equivalence classes must be defined on the same toric variety"
    return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) + polynomial(ac2))
end


function Base.:-(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    @req toric_variety(ac1) === toric_variety(ac2) "The rational equivalence classes must be defined on the same toric variety"
    return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) - polynomial(ac2))
end


Base.:*(c::QQFieldElem, ac::RationalEquivalenceClass) = RationalEquivalenceClass(toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac))
Base.:*(c::Rational{Int64}, ac::RationalEquivalenceClass) = RationalEquivalenceClass(toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac))
Base.:*(c::T, ac::RationalEquivalenceClass) where {T <: IntegerUnion} = RationalEquivalenceClass(toric_variety(ac), coefficient_ring(toric_variety(ac))(c) * polynomial(ac))


####################################################
# 5: Intersection product
####################################################

function Base.:*(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    @req toric_variety(ac1) === toric_variety(ac2) "The rational equivalence classes must be defined on the same toric variety"
    return RationalEquivalenceClass(toric_variety(ac1), polynomial(ac1) * polynomial(ac2))
end


Base.:^(ac::RationalEquivalenceClass, p::T) where {T <: IntegerUnion} = RationalEquivalenceClass(toric_variety(ac), polynomial(ac)^p)


function Base.:*(ac::RationalEquivalenceClass, sv::ClosedSubvarietyOfToricVariety)
    @req toric_variety(ac) === toric_variety(sv) "The rational equivalence class and the closed subvariety must be defined on the same toric variety"
    return ac * rational_equivalence_class(sv)
end


function Base.:*(sv::ClosedSubvarietyOfToricVariety, ac::RationalEquivalenceClass)
    @req toric_variety(ac) === toric_variety(sv) "The rational equivalence class and the closed subvariety must be defined on the same toric variety"
    return ac * rational_equivalence_class(sv)
end


function Base.:*(sv1::ClosedSubvarietyOfToricVariety, sv2::ClosedSubvarietyOfToricVariety)
    @req toric_variety(sv1) === toric_variety(sv2) "The closed subvarieties must be defined on the same toric variety"
    return rational_equivalence_class(sv1) * rational_equivalence_class(sv2)
end


####################################################
# 6: Equality and hash
####################################################

function Base.:(==)(ac1::RationalEquivalenceClass, ac2::RationalEquivalenceClass)
    return toric_variety(ac1) === toric_variety(ac2) && polynomial(ac1) == polynomial(ac2)
end

function Base.hash(ac::RationalEquivalenceClass, h::UInt) 
    b = 0xb5d4ac6b9084eb6e  % UInt
    h = hash(toric_variety(ac), h)
    h = hash(polynomial(ac), h)
    return xor(h, b)
end


####################################################
# 7: Display
####################################################

function Base.show(io::IO, ac::RationalEquivalenceClass)
    if is_trivial(ac)
        join(io, "Trivial rational equivalence class on a normal toric variety")
    else
      # otherwise, extract properties to represent the rational equivalence class
      r = representative(ac)
      coeffs = [c for c in AbstractAlgebra.coefficients(r)]
      expos = [matrix(ZZ, [k for k in AbstractAlgebra.exponent_vectors(m)]) for m in AbstractAlgebra.monomials(r)]
      indets = gens(chow_ring(toric_variety(ac)))

      # form string to be printed
      properties_string = String[]
      for i in 1:length(coeffs)
          m = String[]
          for j in 1:ncols(expos[i])
              for k in 1:expos[i][1, j]
                  push!(m, string(indets[j]))
              end
          end
          tmp = join(m, ",")
          if i == 1 && coeffs[i] == 1
              push!(properties_string, "Rational equivalence classon a normal toric variety represented by V($tmp)")
          elseif i == 1 && coeffs[i] != 1
              push!(properties_string, "Rational equivalence class on a normal toric variety represented by $(coeffs[i])V($tmp)")
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
