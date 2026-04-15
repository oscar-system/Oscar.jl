########################
# 1: The Julia type for toric line bundles
########################

abstract type ToricCoherentSheaf end

@attributes mutable struct ToricLineBundle <: ToricCoherentSheaf
  toric_variety::NormalToricVarietyType
  picard_class::FinGenAbGroupElem
  function ToricLineBundle(
    toric_variety::NormalToricVarietyType, picard_class::FinGenAbGroupElem
  )
    @req parent(picard_class) === picard_group_with_map(toric_variety)[1] "The class must belong to the Picard group of the toric variety"
    return new(toric_variety, picard_class)
  end
end

########################
# 2: Generic constructors
########################

@doc raw"""
    toric_line_bundle(v::NormalToricVarietyType, picard_class::FinGenAbGroupElem)

Construct the line bundle on the abstract normal toric variety with given class
in the Picard group of the toric variety in question.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> pc = picard_group_with_map(P2)[1];

julia> l = toric_line_bundle(P2, pc([1]))
Toric line bundle on a normal toric variety
```
"""
toric_line_bundle(v::NormalToricVarietyType, picard_class::FinGenAbGroupElem) =
  ToricLineBundle(
    v, picard_class
  )

@doc raw"""
    toric_line_bundle(v::NormalToricVarietyType, picard_class::Vector{T}) where {T <: IntegerUnion}

Construct the line bundle on the abstract normal toric variety `v` with class `c`
in the Picard group of `v`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, [ZZRingElem(2)])
Toric line bundle on a normal toric variety
```
"""
function toric_line_bundle(
  v::NormalToricVarietyType, picard_class::Vector{T}
) where {T<:IntegerUnion}
  return ToricLineBundle(v, picard_group_with_map(v)[1](picard_class))
end

########################
# 3: Special constructor
########################

@doc raw"""
    toric_line_bundle(v::NormalToricVarietyType, d::ToricDivisor)

Construct the toric variety associated to a (Cartier) torus-invariant divisor `d` on the normal toric variety `v`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l = toric_line_bundle(v, toric_divisor(v, [1, 2, 3]))
Toric line bundle on a normal toric variety
```
"""
function toric_line_bundle(v::NormalToricVarietyType, d::ToricDivisor)
  @req is_cartier(d) "The toric divisor must be Cartier to define a toric line bundle"
  f = map_from_torusinvariant_cartier_divisor_group_to_picard_group(v)
  g = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v)
  cartier_d = preimage(
    g, sum(coefficients(d) .* gens(torusinvariant_weil_divisor_group(v)))
  )
  picard_class = f(cartier_d)
  l = ToricLineBundle(v, picard_class)
  set_attribute!(l, :toric_divisor, d)
  return l
end

@doc raw"""
    toric_line_bundle(d::ToricDivisor)

Construct the toric variety associated to a (Cartier) torus-invariant divisor `d`.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> d = toric_divisor(v, [1, 2, 3]);

julia> l = toric_line_bundle(d)
Toric line bundle on a normal toric variety
```
"""
toric_line_bundle(d::ToricDivisor) = toric_line_bundle(toric_variety(d), d)

@doc raw"""
    toric_line_bundle(v::NormalToricVarietyType, dc::ToricDivisorClass)

Construct the toric variety associated to a divisor class in the class group
of a toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> d = toric_divisor(v, [1, 2, 3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> dc = toric_divisor_class(d)
Divisor class on a normal toric variety

julia> l = toric_line_bundle(v, dc)
Toric line bundle on a normal toric variety
```
"""
function toric_line_bundle(v::NormalToricVarietyType, dc::ToricDivisorClass)
  f = map_from_torusinvariant_weil_divisor_group_to_class_group(v)
  g = map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v)
  h = map_from_torusinvariant_cartier_divisor_group_to_picard_group(v)
  cartier_class = preimage(g * f, divisor_class(dc))
  td = toric_divisor(v, _vec(g(cartier_class).coeff))
  l = ToricLineBundle(v, h(cartier_class))
  set_attribute!(td, :is_cartier, true)
  set_attribute!(l, :toric_divisor, td)
  return l
end

@doc raw"""
    toric_line_bundle(dc::ToricDivisorClass)

Construct the toric variety associated to a divisor class in the class group
of a toric variety.

# Examples
```jldoctest
julia> v = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> d = toric_divisor(v, [1, 2, 3])
Torus-invariant, non-prime divisor on a normal toric variety

julia> dc = toric_divisor_class(d)
Divisor class on a normal toric variety

julia> l = toric_line_bundle(dc)
Toric line bundle on a normal toric variety
```
"""
toric_line_bundle(dc::ToricDivisorClass) = toric_line_bundle(toric_variety(dc), dc)

########################
# 4: Tensor products
########################

function Base.:*(l1::ToricLineBundle, l2::ToricLineBundle)
  @req toric_variety(l1) === toric_variety(l2) "The line bundles must be defined on the same toric variety"
  return toric_line_bundle(toric_variety(l1), picard_class(l1) + picard_class(l2))
end
Base.:inv(l::ToricLineBundle) = toric_line_bundle(toric_variety(l), (-1) * picard_class(l))
Base.:^(l::ToricLineBundle, p::ZZRingElem) = toric_line_bundle(
  toric_variety(l), p * picard_class(l)
)
Base.:^(l::ToricLineBundle, p::Int) = l^ZZRingElem(p)

########################
# 5: Equality and hash
########################

function Base.:(==)(l1::ToricLineBundle, l2::ToricLineBundle)
  return toric_variety(l1) === toric_variety(l2) && picard_class(l1) == picard_class(l2)
end

function Base.hash(l::ToricLineBundle, h::UInt)
  b = 0xa2b0a2cd60a8ffbf % UInt
  h = hash(toric_variety(l), h)
  h = hash(picard_class(l), h)
  return xor(h, b)
end

########################
# 6: Display
########################

function Base.show(io::IO, line_bundle::ToricLineBundle)

  # initiate properties string
  properties_string = ["Toric"]

  # collect known properties
  if has_attribute(line_bundle, :toric_divisor)
    td = toric_divisor(line_bundle)
    push_attribute_if_exists!(properties_string, td, :is_principal, "trivial")
    push_attribute_if_exists!(properties_string, td, :is_basepoint_free, "basepoint-free")
    ample_cb!(a, b) = push_attribute_if_exists!(a, b, :is_ample, "ample")
    push_attribute_if_exists!(
      properties_string, td, :is_very_ample, "very-ample"; callback=(ample_cb!)
    )
  end

  # print
  push!(properties_string, "line bundle on a normal toric variety")
  join(io, properties_string, ", ", " ")
end
