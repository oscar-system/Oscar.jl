# Deprecated in 0.10.*
@deprecate automorphisms(x::Graph) automorphism_group_generators(x)

# Deprecated after 0.11.3
@deprecate is_finite_order(x::GAPGroupElem) is_finiteorder(x)
@deprecate elements(C::GroupConjClass) collect(C)
@deprecate conjugate_subgroup(G::GAPGroup, x::GAPGroupElem) conjugate_group(G, x)
@deprecate gap_perm(L::AbstractVector{<:IntegerUnion}) perm(L::AbstractVector{<:IntegerUnion})
@deprecate elements(C::GroupCoset) collect(C)
@deprecate elements(C::GroupDoubleCoset) collect(C)
@deprecate map_from_character_to_principal_divisors(v::AbstractNormalToricVariety) map_from_character_lattice_to_torusinvariant_weil_divisor_group(v)
@deprecate map_from_weil_divisors_to_class_group(v::AbstractNormalToricVariety) map_from_torusinvariant_weil_divisor_group_to_class_group(v)
@deprecate map_from_cartier_divisor_group_to_torusinvariant_divisor_group(v::AbstractNormalToricVariety) map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v)
@deprecate map_from_cartier_divisor_group_to_picard_group(v::AbstractNormalToricVariety) map_from_torusinvariant_cartier_divisor_group_to_picard_group(v)
@deprecate cartier_divisor_group(v::AbstractNormalToricVariety) torusinvariant_cartier_divisor_group(v)
@deprecate torusinvariant_divisor_group(v::AbstractNormalToricVariety) torusinvariant_weil_divisor_group(v)
@deprecate StructureSheaf(v::AbstractNormalToricVariety) structure_sheaf
@deprecate morphism_on_cartier_divisor_group(tm::ToricMorphism) morphism_on_torusinvariant_cartier_divisor_group(tm)
@deprecate bounded(Obj::Polyhedron) is_bounded(Obj)
@deprecate vf_group(P::Polyhedron) automorphism_group(P; action = :on_facets)
@deprecate birkhoff(n::Integer; even::Bool = false) birkhoff_polytope(n; even=even)
@deprecate cross(args...) cross_polytope(args...)
@deprecate gelfand_tsetlin(lambda::AbstractVector) gelfand_tsetlin_polytope(lambda)


# Deprecate after 0.11.4
function AffineNormalToricVariety(C::Cone; set_attributes::Bool = true)
    Base.depwarn("'AffineNormalToricVariety(C::Cone; set_attributes::Bool = true)' is deprecated, use "*
    "'affine_normal_toric_variety(C::Cone; set_attributes::Bool = true)' instead.", :AffineNormalToricVariety)
    affine_normal_toric_variety(C; set_attributes = set_attributes)
end

function AffineNormalToricVariety(v::NormalToricVariety; set_attributes::Bool = true)
    Base.depwarn("'AffineNormalToricVariety(v::NormalToricVariety; set_attributes::Bool = true)' is deprecated, use "*
    "'affine_normal_toric_variety(v::NormalToricVariety; set_attributes::Bool = true)' instead.", :AffineNormalToricVariety)
    affine_normal_toric_variety(v; set_attributes = set_attributes)
end

function NormalToricVariety(C::Cone; set_attributes::Bool = true)
    Base.depwarn("'NormalToricVariety(C::Cone; set_attributes::Bool = true)' is deprecated, use "*
    "'normal_toric_variety(C::Cone; set_attributes::Bool = true)' instead.", :NormalToricVariety)
    normal_toric_variety(C; set_attributes = set_attributes)
end

function NormalToricVariety(rays::Vector{Vector{Int64}}, max_cones::Vector{Vector{Int64}}; non_redundant::Bool = false, set_attributes::Bool = true)
    Base.depwarn("'NormalToricVariety(rays::Vector{Vector{Int64}}, max_cones::Vector{Vector{Int64}}; "*
    "non_redundant::Bool = false, set_attributes::Bool = true)' is deprecated, use "*
    "'normal_toric_variety(rays::Vector{Vector{Int64}}, max_cones::Vector{Vector{Int64}}; "*
    "non_redundant::Bool = false, set_attributes::Bool = true)' instead.", :NormalToricVariety)
    normal_toric_variety(rays, max_cones; non_redundant = non_redundant, set_attributes = set_attributes)
end

function NormalToricVariety(PF::PolyhedralFan; set_attributes::Bool = true)
    Base.depwarn("'NormalToricVariety(PF::PolyhedralFan; set_attributes::Bool = true)' is deprecated, use"*
    "'normal_toric_variety(PF::PolyhedralFan; set_attributes::Bool = true)' instead.", :NormalToricVariety)
    normal_toric_variety(PF; set_attributes = set_attributes)
end

function NormalToricVariety(P::Polyhedron; set_attributes::Bool = true)
    Base.depwarn("'NormalToricVariety(P::Polyhedron; set_attributes::Bool = true)' is deprecated, use"*
    "'normal_toric_variety(P::Polyhedron; set_attributes::Bool = true)' instead.", :NormalToricVariety)
    normal_toric_variety(P; set_attributes = set_attributes)
end

@deprecate NormalToricVarietiesFromStarTriangulations(P::Polyhedron; set_attributes::Bool = true) normal_toric_varieties_from_star_triangulations(P; set_attributes = set_attributes)
@deprecate NormalToricVarietyFromGLSM(charges::fmpz_mat; set_attributes::Bool = true) normal_toric_varieties_from_glsm(charges; set_attributes = set_attributes)

function RationalEquivalenceClass(v::AbstractNormalToricVariety, coefficients::Vector{T}) where {T <: IntegerUnion}
    Base.depwarn("'RationalEquivalenceClass(v::AbstractNormalToricVariety, coefficients::Vector{T}) where {T <: IntegerUnion}'"*
    " is deprecated, use 'rational_equivalence_class(v::AbstractNormalToricVariety, coefficients::Vector{T}) "*
    "where {T <: IntegerUnion}' instead.", RationalEquivalenceClass)
    rational_equivalence_class(v, coefficients)
end

function RationalEquivalenceClass(d::ToricDivisor)
    Base.depwarn("'RationalEquivalenceClass(d::ToricDivisor)' is deprecated, use "*
    "'rational_equivalence_class(d::ToricDivisor)' instead.", :RationalEquivalenceClass)
    rational_equivalence_class(d)
end

function RationalEquivalenceClass(c::ToricDivisorClass)
    Base.depwarn("'RationalEquivalenceClass(c::ToricDivisorClass)' is deprecated, use "*
    "'rational_equivalence_class(c::ToricDivisorClass)' instead.", :RationalEquivalenceClass)
    rational_equivalence_class(c)
end

function RationalEquivalenceClass(l::ToricLineBundle)
    Base.depwarn("'RationalEquivalenceClass(l::ToricLineBundle)' is deprecated, use "*
    "'rational_equivalence_class(l::ToricLineBundle)' instead.", :RationalEquivalenceClass)
    rational_equivalence_class(l)
end

function RationalEquivalenceClass(cc::CohomologyClass)
    Base.depwarn("'RationalEquivalenceClass(cc::CohomologyClass)' is deprecated, use "*
    "'rational_equivalence_class(cc::CohomologyClass)' instead.", :RationalEquivalenceClass)
    rational_equivalence_class(cc)
end

function RationalEquivalenceClass(sv::ClosedSubvarietyOfToricVariety)
    Base.depwarn("'RationalEquivalenceClass(sv::ClosedSubvarietyOfToricVariety)' is deprecated, use "*
    "'rational_equivalence_class(sv::ClosedSubvarietyOfToricVariety)' instead.", :RationalEquivalenceClass)
    rational_equivalence_class(sv)
end

function CohomologyClass(d::ToricDivisor)
    Base.depwarn("'CohomologyClass(d::ToricDivisor)' is deprecated, use "*
    "'cohomology_class(d::ToricDivisor)' instead.", :CohomologyClass)
    cohomology_class(d)
end

function CohomologyClass(c::ToricDivisorClass)
    Base.depwarn("'CohomologyClass(c::ToricDivisorClass)' is deprecated, use "*
    "'cohomology_class(c::ToricDivisorClass)' instead.", :CohomologyClass)
    cohomology_class(c)
end

function CohomologyClass(l::ToricLineBundle)
    Base.depwarn("'CohomologyClass(l::ToricLineBundle)' is deprecated, use "*
    "'cohomology_class(l::ToricLineBundle)' instead.", :CohomologyClass)
    cohomology_class(l)
end

function CyclicQuotientSingularity(n::T, q::T) where {T <: IntegerUnion}
    Base.depwarn("'CyclicQuotientSingularity(n::fmpz, q::fmpz)' is deprecated, use "*
    "'cyclic_quotient_singularity(n::fmpz, q::fmpz)' instead.", :CyclicQuotientSingularity)
    cyclic_quotient_singularity(n, q)
end

function ClosedSubvarietyOfToricVariety(toric_variety::AbstractNormalToricVariety, defining_polynomials::Vector{MPolyDecRingElem{fmpq, fmpq_mpoly}})
    Base.depwarn("ClosedSubvarietyOfToricVariety(toric_variety::AbstractNormalToricVariety, "*
    "defining_polynomials::Vector{MPolyDecRingElem{fmpq, fmpq_mpoly}}) is deprecated, use "*
    "'closed_subvariety_of_toric_variety(toric_variety::AbstractNormalToricVariety, "*
    "defining_polynomials::Vector{MPolyDecRingElem{fmpq, fmpq_mpoly}})' instead.", :ClosedSubvarietyOfToricVariety)
    closed_subvariety_of_toric_variety(toric_variety, defining_polynomials)
end

function ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion}
    Base.depwarn("'ToricDivisor(v::AbstractNormalToricVariety, coeffs::Vector{T}) where {T <: IntegerUnion}' "*
    "is deprecated, use 'toric_divisor(v::AbstractNormalToricVariety, coeffs::Vector{T}) "*
    "where {T <: IntegerUnion}' instead.", :ToricDivisor)
    toric_divisor(v, coeffs)
end
