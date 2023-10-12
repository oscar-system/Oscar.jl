# Deprecated in 0.10.*
@deprecate automorphisms(x::Graph) automorphism_group_generators(x)

# Deprecated after 0.11.3
@deprecate is_finite_order(x::GAPGroupElem) is_finiteorder(x)
@deprecate elements(C::GroupConjClass) collect(C)
@deprecate conjugate_subgroup(G::GAPGroup, x::GAPGroupElem) conjugate_group(G, x)
@deprecate gap_perm(L::AbstractVector{<:IntegerUnion}) perm(L::AbstractVector{<:IntegerUnion})
@deprecate elements(C::GroupCoset) collect(C)
@deprecate elements(C::GroupDoubleCoset) collect(C)
@deprecate map_from_character_to_principal_divisors(v::NormalToricVarietyType) map_from_character_lattice_to_torusinvariant_weil_divisor_group(v)
@deprecate map_from_weil_divisors_to_class_group(v::NormalToricVarietyType) map_from_torusinvariant_weil_divisor_group_to_class_group(v)
@deprecate map_from_cartier_divisor_group_to_torusinvariant_divisor_group(v::NormalToricVarietyType) map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(v)
@deprecate map_from_cartier_divisor_group_to_picard_group(v::NormalToricVarietyType) map_from_torusinvariant_cartier_divisor_group_to_picard_group(v)
@deprecate cartier_divisor_group(v::NormalToricVarietyType) torusinvariant_cartier_divisor_group(v)
@deprecate torusinvariant_divisor_group(v::NormalToricVarietyType) torusinvariant_weil_divisor_group(v)
@deprecate StructureSheaf(v::NormalToricVarietyType) structure_sheaf
@deprecate morphism_on_cartier_divisor_group(tm::ToricMorphism) morphism_on_torusinvariant_cartier_divisor_group(tm)
@deprecate bounded(Obj::Polyhedron) is_bounded(Obj)
@deprecate vf_group(P::Polyhedron) automorphism_group(P; action = :on_facets)
@deprecate birkhoff(n::Integer; even::Bool = false) birkhoff_polytope(n; even=even)
@deprecate cross(args...) cross_polytope(args...)
@deprecate gelfand_tsetlin(lambda::AbstractVector) gelfand_tsetlin_polytope(lambda)


# Deprecated after 0.11.4
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
    "'normal_toric_variety(PF::PolyhedralFan)' instead.", :NormalToricVariety)
    normal_toric_variety(PF)
end

function NormalToricVariety(P::Polyhedron; set_attributes::Bool = true)
    Base.depwarn("'NormalToricVariety(P::Polyhedron; set_attributes::Bool = true)' is deprecated, use"*
    "'normal_toric_variety(P::Polyhedron; set_attributes::Bool = true)' instead.", :NormalToricVariety)
    normal_toric_variety(P; set_attributes = set_attributes)
end

@deprecate NormalToricVarietiesFromStarTriangulations(P::Polyhedron; set_attributes::Bool = true) normal_toric_varieties_from_star_triangulations(P; set_attributes = set_attributes)
@deprecate NormalToricVarietyFromGLSM(charges::ZZMatrix; set_attributes::Bool = true) normal_toric_varieties_from_glsm(charges; set_attributes = set_attributes)

function RationalEquivalenceClass(v::NormalToricVarietyType, coefficients::Vector{T}) where {T <: IntegerUnion}
    Base.depwarn("'RationalEquivalenceClass(v::NormalToricVarietyType, coefficients::Vector{T}) where {T <: IntegerUnion}'"*
    " is deprecated, use 'rational_equivalence_class(v::NormalToricVarietyType, coefficients::Vector{T}) "*
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
    Base.depwarn("'CyclicQuotientSingularity(n::ZZRingElem, q::ZZRingElem)' is deprecated, use "*
    "'cyclic_quotient_singularity(n::ZZRingElem, q::ZZRingElem)' instead.", :CyclicQuotientSingularity)
    cyclic_quotient_singularity(n, q)
end

function ClosedSubvarietyOfToricVariety(toric_variety::NormalToricVarietyType, defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}})
    Base.depwarn("ClosedSubvarietyOfToricVariety(toric_variety::NormalToricVarietyType, "*
    "defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}) is deprecated, use "*
    "'closed_subvariety_of_toric_variety(toric_variety::NormalToricVarietyType, "*
    "defining_polynomials::Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}})' instead.", :ClosedSubvarietyOfToricVariety)
    closed_subvariety_of_toric_variety(toric_variety, defining_polynomials)
end

function ToricDivisor(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion}
    Base.depwarn("'ToricDivisor(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion}' "*
    "is deprecated, use 'toric_divisor(v::NormalToricVarietyType, coeffs::Vector{T}) "*
    "where {T <: IntegerUnion}' instead.", :ToricDivisor)
    toric_divisor(v, coeffs)
end

@deprecate DivisorOfCharacter(v::NormalToricVarietyType, character::Vector{T}) where {T <: IntegerUnion} divisor_of_character(v, character)

function ToricDivisorClass(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion}
    Base.depwarn("'ToricDivisorClass(v::NormalToricVarietyType, coeffs::Vector{T}) where {T <: IntegerUnion}' "*
    "is deprecated, use 'toric_divisor_class(v::NormalToricVarietyType, coeffs::Vector{T}) "*
    "where {T <: IntegerUnion}' instead.", :ToricDivisorClass)
    toric_divisor_class(v, coeffs)
end

function ToricDivisorClass(td::ToricDivisor)
    Base.depwarn("'ToricDivisorClass(td::ToricDivisor)' is deprecated, use "*
    "'toric_divisor_class(td::ToricDivisor) instead.", :ToricDivisorClass)
    toric_divisor_class(td)
end

function ToricLineBundle(v::NormalToricVarietyType, c::Vector{T}) where {T <: IntegerUnion}
    Base.depwarn("'ToricLineBundle(v::NormalToricVarietyType, c::Vector{T}) where {T <: IntegerUnion}'"*
    " is deprecated, use 'toric_line_bundle(v::NormalToricVarietyType, c::Vector{T}) where {T <: IntegerUnion}' "*
    "instead.", :ToricLineBundle)
    toric_line_bundle(v, c)
end

function ToricLineBundle(v::NormalToricVarietyType, d::ToricDivisor)
    Base.depwarn("'ToricLineBundle(v::NormalToricVarietyType, d::ToricDivisor)'"*
    " is deprecated, use 'toric_line_bundle(v::NormalToricVarietyType, d::ToricDivisor)' "*
    "instead.", :ToricLineBundle)
    toric_line_bundle(v, d)
end

function ToricLineBundle(d::ToricDivisor)
    Base.depwarn("'ToricLineBundle(d::ToricDivisor)'"*
    " is deprecated, use 'toric_line_bundle(d::ToricDivisor)' instead", :ToricLineBundle)
    toric_line_bundle(d)
end

function ToricMorphism(domain::NormalToricVarietyType, mapping_matrix::Vector{Vector{T}}, codomain::T2=nothing) where {T <: IntegerUnion, T2 <: Union{NormalToricVarietyType, Nothing}}
    Base.depwarn("'ToricMorphism(domain::NormalToricVarietyType, mapping_matrix::Vector{Vector{T}}, codomain::T2=nothing) "*
    "where {T <: IntegerUnion, T2 <: Union{NormalToricVarietyType, Nothing}}' is depcreated, use "*
    "'toric_morphism(domain::NormalToricVarietyType, mapping_matrix::Vector{Vector{T}}, codomain::T2=nothing) "*
    "where {T <: IntegerUnion, T2 <: Union{NormalToricVarietyType, Nothing}})' instead", :ToricMorphism)
    toric_morphism(domain, mapping_matrix, codomain)
end

function ToricMorphism(domain::NormalToricVarietyType, mapping_matrix::Matrix{T}, codomain::T2=nothing) where {T <: IntegerUnion, T2 <: Union{NormalToricVarietyType, Nothing}}
    Base.depwarn("'ToricMorphism(domain::NormalToricVarietyType, mapping_matrix::Matrix{T}, codomain::T2=nothing) "*
    "where {T <: IntegerUnion, T2 <: Union{NormalToricVarietyType, Nothing}}' is deprecated, use "*
    "'toric_morphism(domain::NormalToricVarietyType, mapping_matrix::Matrix{T}, codomain::T2=nothing) "*
    "where {T <: IntegerUnion, T2 <: Union{NormalToricVarietyType, Nothing}}' instead", :ToricMorphism)
    toric_morphism(domain, mapping_matrix, codomain)
end

function ToricMorphism(domain::NormalToricVarietyType, mapping_matrix::ZZMatrix, codomain::T=nothing) where {T <: Union{NormalToricVarietyType, Nothing}}
    Base.depwarn("'ToricMorphism(domain::NormalToricVarietyType, mapping_matrix::ZZMatrix, codomain::T=nothing) "*
    "where {T <: Union{NormalToricVarietyType, Nothing}}' is deprecated, use "*
    "'toric_morphism(domain::NormalToricVarietyType, mapping_matrix::ZZMatrix, codomain::T=nothing) "*
    "where {T <: Union{NormalToricVarietyType, Nothing}}' instead", :ToricMorphism)
    toric_morphism(domain, mapping_matrix, codomain)
end


@deprecate ToricIdentityMorphism(v::NormalToricVarietyType) toric_identity_morphism(v)

@deprecate induced_class_function induce
@deprecate radical_subgroup solvable_radical
@deprecate has_radical_subgroup has_solvable_radical
@deprecate set_radical_subgroup set_solvable_radical

# `is_characteristic` had the wrong order of arguments,
# see https://github.com/oscar-system/Oscar.jl/issues/1793
@deprecate is_characteristic(G::T, H::T) where T <: GAPGroup is_characteristic_subgroup(H, G)

# `is_maximal` had the wrong order of arguments,
# see https://github.com/oscar-system/Oscar.jl/issues/1793
@deprecate is_maximal(G::T, H::T) where T <: GAPGroup is_maximal_subgroup(H, G)

# `is_normal` had the wrong order of arguments,
# see https://github.com/oscar-system/Oscar.jl/issues/1793
@deprecate is_normal(G::T, H::T) where T <: GAPGroup is_normalized_by(H, G)


# Deprecated after 0.12.0
@deprecate contains(P::Polyhedron, v::AbstractVector) Base.in(v, P)
@deprecate contains(C::Cone, v::AbstractVector) Base.in(v, C)
@deprecate blowup_on_ith_minimal_torus_orbit(v::NormalToricVarietyType, n::Int, coordinate_name::String; set_attributes::Bool = true) blow_up(v, n; coordinate_name = coordinate_name, set_attributes = set_attributes)
@deprecate starsubdivision(PF::_FanLikeType, n::Int) star_subdivision(PF, n)

# Deprecated after 0.12.1
@deprecate hirzebruch_surface(r::Int; set_attributes::Bool = true) hirzebruch_surface(NormalToricVariety, r; set_attributes = set_attributes)
@deprecate del_pezzo_surface(b::Int; set_attributes::Bool = true) del_pezzo_surface(NormalToricVariety, b; set_attributes = set_attributes)

# PolyhedralComplex -> polyhedral_complex
function PolyhedralComplex{T}(
                polyhedra::IncidenceMatrix, 
                vr::AbstractCollection[PointVector], 
                far_vertices::Union{Vector{Int}, Nothing} = nothing, 
                L::Union{AbstractCollection[RayVector], Nothing} = nothing;
                non_redundant::Bool = false
            ) where T<:scalar_types
  Base.depwarn("'PolyhedralComplex{$T}(x...)' is deprecated, use 'polyhedral_complex($T, x...)' instead.", :PolyhedralComplex)
  return polyhedral_complex(T, polyhedra, vr, far_vertices, L; non_redundant=non_redundant)
end
function PolyhedralComplex(
                polyhedra::IncidenceMatrix, 
                vr::AbstractCollection[PointVector], 
                far_vertices::Union{Vector{Int}, Nothing} = nothing, 
                L::Union{AbstractCollection[RayVector], Nothing} = nothing;
                non_redundant::Bool = false
            )
  Base.depwarn("'PolyhedralComplex' is deprecated, use 'polyhedral_complex' instead.", :PolyhedralComplex)
  return polyhedral_complex(QQFieldElem, polyhedra, vr, far_vertices, L; non_redundant=non_redundant)
end
function PolyhedralComplex(iter::SubObjectIterator{Polyhedron{T}}) where T<:scalar_types
  Base.depwarn("'PolyhedralComplex' is deprecated, use 'polyhedral_complex' instead.", :PolyhedralComplex)
  return polyhedral_complex(iter)
end
@deprecate PolyhedralComplex(p::Polymake.BigObject) polyhedral_complex(p)

# PolyhedralFan -> polyhedral_fan
@deprecate PolyhedralFan(p::Polymake.BigObject) polyhedral_fan(p)
function PolyhedralFan{T}(Rays::AbstractCollection[RayVector], 
                          LS::Union{AbstractCollection[RayVector], Nothing}, 
                          Incidence::IncidenceMatrix; 
                          non_redundant::Bool = false) where T<:scalar_types
  Base.depwarn("'PolyhedralFan{$T}(x...)' is deprecated, use 'polyhedral_fan($T, x...)' instead.", :PolyhedralFan)
  return polyhedral_fan(T, Rays, LS, Incidence; non_redundant=non_redundant)
end
function PolyhedralFan{T}(Rays::AbstractCollection[RayVector], Incidence::IncidenceMatrix; non_redundant::Bool = false) where T<:scalar_types
  Base.depwarn("'PolyhedralFan{$T}(x...)' is deprecated, use 'polyhedral_fan($T, x...)' instead.", :PolyhedralFan)
  return polyhedral_fan(T, Rays, Incidence; non_redundant=non_redundant)
end
@deprecate PolyhedralFan(Rays::AbstractCollection[RayVector], LS::Union{AbstractCollection[RayVector], Nothing}, Incidence::IncidenceMatrix; non_redundant::Bool = false) polyhedral_fan(QQFieldElem, Rays, LS, Incidence; non_redundant = non_redundant)
@deprecate PolyhedralFan(Rays::AbstractCollection[RayVector], Incidence::IncidenceMatrix; non_redundant::Bool = false) polyhedral_fan(QQFieldElem, Rays, Incidence; non_redundant = non_redundant)
function PolyhedralFan(itr::AbstractVector{Cone{T}}) where T<:scalar_types
  Base.depwarn("'PolyhedralFan' is deprecated, use 'polyhedral_fan' instead.", :PolyhedralFan)
  return polyhedral_fan(itr)
end
function PolyhedralFan(C::Cone{T}) where T<:scalar_types
  Base.depwarn("'PolyhedralFan' is deprecated, use 'polyhedral_fan' instead.", :PolyhedralFan)
  return polyhedral_fan(C)
end
function PolyhedralFan{T}(Rays::AbstractCollection[RayVector], LS::AbstractCollection[RayVector], Incidence::Matrix{Bool}) where T<:scalar_types
  Base.depwarn("'PolyhedralFan{$T}(x...)' is deprecated, use 'polyhedral_fan($T, x...)' instead.", :PolyhedralFan)
  return polyhedral_fan(T, Rays, LS, Incidence)
end
function PolyhedralFan{T}(Rays::AbstractCollection[RayVector], Incidence::Matrix{Bool}) where T<:scalar_types
  Base.depwarn("'PolyhedralFan{$T}(x...)' is deprecated, use 'polyhedral_fan($T, x...)' instead.", :PolyhedralFan)
  return polyhedral_fan(T, Rays, Incidence)
end

# SubdivisionOfPoints -> subdivision_of_points
@deprecate SubdivisionOfPoints(points, C) subdivision_of_points(points, C)
function SubdivisionOfPoints{T}(points, C) where T<:scalar_types
  Base.depwarn("'SubdivisionOfPoints{$T}(x...)' is deprecated, use 'subdivision_of_points($T, x...)' instead.", :SubdivisionOfPoints)
  return subdivision_of_points(T, points, C)
end
@deprecate SubdivisionOfPoints(p::Polymake.BigObject) subdivision_of_points(p)

# Polyhedron -> polyhedron
function Polyhedron{T}(first, second) where T<:scalar_types
  Base.depwarn("'Polyhedron{$T}(x...)' is deprecated, use 'polyhedron($T, x...)' instead.", :Polyhedron)
  return polyhedron(T, first, second)
end
@deprecate Polyhedron(A) polyhedron(A)
@deprecate Polyhedron(A, b) polyhedron(A, b)

# Cone -> positive_hull
@deprecate Cone(R;kwargs...) positive_hull(R;kwargs...)
@deprecate Cone(R,L;kwargs...) positive_hull(R,L;kwargs...)
function Cone{T}(R, L; kwargs...) where T<:scalar_types
  Base.depwarn("'Cone{$T}(x...)' is deprecated, use 'positive_hull($T, x...)' instead.", :Polyhedron)
  return positive_hull(T, R, L; kwargs...)
end
function Cone{T}(R; kwargs...) where T<:scalar_types
  Base.depwarn("'Cone{$T}(x...)' is deprecated, use 'positive_hull($T, x...)' instead.", :Polyhedron)
  return positive_hull(T, R; kwargs...)
end

# LinearProgram -> linear_program
function LinearProgram(p::Polyhedron{T}, lp, c) where T<:scalar_types
  Base.depwarn("'LinearProgram(x...)' is deprecated, use 'linear_program(x...)' instead.", :LinearProgram)
  return linear_program(p, lp, c)
end
function LinearProgram{T}(P::Polyhedron{T}, objective::AbstractVector; kwargs...) where T<:scalar_types
  Base.depwarn("'LinearProgram(x...)' is deprecated, use 'linear_program(x...)' instead.", :LinearProgram)
  return linear_program(P, objective; kwargs...)
end
function LinearProgram(P::Polyhedron{T}, objective::AbstractVector; kwargs...) where T<:scalar_types
  Base.depwarn("'LinearProgram(x...)' is deprecated, use 'linear_program(x...)' instead.", :LinearProgram)
  return linear_program(P, objective; kwargs...)
end
function LinearProgram{T}(A::Union{Oscar.MatElem,AbstractMatrix}, b, c::AbstractVector; kwargs...)  where T<:scalar_types
  Base.depwarn("'LinearProgram(x...)' is deprecated, use 'linear_program(x...)' instead.", :LinearProgram)
  return linear_program(T, A, b, c; kwargs...)
end


# MixedIntegerLinearProgram -> mixed_integer_linear_program
function MixedIntegerLinearProgram(p::Polyhedron{T}, lp, c) where T<:scalar_types
  Base.depwarn("'MixedIntegerLinearProgram(x...)' is deprecated, use 'mixed_integer_linear_program(x...)' instead.", :MixedIntegerLinearProgram)
  return mixed_integer_linear_program(p, lp, c)
end
function MixedIntegerLinearProgram{T}(P::Polyhedron{T}, objective::AbstractVector; kwargs...) where T<:scalar_types
  Base.depwarn("'MixedIntegerLinearProgram(x...)' is deprecated, use 'mixed_integer_linear_program(x...)' instead.", :MixedIntegerLinearProgram)
  return mixed_integer_linear_program(P, objective; kwargs...)
end
function MixedIntegerLinearProgram(P::Polyhedron{T}, objective::AbstractVector; kwargs...) where T<:scalar_types
  Base.depwarn("'MixedIntegerLinearProgram(x...)' is deprecated, use 'mixed_integer_linear_program(x...)' instead.", :MixedIntegerLinearProgram)
  return mixed_integer_linear_program(P, objective; kwargs...)
end
function MixedIntegerLinearProgram{T}(A::Union{Oscar.MatElem,AbstractMatrix}, b, c::AbstractVector; kwargs...)  where T<:scalar_types
  Base.depwarn("'MixedIntegerLinearProgram(x...)' is deprecated, use 'mixed_integer_linear_program(x...)' instead.", :MixedIntegerLinearProgram)
  return mixed_integer_linear_program(T, A, b, c; kwargs...)
end

# see https://github.com/oscar-system/Oscar.jl/pull/2368
@deprecate FreeModElem(coords::SRow{T}, parent::FreeMod_dec{T}) where T <: CRingElem_dec FreeModElem_dec(coords, parent)

# see https://github.com/oscar-system/Oscar.jl/pull/2519
@deprecate group_class_function(tbl::GAPGroupCharacterTable, values::GapObj) class_function(tbl, values)
@deprecate group_class_function(tbl::GAPGroupCharacterTable, values::Vector{<:QQAbElem}) class_function(tbl, values)
@deprecate group_class_function(G::GAPGroup, values::GapObj) class_function(G, values)
@deprecate group_class_function(G::GAPGroup, values::Vector{<:QQAbElem}) class_function(G, values)

# Deprecated after 0.13.0
@deprecate fan(v::NormalToricVarietyType) polyhedral_fan(v)
@deprecate restrict_automorphism_group(G::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMor, check::Bool) restrict_automorphism_group(G, i; check)

# Polyhedral object wrappers now require a parent field
function Cone{T}(obj::Polymake.BigObject) where T<:scalar_types
  Base.depwarn("'Cone{$T}(obj::Polymake.BigObject)' is deprecated, use 'Cone{$T}(obj, f)' with 'f::Field', or 'cone(obj)' instead.", :Cone)
  return Cone{T}(obj, _detect_default_field(T, obj))
end

function PolyhedralComplex{T}(obj::Polymake.BigObject) where T<:scalar_types
  Base.depwarn("'PolyhedralComplex{$T}(obj::Polymake.BigObject)' is deprecated, use 'PolyhedralComplex{$T}(obj, f)' with 'f::Field', or 'polyhedral_complex(obj)' instead.", :PolyhedralComplex)
  return PolyhedralComplex{T}(obj, _detect_default_field(T, obj))
end

function PolyhedralFan{T}(obj::Polymake.BigObject) where T<:scalar_types
  Base.depwarn("'PolyhedralFan{$T}(obj::Polymake.BigObject)' is deprecated, use 'PolyhedralFan{$T}(obj, f)' with 'f::Field', or 'polyhedral_fan(obj)' instead.", :PolyhedralFan)
  return PolyhedralFan{T}(obj, _detect_default_field(T, obj))
end

function Polyhedron{T}(obj::Polymake.BigObject) where T<:scalar_types
  Base.depwarn("'Polyhedron{$T}(obj::Polymake.BigObject)' is deprecated, use 'Polyhedron{$T}(obj, f)' with 'f::Field', or 'polyhedron(obj)' instead.", :Polyhedron)
  return Polyhedron{T}(obj, _detect_default_field(T, obj))
end

function SubdivisionOfPoints{T}(obj::Polymake.BigObject) where T<:scalar_types
  Base.depwarn("'SubdivisionOfPoints{$T}(obj::Polymake.BigObject)' is deprecated, use 'SubdivisionOfPoints{$T}(obj, f)' with 'f::Field', or 'subdivision_of_points(obj)' instead.", :SubdivisionOfPoints)
  return SubdivisionOfPoints{T}(obj, _detect_default_field(T, obj))
end

@deprecate is_hermitian_matrix(B::MatElem{T}) where T <: FinFieldElem is_hermitian(B)
@alias ishermitian_matrix is_hermitian_matrix

@deprecate is_alternating_form(f::SesquilinearForm) is_alternating(f)
@deprecate is_hermitian_form(f::SesquilinearForm) is_hermitian(f)
@deprecate is_quadratic_form(f::SesquilinearForm) is_quadratic(f)
@deprecate is_symmetric_form(f::SesquilinearForm) is_symmetric(f)

@alias isalternating_form is_alternating_form
@alias ishermitian_form is_hermitian_form
@alias isquadratic_form is_quadratic_form
@alias issymmetric_form is_symmetric_form

@deprecate dual_cone(v::AffineNormalToricVariety) weight_cone(v)
@deprecate cohomology_index(tvs::ToricVanishingSet) cohomology_indices(tvs)
