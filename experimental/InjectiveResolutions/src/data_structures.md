# Data structures for injective resolutions and local cohomology modules

## monoid algebras
> struct FaceQ # face of semigroup  
    prime::Union{MPolyIdeal,MPolyQuoIdeal} 
    poly::Polyhedron  #face of Q corresponding to prime
> end

> struct HyperplaneQ # a hyperplane bounding the cone RR_{\geq 0}Q
    hyperplane::Polyhedron
    A::Matrix{Int}
    b::Vector{Int}
> end

> struct MonoidAlgebra # monoid algebra with associated data
    algebra::Union{MPolyRing, MPolyQuoRing}
    cone::Polyhedron
    faces::Vector{FaceQ}
    hyperplanes::Vector{HyperplaneQ}
    # normal::Bool
    pointed::Bool #cone pointed
    zonotope::Tuple{Polyhedron,Vector{Int}}
> end

> struct MonoidAlgebraIdeal
    monoidAlgebra::MonoidAlgebra 
    ideal::Ideal
> end


### modules over monoid algebras
> struct MonoidAlgebraModule
    monoidAlgebra::MonoidAlgebra
    mod::SubquoModule
> end

> struct InjMod
    monoidAlgebra::MonoidAlgebra
    indecInjectives::Vector{IndecInj}
> end

## injective and irreducible resolutions
#### changed!
> struct IndecInj # indecomposable injective
    <!-- prime::Ideal -->
    <!-- face::Polyhedron -->
    face::FaceQ
    vector::Vector{Int}
> end

> struct IrrSum # irreducible sum
    mod::SubquoModule
    components::Vector{IndecInj}
> end

> struct IrrRes # irreducible resolution (including all computed data and the cochain complex)
    irrSums::Vector{IrrSum}
    cochainMaps::Vector{SubQuoHom}
    surjections::Vector{SubQuoHom}
    inclusions::Vector{SubQuoHom}
    cokernels::Vector{SubquoModule}
    cochainComplex::ComplexOfMorphisms # if sequence not exact return trivial cochain_complex (M0 -> M0)
> end

> struct InjRes
    injMods::Vector{InjMod}
    cochainMaps::Vector{MatElem}
    upto::Int
    irrRes::IrrRes
    <!-- shift::Vector{Int} -->
> end

**not needed**
> struct IrrHull # irreducible hull
    faces::Vector{FaceQ}
    vectors::Vector{Tuple{SubquoModuleElem,Vector{Int}}}
    lambda::AbstractAlgebra.Generic.MatSpaceElem
> end

## local cohomology
> struct Sector 
    A::Vector{Vector{Int}}
    sector::Polyhedron
    indexVector::Vector{Int}
    J_A::FreeModule
    H::QuotientModule
> end 

> struct SectorPartition
    sectors::Vector{Sectors}
    maps::Vector{Tuple{Sector,Sector, ModuleHomomorphism}}
> end

> 
