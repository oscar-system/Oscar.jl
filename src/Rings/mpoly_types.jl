export mpoly_type, mpoly_ring_type

mpoly_ring_type(R::T) where {T<:AbstractAlgebra.Ring} = MPolyRing{elem_type(T)}
mpoly_ring_type(::Type{T}) where {T<:AbstractAlgebra.Ring} = MPolyRing{elem_type(T)}
mpoly_type(R::T) where {T<:AbstractAlgebra.Ring} = MPolyElem{elem_type(T)}
mpoly_type(::Type{T}) where {T<:AbstractAlgebra.Ring} = MPolyElem{elem_type(T)}

mpoly_ring_type(R::T) where {T<:AbstractAlgebra.Field} = MPolyRing{elem_type(T)}
mpoly_ring_type(::Type{T}) where {T<:AbstractAlgebra.Field} = MPolyRing{elem_type(T)}
mpoly_ring_type(a::T) where {T<:AbstractAlgebra.FieldElem} = MPolyRing{T}
mpoly_ring_type(::Type{T}) where {T<:AbstractAlgebra.FieldElem} = MPolyRing{T}

mpoly_type(R::T) where {T<:AbstractAlgebra.Field} = MPolyElem{elem_type(T)}
mpoly_type(::Type{T}) where {T<:AbstractAlgebra.Field} = MPolyElem{elem_type(T)}
mpoly_type(a::T) where {T<:AbstractAlgebra.FieldElem} = MPolyElem{T}
mpoly_type(::Type{T}) where {T<:AbstractAlgebra.FieldElem} = MPolyElem{T}

mpoly_ring_type(R::FlintRationalField) = FmpqMPolyRing
mpoly_ring_type(::Type{FlintRationalField}) = FmpqMPolyRing
mpoly_ring_type(a::fmpq) = FmpqMPolyRing
mpoly_ring_type(::Type{fmpq}) = FmpqMPolyRing

mpoly_type(R::FlintRationalField) = fmpq_mpoly
mpoly_type(::Type{FlintRationalField}) = fmpq_mpoly
mpoly_type(a::fmpq) = fmpq_mpoly
mpoly_type(::Type{fmpq}) = fmpq_mpoly

mpoly_ring_type(R::Nemo.GaloisField) = GFPMPolyRing
mpoly_ring_type(::Type{Nemo.GaloisField}) = GFPMPolyRing
mpoly_ring_type(a::gfp_elem) = GFPMPolyRing
mpoly_ring_type(::Type{gfp_elem}) = GFPMPolyRing

mpoly_type(R::Nemo.GaloisField) = gfp_mpoly
mpoly_type(::Type{Nemo.GaloisField}) = gfp_mpoly
mpoly_type(a::gfp_elem) = gfp_mpoly
mpoly_type(::Type{gfp_elem}) = gfp_mpoly
