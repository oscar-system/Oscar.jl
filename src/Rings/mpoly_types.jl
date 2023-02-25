export mpoly_type, mpoly_ring_type, mpoly_dec_ring_type, mpoly_dec_type

mpoly_ring_type(R::T) where {T<:AbstractAlgebra.Ring} = Generic.MPolyRing{elem_type(T)}
mpoly_ring_type(::Type{T}) where {T<:AbstractAlgebra.Ring} = Generic.MPolyRing{elem_type(T)}
mpoly_type(R::T) where {T<:AbstractAlgebra.Ring} = Generic.MPolyRingElem{elem_type(T)}
mpoly_type(::Type{T}) where {T<:AbstractAlgebra.Ring} = Generic.MPolyRingElem{elem_type(T)}

mpoly_ring_type(a::T) where {T<:AbstractAlgebra.RingElem} = Generic.MPolyRing{T}
mpoly_ring_type(::Type{T}) where {T<:AbstractAlgebra.RingElem} = Generic.MPolyRing{T}

mpoly_type(a::T) where {T<:AbstractAlgebra.RingElem} = Generic.MPolyRingElem{T}
mpoly_type(::Type{T}) where {T<:AbstractAlgebra.RingElem} = Generic.MPolyRingElem{T}

mpoly_ring_type(R::QQField) = QQMPolyRing
mpoly_ring_type(::Type{QQField}) = QQMPolyRing
mpoly_ring_type(a::QQFieldElem) = QQMPolyRing
mpoly_ring_type(::Type{QQFieldElem}) = QQMPolyRing

mpoly_type(R::QQField) = QQMPolyRingElem
mpoly_type(::Type{QQField}) = QQMPolyRingElem
mpoly_type(a::QQFieldElem) = QQMPolyRingElem
mpoly_type(::Type{QQFieldElem}) = QQMPolyRingElem

mpoly_ring_type(R::Nemo.fpField) = fpMPolyRing
mpoly_ring_type(::Type{Nemo.fpField}) = fpMPolyRing
mpoly_ring_type(a::fpFieldElem) = fpMPolyRing
mpoly_ring_type(::Type{fpFieldElem}) = fpMPolyRing

mpoly_type(R::Nemo.fpField) = fpMPolyRingElem
mpoly_type(::Type{Nemo.fpField}) = fpMPolyRingElem
mpoly_type(a::fpFieldElem) = fpMPolyRingElem
mpoly_type(::Type{fpFieldElem}) = fpMPolyRingElem

mpoly_dec_ring_type(A::T) where {T<:MPolyRing} = mpoly_dec_ring_type(typeof(A))
mpoly_dec_type(A::T) where {T<:MPolyRing} = mpoly_dec_type(typeof(A))
mpoly_dec_ring_type(::Type{T}) where {S<:RingElem, T<:MPolyRing{S}} = MPolyDecRing{S, T}
mpoly_dec_type(::Type{T}) where {S<:RingElem, T<:MPolyRing{S}} = MPolyDecRingElem{S, elem_type(T)}
