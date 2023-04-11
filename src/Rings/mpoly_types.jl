mpoly_dec_ring_type(A::T) where {T<:MPolyRing} = mpoly_dec_ring_type(typeof(A))
mpoly_dec_type(A::T) where {T<:MPolyRing} = mpoly_dec_type(typeof(A))
mpoly_dec_ring_type(::Type{T}) where {S<:RingElem, T<:MPolyRing{S}} = MPolyDecRing{S, T}
mpoly_dec_type(::Type{T}) where {S<:RingElem, T<:MPolyRing{S}} = MPolyDecRingElem{S, elem_type(T)}
