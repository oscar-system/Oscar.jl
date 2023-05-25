@attributes mutable struct NormalToricVariety{BRT} <: AbsCoveredScheme{BRT}
           polymakeNTV::Polymake.BigObject
           NormalToricVariety(polymakeNTV::Polymake.BigObject) = new{QQFieldElem}(polymakeNTV)
end

@attributes mutable struct AffineNormalToricVariety{BRT, RT} <: AbsSpec{BRT, RT}
           polymakeNTV::Polymake.BigObject
           AffineNormalToricVariety(polymakeNTV::Polymake.BigObject) = new{QQFieldElem, MPolyQuo}(polymakeNTV)
end

@attributes mutable struct CyclicQuotientSingularity <: AbsCoveredScheme{BRT}
    polymakeNTV::Polymake.BigObject
    CyclicQuotientSingularity(polymakeNTV::Polymake.BigObject) = new{QQFieldElem}(polymakeNTV)
end

NormalToricVarietyType = Union{NormalToricVariety, AffineNormalToricVariety, CyclicQuotientSingularity}
