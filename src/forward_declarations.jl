@attributes mutable struct NormalToricVariety
           polymakeNTV::Polymake.BigObject
           NormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct AffineNormalToricVariety
           polymakeNTV::Polymake.BigObject
           AffineNormalToricVariety(polymakeNTV::Polymake.BigObject) = new(polymakeNTV)
end

@attributes mutable struct CyclicQuotientSingularity
    polymakeNTV::Polymake.BigObject
end

NormalToricVarietyType = Union{NormalToricVariety, AffineNormalToricVariety, CyclicQuotientSingularity}
