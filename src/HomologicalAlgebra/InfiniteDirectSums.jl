export InfiniteDirectSum, index_set
export InfiniteDirectSumElem, parent, support, components
export InfiniteDirectSumHom, domain, codomain, map_on_generators

mutable struct InfiniteDirectSum{RingType, RingElemType, IndexSetType, IndexType}
  base_ring::RingType
  index_set::IndexSetType

  #fields for caching
  cached_keys::Vector{IndexType}

  function InfiniteDirectSum(R::RingType, index_set::IndexSetType) where {RingType, IndexSetType}
    return new{RingType, elem_type(R), IndexSetType, elem_type(index_set)}(R, index_set)
  end
end

base_ring(F::InfiniteDirectSum) = F.base_ring
index_set(F::InfiniteDirectSum) = F.index_set

mutable struct InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType}
  parent::InfiniteDirectSum{RingType, RingElemType, IndexSetType, IndexType}
  components::Dict{IndexType, RingElemType}

  function InfiniteDirectSumElem(
      F::InfiniteDirectSum{RingType, RingElemType, IndexSetType, IndexType}, 
      v::Dict{IndexType, RingElemType}
    ) where {RingType, RingElemType, IndexSetType, IndexType}
    for k in keys(v)
      parent(v[k]) == base_ring(F) || error("elements do not belong to the correct ring")
    end
    return new{RingType, RingElemType, IndexSetType, IndexType}(F, v)
  end
end

parent(v::InfiniteDirectSumElem) = v.parent
support(v::InfiniteDirectSumElem) = keys(v.components)
components(v::InfiniteDirectSumElem) = v.components
base_ring(v::InfiniteDirectSumElem) = base_ring(parent(v))

function getindex(
    v::InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType}, 
    i::IndexType
  ) where {RingType, RingElemType, IndexSetType, IndexType}
  haskey(components(v), i) && return components(v)[i]
  return zero(base_ring(v))
end

function getindex(
    v::InfiniteDirectSumElem,
    i::Integer
  ) 
  return v[index_set(v)(i)]
end


function getindex(
    F::InfiniteDirectSum{RingType, RingElemType, IndexSetType, IndexType}, 
    i::IndexType
  ) where {RingType, RingElemType, IndexSetType, IndexType}
  D = Dict{IndexType, RingElemType}()
  D[i] = one(base_ring(F))
  return InfiniteDirectSumElem(F, D)
end

function getindex(
    F::InfiniteDirectSum,
    i::Integer
  ) 
  return F[index_set(F)(i)]
end

zero(F::InfiniteDirectSum) = InfiniteDirectSumElem(F, Dict{elem_type(index_set(F)), elem_type(base_ring(F))}())

function +(
	   u::InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType},
	   v::InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType}
	   ) where {RingType, RingElemType, IndexSetType, IndexType}
  parent(u) == parent(v) || error("elements do not belong to the same module")
  values = Dict{IndexType, RingElemType}()
  for k in union(support(u), support(v))
    values[k] = u[k] + v[k]
  end
  return InfiniteDirectSumElem(parent(u), values)
end

function *(
	   a::RingElemType,
	   u::InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType}
	   ) where {RingType, RingElemType, IndexSetType, IndexType}
  values = Dict{IndexType, RingElemType}()
  for i in support(u)
    values[i] = a*u[i]
  end
  return InfiniteDirectSumElem(parent(u), values)
end

function *(
	   a::Integer,
	   u::InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType}
	   ) where {RingType, RingElemType, IndexSetType, IndexType}
  return base_ring(u)(a)*u
end

function add!(
    u::InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType},
    v::InfiniteDirectSumElem{RingType, RingElemType, IndexSetType, IndexType}
  ) where {RingType, RingElemType, IndexSetType, IndexType}
  parent(u) == parent(v) || error("elements do not belong to the same module")
  for k in support(v)
    components(u)[k] = u[k] + v[k]
  end
  return u
end
  
mutable struct InfiniteDirectSumHom{RingType, RingElemType, DomainIndexSetType, DomainIndexType, CodomainIndexSetType, CodomainIndexType}
  domain::InfiniteDirectSum{RingType, RingElemType, DomainIndexSetType, DomainIndexType}
  codomain::InfiniteDirectSum{RingType, RingElemType, CodomainIndexSetType, CodomainIndexType}
  map_on_generators::Function

  function InfiniteDirectSumHom(
      domain::InfiniteDirectSum{RingType, RingElemType, DomainIndexSetType, DomainIndexType},
      codomain::InfiniteDirectSum{RingType, RingElemType, CodomainIndexSetType, CodomainIndexType},
      map_on_generators::Function
    ) where {RingType, RingElemType, DomainIndexSetType, DomainIndexType, CodomainIndexSetType, CodomainIndexType}
    base_ring(domain) == base_ring(codomain) || error("modules are not defined over the same ring")
    return new{RingType, RingElemType, DomainIndexSetType, DomainIndexType, CodomainIndexSetType, CodomainIndexType}(domain, codomain, map_on_generators)
  end
end

domain(f::InfiniteDirectSumHom) = f.domain
codomain(f::InfiniteDirectSumHom) = f.codomain
map_on_generators(f::InfiniteDirectSumHom) = f.map_on_generators
getindex(f::InfiniteDirectSumHom{RT, RET, DIST, DIT, CIST, CIT}, i::DIT) where {RT, RET, DIST, DIT, CIST, CIT} = map_on_generators(f)(i)
getindex(f::InfiniteDirectSumHom, i::Integer) = f[index_set(domain(f))(i)]

function (f::InfiniteDirectSumHom{RT, RET, DIST, DIT, CIST, CIT})(v::InfiniteDirectSumElem{RT, RET, DIST, DIT}) where {RT, RET, DIST, DIT, CIST, CIT} 
  w = zero(codomain(f))
  for i in support(v)
    add!(w, v[i]*f[i])
  end
  return w
end
    


