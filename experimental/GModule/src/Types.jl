mutable struct GModuleHom{ G, T1, T2} <: Map{GModule{G, T1}, GModule{G, T2}, OscarMap, GModuleHom}

  GM1::GModule{G, T1}
  GM2::GModule{G, T2}
  module_map::Map

  function GModuleHom(
    M1::GModule,
    M2::GModule,
    mp::Map;
    check::Bool = false
    )
    # Need to require that
    #   1. Both GModules have the same group
    #   2. The group action is respected
    @req M1.G === M2.G "groups need to be identical"
    @req domain(mp) === M1.M && codomain(mp) === M2.M "map need to map 1st module into 2nd"
    #not every hom is a G-Hom...that is what check is supposed to do - eventually
    #see 2.
    if check #only works if mp is a morphism so that "*" and "==" are doing
      #s.th. useful
      @assert all(g->action(M1, g)*mp == mp*action(M2, g), gens(M1.G))
    end

    return new{typeof(M1.G), typeof(M1.M), typeof(M2.M)}(M1, M2, mp)
  end
end

function hom(M1::GModule{T}, M2::GModule{T}, mp::Map; check::Bool = true) where T <: AbstractAlgebra.Group
  return GModuleHom(M1, M2, mp; check)
end

function hom(M1::GModule{T}, M2::GModule{T}, mp::MatElem; check::Bool = true) where T <: AbstractAlgebra.Group
  return GModuleHom(M1, M2, hom(M1.M, M2.M, mp); check)
end

domain(M::GModuleHom) = M.GM1
codomain(M::GModuleHom) = M.GM2
parent(M::GModuleHom) = Hecke.MapParent(domain(M), codomain(M), "homomorphisms")

mutable struct GModuleElem{T}
  parent::GModule
  data::T
end

parent(a::GModuleElem) = a.parent

function (C::GModule)(a::Union{ModuleElem, FinGenAbGroupElem})
  @req parent(a) === C.M "wrong parent for $a"
  return GModuleElem(C, a)
end

function ==(a::GModuleElem, b::GModuleElem)
  @req parent(a) === parent(b) "parents differ"
  return a.data == b.data
end

function hash(a::GModuleElem, u::UInt)
  return hash(a.data, u)
end

function +(a::GModuleElem, b::GModuleElem)
  @req parent(a) === parent(b) "parents differ"
  return GModuleElem(parent(a), a.data + b.data)
end

function -(a::GModuleElem, b::GModuleElem)
  @req parent(a) === parent(b) "parents differ"
  return GModuleElem(parent(a), a.data - b.data)
end

function -(a::GModuleElem)
  return GModuleElem(parent(a), -a.data)
end

function *(a::GModuleElem, g::GroupElem)
  @req parent(a).G === parent(g) "group element has wrong parent"
  return GModuleElem(parent(a), action(parent(a), g, a.data))
end

function (A::GModuleHom)(a::GModuleElem)
  @req parent(a) === domain(A) "element has wrong parent"
  return GModuleElem(codomain(A), A.module_map(a.data))
end

function kernel(A::GModuleHom)
  return sub(domain(A), kernel(A.module_map)[2])
end

function image(A::GModuleHom)
  return sub(codomain(A), image(A.module_map)[2])
end
