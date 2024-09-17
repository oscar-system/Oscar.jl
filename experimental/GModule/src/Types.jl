mutable struct GModuleHom{
    G<:Any,
    T1<:Any,
    T2<:Any,
    RingMapType<:Any} <: Map{GModule{G, T1}, GModule{G, T2}}

    GM1::GModule{G, T1}
    GM2::GModule{G, T2}
    module_map::Map{T1, T2}

    function GModuleHom(
        M1::GModule{G, AbstractFreeMod},
        M2::GModule{G, S},
        a::Vector{ModuleElemType}
        ) where {G, S<:ModuleFP, ModuleElemtype<:ModuleFPElem}
        # Need to require that
        #   1. Both GModules have the same group
        #   2. The group action is respected
        @assert M1.gT == M2.gT
        r = new{G, typeof(M1.mT), S, Nothing}(M1, M2, FreeModuleHom(M1.mT, M2.mT, a))
    end

    function GModuleHom(
        M1::GModule{G, AbstractFreeMod},
        M2::GModule{G, S},
        a::Vector{ModuleElemType},
        h::RingMapType
      ) where {G, S, ModuleElemType<:ModuleFPElem, RingMapType}
      @assert M1.gT == M2.gT
      r = new{G, typeof(M1.mT), S, RingMapType}(M1, M2, FreeModuleHom(M1.mT, M2.mT, a, h))
    end
end


function GModuleHom(
    M1::GModule{G, AbstractFreeMod{T}},
    M2::GModule{G, S},
    mat::MatElem{T}) where {T<:RingElem, S<:AbstractFreeMod}
    @assert M1.gT == M2.gT
    r = new{G, typeof(M1.mT), S, Nothing}(M1, M2, FreeModuleHom(M1.mT, m2.mT, mat))
end

function GModuleHom(
    M1::GModule{G, AbstractFreeMod{T}},
    M2::GModule{G, S},
    mat::MatElem{T}) where {T<:RingElem, S<:ModuleFP}
    @assert M1.gT == M2.gT
    r = new{G, typeof(M1.mT), S, Nothing}(M1, M2, FreeModuleHom(M1.mT, m2.mT, mat))
end

function GModuleHom(
    M1::GModule{G, AbstractFreeMod{T}},
    M2::GModule{G, S},
    mat::MatElem{T}, h::RingMapType) where {T<:RingElem, S<:AbstractFreeMod, RingMapType}
    @assert M1.gT == M2.gT
    r = new{G, typeof(M1.mT), S, Nothing}(M1, M2, FreeModuleHom(M1.mT, m2.mT, mat, h))
end

function GModuleHom(
    M1::GModule{G, AbstractFreeMod{T}},
    M2::GModule{G, S},
    mat::MatElem{T}, h::RingMapType) where {T<:RingElem, S<:ModuleFP, RingMapType}
    @assert M1.gT == M2.gT
    r = new{G, typeof(M1.mT), S, Nothing}(M1, M2, FreeModuleHom(M1.mT, m2.mT, mat, h))
end


domain(M::Map(GModuleHom)) = M.Gm1
codomain(M::Map(GModuleHom)) = M.Gm2
