#=
# Since some code coverage tool blocks the tests from passing if this file is 
# present with proper julia code inside, the code below is commented out. 
#
# If you want to implement your own hyper complex class, you may start with 
# the template below and replace every occurrence of `SimplicialCo` with your favourite
# name for your new class. Then you have to fill in the gaps according to your 
# needs. We provide a sample implementation below.
=#


### Production of the chains
struct SimplicialCoChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  R::Ring
  K::SimplicialComplex

  function SimplicialCoChainFactory(R::Ring, K::SimplicialComplex)
    return new{FreeMod{elem_type(R)}}(R, K)
  end
end

function (fac::SimplicialCoChainFactory)(self::AbsHyperComplex, i::Tuple)
  return FreeMod(fac.R, length(faces(fac.K, only(i))))
end

function can_compute(fac::SimplicialCoChainFactory, self::AbsHyperComplex, i::Tuple)
  return 0 <= only(i) <= dim(fac.K)
end

### Production of the morphisms 
struct SimplicialCoMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  function SimplicialCoMapFactory()
    return new{FreeModuleHom}()
  end
end

function (::SimplicialCoMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  fac = chain_factory(self)
  i = only(I)
  dom = self[i]
  cod = self[i+1]
  return hom(dom, cod, transpose(matrix(fac.R, Polymake.topaz.boundary_matrix(fac.K.pm_simplicialcomplex, i+1))))
end

function can_compute(fac::SimplicialCoMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  fac = chain_factory(self)
  p == 1 || return false
  return 0 <= only(i) < dim(fac.K)
end

### The concrete struct
@attributes mutable struct SimplicialCoComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function SimplicialCoComplex(R::Ring, K::SimplicialComplex)
    chain_fac = SimplicialCoChainFactory(R::Ring, K::SimplicialComplex)
    map_fac = SimplicialCoMapFactory()

    internal_complex = HyperComplex(1, chain_fac, map_fac, [:cochain])
    ChainType = FreeMod{elem_type(R)}
    MorphismType = FreeModuleHom
    return new{ChainType, MorphismType}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::SimplicialCoComplex) = c.internal_complex

