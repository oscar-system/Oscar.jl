#=
# Since some code coverage tool blocks the tests from passing if this file is 
# present with proper julia code inside, the code below is commented out. 
#
# If you want to implement your own hyper complex class, you may start with 
# the template below and replace every occurrence of `MyNew` with your favourite
# name for your new class. Then you have to fill in the gaps according to your 
# needs. We provide a sample implementation below.
=#


#========================== start template  ===================================

### Production of the chains
struct MyNewChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  # Fields needed for production

  function MyNewChainFactory(...)
    # Fill in the constructor
  end
end

function (fac::MyNewChainFactory)(self::AbsHyperComplex, i::Tuple)
  # Production of the chains at index i
end

function can_compute(fac::MyNewChainFactory, self::AbsHyperComplex, i::Tuple)
  # Deciding whether the entry at index i can be produced
end

### Production of the morphisms 
struct MyNewMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  # Fields needed for production

  function MyNewMapFactory(...)
    # Fill in the constructor
  end
end

function (fac::MyNewMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  # Production of the outgoing morphism at index i in the p-th direction
end

function can_compute(fac::MyNewMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  # Deciding whether the outgoing map at index i in the p-th direction can be produced
end

### The concrete struct
@attributes mutable struct MyNewComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}

  function MyNewComplex(...)
    chain_fac = MyNewChainFactory(...)
    map_fac = MyNewMapFactory(...)

    # Assuming d is the dimension of the new complex
    internal_complex = HyperComplex(d, chain_fac, map_fac, [:chain for i in 1:d])
    # Assuming that ChainType and MorphismType are provided by the input
    return new{ChainType, MorphismType}(internal_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::MyNewComplex) = c.internal_complex

========================== end template ===================================#

#=
# A dummy implementation of a "reflected hyper complex".
#
# For a hypercomplex `hc` the `ReflectedComplex(hc)` `rhc` has entry 
# `hc[i]` at index `-i` and, accordingly, all maps in the other direction.
#
=#

### Production of the chains
struct ReflectedChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  # Fields needed for production

  # In this case we only need the original complex whose reflected complex we are building.
  original_complex::AbsHyperComplex

  function ReflectedChainFactory(original_complex::AbsHyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType}
    return new{ChainType}(original_complex)
  end
end

function (fac::ReflectedChainFactory)(self::AbsHyperComplex, i::Tuple)
  # `self` is always the instance of the concrete hyper complex you are 
  # creating here, `i` is the index to be created. The result will in general
  # be cached (unless the internal HyperComplex below is created with `cached=false`)
  # and this function will only be called at most once for every index `i`.

  # Create the reflected index
  j = Tuple(-collect(i))
  # Access the required data for production from the internals of the factory:
  return fac.original_complex[j] 
end

function can_compute(fac::ReflectedChainFactory, self::AbsHyperComplex, i::Tuple)
  # Determine whether or not this index can be computed. This is used to inform 
  # the outside program as to whether or not asking for the `i`-th index is 
  # legitimate.

  # Create the reflected index
  j = Tuple(-collect(i))
  # Access the required data for production from the internals of the factory:
  return can_compute_index(fac.original_complex, j)
end

### Production of the morphisms 
struct ReflectedMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  # Fields needed for production

  # In this case we only need the original complex whose reflected complex we are building.
  original_complex::AbsHyperComplex

  function ReflectedMapFactory(original_complex::AbsHyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType}
    return new{MorphismType}(original_complex)
  end
end

function (fac::ReflectedMapFactory)(self::AbsHyperComplex, p::Int, i::Tuple)
  # Production of the outgoing morphism at index i in the p-th direction

  # Create the reflected index
  j = Tuple(-collect(i))

  # Return the corresponding map from the original complex; 
  # access to the relevant data from within the factory
  return map(fac.original_complex, p, j)
end

function can_compute(fac::ReflectedMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  # Decide whether the outgoing map at index i in the p-th direction can be produced;
  # see `can_compute` above.
  
  # Create the reflected index
  j = Tuple(-collect(i))
  return can_compute_map(fac.original_complex, p, j)
end

### The concrete struct
# 
# It is advised to use another concrete instance of `HyperComplex` in the internals 
# and forward its functionality via `underlying_complex` below. This will take care 
# of the caching of elements and implementing the whole interface for hyper complexes 
# to the outside using only your factories above. 
@attributes mutable struct ReflectedComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  # Internal data for your new complex
  internal_complex::HyperComplex{ChainType, MorphismType}

  # There could be more. In this case, a reflected complex should 
  # know that there is some original complex where it comes from.
  original_complex::AbsHyperComplex{ChainType, MorphismType}

  function ReflectedComplex(original_complex::AbsHyperComplex{ChainType, MorphismType}) where {ChainType, MorphismType}

    # Create the factories as they are required as arguments for the `HyperComplex` constructor.
    chain_fac = ReflectedChainFactory(original_complex)
    map_fac = ReflectedMapFactory(original_complex)

    # directions are inverted by reflection
    dirs = [(direction(original_complex, i) == :chain ? (:cochain) : (:chain)) for i in 1:dim(original_complex)]
    upper_bounds = Union{Int, Nothing}[has_lower_bound(original_complex, i) ? -lower_bound(original_complex, i) : nothing for i in 1:dim(original_complex)]
    lower_bounds = Union{Int, Nothing}[has_upper_bound(original_complex, i) ? -upper_bound(original_complex, i) : nothing for i in 1:dim(original_complex)]

    # Create an instance of `HyperComplex` using your factories.
    #
    # Actually, this hyper complex could have been created with `cached=false` below, 
    # because caching is already happening in the `original_complex`. However, this 
    # is not advised in general!
    internal_complex = HyperComplex(
                                    dim(original_complex), # the reflected complex has the same dimension
                                    chain_fac, map_fac,    # use your factories
                                    dirs;
                                    upper_bounds,
                                    lower_bounds #, cached=true # an optional argument to switch off caching 
                                                        # and only use the factories.
                                                        # Use only if you know what you're doing!
                                   )
    # `ChainType` and `MorphismType` are the same as the original complex.
    return new{ChainType, MorphismType}(internal_complex, original_complex)
  end
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::ReflectedComplex) = c.internal_complex

### Implement specialized functionality for your new instance of ReflectedComplex
original_complex(c::ReflectedComplex) = c.original_complex

