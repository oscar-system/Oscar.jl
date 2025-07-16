export ActionPolyRing,
       ActionPolyRingElem,
       DifferencePolyRing,
       DifferencePolyRingElem,
       difference_polynomial_ring,
       ndiffs,
       nelementary_symbols,
       elementary_symbols,
       diff_action

#######################################
#
#  Data types 
#
#######################################

abstract type ActionPolyRing{T} <: Ring end #Perhaps a mutual struct later

abstract type ActionPolyRingElem{T} <: RingElem end

### Difference ###
mutable struct DifferencePolyRing{T} <: ActionPolyRing{T}
  upoly_ring::AbstractAlgebra.Generic.UniversalPolyRing{T}
  elementary_symbols::Vector{Symbol}
  ndiffs::Int
  internal_ordering::Tuple{Symbol, Symbol}
  jet_to_var::Any #Always of type Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}
  var_to_jet::Any #Always of type Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}


  function DifferencePolyRing{T}(R::Ring, nelementary_symbols::Int, ndiffs::Int, internal_ordering::Tuple{Symbol, Symbol}) where {T}
    @req nelementary_symbols >= 0 "The number of elementary symbols must be nonnegative"
    @req ndiffs >= 0 "The number of endomorphisms must be nonnegative"
    @req internal_ordering[1] in [:lex, :deglex, :degrevlex] "ordering of the elementary variables must be one of :lex, :deglex, :degrevlex"
    @req internal_ordering[2] in [:top, :pot] "extension must be one of :top (term-over-position) or :pot (position-over-term)"
    zeroind = fill(0, ndiffs)
    elementary_symbols = map(x -> Symbol("u" * string(x)), 1:nelementary_symbols)
    elem_syms_index = map(x -> Symbol("u" * string(x) * "["* join(zeroind) * "]"), 1:nelementary_symbols) 
    upoly_ring, elemvars = universal_polynomial_ring(R, elem_syms_index)
    
    dpr = new{T}(upoly_ring, elementary_symbols, ndiffs, internal_ordering)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()

    for i in 1:nelementary_symbols
      jet, var = (i, zeroind), dpr(elemvars[i])
      jet_to_var[jet], var_to_jet[var] = var, jet
    end

    dpr.jet_to_var = jet_to_var
    dpr.var_to_jet = var_to_jet
    
    return dpr
  end
 
  function DifferencePolyRing{T}(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int, internal_ordering::Tuple{Symbol, Symbol}) where {T}
    @req ndiffs >= 0 "The number of endomorphisms must be nonnegative"
    @req internal_ordering[1] in [:lex, :deglex, :degrevlex] "ordering of the elementary variables must be one of :lex, :deglex, :degrevlex"
    @req internal_ordering[2] in [:top, :pot] "extension must be one of :top (term-over-position) or :pot (position-over-term)"
    zeroind = fill(0, ndiffs)
    nelementary_symbols = length(elementary_symbols)
    elem_syms_index = map(s -> string(s) * "[" * join(zeroind) * "]", elementary_symbols)
    upoly_ring, elemvars = universal_polynomial_ring(R, elem_syms_index)
    
    dpr = new{T}(upoly_ring, elementary_symbols, ndiffs, internal_ordering)

    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()

    for i in 1:nelementary_symbols
      jet, var = (i, zeroind), dpr(elemvars[i])
      jet_to_var[jet], var_to_jet[var] = var, jet
    end
      
    dpr.jet_to_var = jet_to_var
    dpr.var_to_jet = var_to_jet
    
    return dpr
  end

end

mutable struct DifferencePolyRingElem{T} <: ActionPolyRingElem{T}
  p::AbstractAlgebra.Generic.UniversalPolyRingElem{T}
  parent::DifferencePolyRing{T}
  #leader::AbstractAlgebra.Generic.UniversalPolyRingElem{T}

  DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}) where {T} = new{T}(zero(dpr.upoly_ring), dpr)

  function DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}, upre::AbstractAlgebra.Generic.UniversalPolyRingElem{T}) where {T}
    @req dpr.upoly_ring === parent(upre) "The parent does not match"
    new{T}(upre, dpr)
  end

  function DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}, mpre::MPolyRingElem{T}) where {T}
    @req dpr.upoly_ring.mpoly_ring === parent(mpre) "The parent does not match"
    new{T}(mpre, dpr)
  end

end

elem_type(::Type{DifferencePolyRing{T}}) where {T} = DifferencePolyRingElem{T}

parent_type(::Type{DifferencePolyRingElem{T}}) where {T} = DifferencePolyRing{T}

base_ring_type(::Type{<:ActionPolyRing{T}}) where {T} = parent_type(T)

is_domain_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_domain_type(T)

is_exact_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_exact_type(T)

#######################################
#
#  Construction 
#
#######################################

##### Algebras #####

### Difference ###
@doc raw"""
    difference_polynomial_ring(R::Ring, nelementary_symbols::Int, ndiffs::Int) -> Tuple{DifferencePolyRing, Vector{DifferencePolyRingElem}}

Return a tuple consisting of the difference polynomial ring over the ring `R` with the specified number of elementary variables and commuting endomorphisms, and the vector of
these elementary variables.
"""
function difference_polynomial_ring(R::Ring, nelementary_symbols::Int, ndiffs::Int; internal_ordering = (:lex, :top)::Tuple{Symbol, Symbol})
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, nelementary_symbols, ndiffs, internal_ordering)
  zeroind = fill(0, ndiffs)
  return (dpr, map(i -> dpr.jet_to_var[(i, zeroind)], 1:nelementary_symbols))
end

@doc raw"""
    difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int) -> Tuple{DifferencePolyRing, Vector{DifferencePolyRingElem}}

Return a tuple consisting of the difference polynomial ring over the ring `R` with the specified elementary variables and number of commuting endomorphisms, and the vector of
these elementary variables. Note that the multiindex [0..0] of length 'ndiffs' is appended to the variable names provided.
"""
function difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int; internal_ordering = (:lex, :top)::Tuple{Symbol, Symbol})
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, elementary_symbols, ndiffs, internal_ordering)
  zeroind = fill(0, ndiffs)
  return (dpr, map(i -> dpr.jet_to_var[(i, zeroind)], 1:nelementary_symbols(dpr)))
end

### Differential ###

##### Elements #####

### Difference ###

(dpr::DifferencePolyRing)() = DifferencePolyRingElem{elem_type(base_ring(dpr))}(dpr)

(dpr::DifferencePolyRing)(upre::AbstractAlgebra.Generic.UniversalPolyRingElem) = DifferencePolyRingElem{elem_type(base_ring(dpr))}(dpr, upre)

(dpr::DifferencePolyRing)(a::R) where {R <: RingElement} = dpr(dpr.upoly_ring(a))

### Differential ###

#######################################
#
#  Basic field access 
#
#######################################

##### Algebras #####

### Difference ###
base_ring(dpr::DifferencePolyRing) = base_ring(dpr.upoly_ring)

ndiffs(dpr::DifferencePolyRing) = dpr.ndiffs

elementary_symbols(dpr::DifferencePolyRing) = dpr.elementary_symbols

nelementary_symbols(dpr::DifferencePolyRing) = length(elementary_symbols(dpr))

internal_ordering(dpr::DifferencePolyRing) = dpr.internal_ordering

##### Elements #####

parent(dpre::DifferencePolyRingElem) = dpre.parent

#######################################
#
#  Basic ring functionality 
#
#######################################

zero(dpr::DifferencePolyRing) = dpr(0)

one(dpr::DifferencePolyRing) = dpr(1)

#######################################
#
#  Basic action polynomial functionality 
#
#######################################

@doc raw"""
    is_monomial(p::ActionPolyRingElem)

Return `true` if `p` is a monomial, `false` otherwise.
"""
is_monomial(apre::ActionPolyRingElem) = is_monomial(__poly(apre))

@doc raw"""
    total_degree(p::ActionPolyRingElem) -> Int

Return the total degree of `p`.
"""
total_degree(apre::ActionPolyRingElem) = total_degree(__poly(apre))

@doc raw"""
    degree(p::ActionPolyRingElem, i::Int, jet::Vector{Int}) -> Int

Return the degree of the polynomial `p` in the `i`-th elementary variable with
multiindex `jet`. If this jet variable is valid but still untracked, return $0$.
"""
function degree(apre::ActionPolyRingElem, i::Int, jet::Vector{Int})
  apr = parent(apre)
  @req __is_valid_jet(apr, i, jet) "invalid jet variable"
  jtv = apr.jet_to_var
  if haskey(jtv, (i,jet))
    idx = findfirst(var -> var == __poly(jtv[(i,jet)]), gens(apr.upoly_ring))
    return degree(__poly(apre), gen(apr.upoly_ring, idx))
  end
  return 0
end

@doc raw"""
    vars(p::ActionPolyRingElem)

Return the variables actually occuring in `p`.
"""
vars(apre::ActionPolyRingElem) = parent(apre).(vars(__poly(apre)))

@doc raw"""
    gen(apr::ActionPolyRing, i::Int, jet::Vector{Int}) -> Int

Return the the `i`-th elementary variable with multiindex `jet` in the action polynomial
`apr`. If this jet variable was untracked, it is tracked afterwards.
"""
function gen(apr::ActionPolyRing, i::Int, jet::Vector)
  @req __is_valid_jet(apr, i, jet) "invalid jet variable"
  jtv = apr.jet_to_var
  if haskey(jtv, (i,jet))
    return jtv[(i, jet)]
  end
  __add_new_jetvar!(apr, i, jet)
  return jtv[(i,jet)]
end

@doc raw"""
    gens(apr::ActionPolyRing) -> Vector{ActionPolyRingElem}

Return the currently tracked variables of the action polynomial ring `apr` as a vector.
"""
function gens(apr::ActionPolyRing)
  #Probably sort this later
  return collect(values(apr.jet_to_var))
end

@doc raw"""
    number_of_generators(apr::ActionPolyRing)

Return the number of generators of the action polynomial ring `apr`.
"""
number_of_generators(apr::ActionPolyRing) = number_of_generators(apr.upoly_ring)

@doc raw"""
    diff_action(p::ActionPolyRingElem, i::Int) -> ActionPolyRingElem

Apply the `i`-th diff-action to the polynomial `p`, that is, increase the multiindex of each
variable occuring in `p` at position `i` by one and return the resulting polynomial.
"""
function diff_action(apre::ActionPolyRingElem, i::Int)
  apr = parent(apre)
  @req i in 1:ndiffs(apr) "index out of range"
  jet_idxs = map(var -> apr.var_to_jet[var], vars(apre))
  for j in 1:length(jet_idxs)
    jet_idx = jet_idxs[j]
    jet_idx[2][i] += 1 #increase the i-th multiindex by 1
    if !haskey(apr.jet_to_var, jet_idx) #if this leads to an untracked variable...
      __add_new_jetvar!(apr, jet_idx[1], jet_idx[2]) #track it
    end
    
  end
end

#######################################
#
#  Aux 
#
#######################################

__poly(apre::ActionPolyRingElem) = apre.p

#Check if the jet_to_var dictionary of apr could contain the key (i,jet).
__is_valid_jet(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = i in 1:nelementary_symbols(apr) && length(jet) == ndiffs(apr) && all(j -> j>=0, jet)

function __add_new_jetvar!(apr::ActionPolyRing, i::Int, jet::Vector{Int})
  jtv = apr.jet_to_var
  upoly_new_var = gen(apr.upoly_ring, string(apr.elementary_symbols[i]) * "[" * join(jet) * "]")
  jtv[(i, jet)] = apr(upoly_new_var)
  apr.var_to_jet[apr(upoly_new_var)] = (i, jet)
  return jtv[(i, jet)]
end

