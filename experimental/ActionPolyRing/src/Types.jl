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
  jet_to_upoly_idx::Dict{Tuple{Int, Vector{Int}}, Int}

  function DifferencePolyRing{T}(R::Ring, nelementary_symbols::Int, ndiffs::Int, internal_ordering::Tuple{Symbol, Symbol}) where {T}
    @req nelementary_symbols >= 0 "The number of elementary symbols must be nonnegative"
    @req ndiffs >= 0 "The number of endomorphisms must be nonnegative"
    @req internal_ordering[1] in [:lex, :deglex, :degrevlex] "ordering of the elementary variables must be one of :lex, :deglex, :degrevlex"
    @req internal_ordering[2] in [:top, :pot] "extension must be one of :top (term-over-position) or :pot (position-over-term)"
    elementary_symbols = map(x -> Symbol("u" * string(x)), 1:nelementary_symbols)
    upoly_ring = universal_polynomial_ring(R)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
    
    return new{T}(upoly_ring, elementary_symbols, ndiffs, internal_ordering, jet_to_var, var_to_jet, jet_to_upoly_idx)
  end
 
  function DifferencePolyRing{T}(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int, internal_ordering::Tuple{Symbol, Symbol}) where {T}
    @req ndiffs >= 0 "The number of endomorphisms must be nonnegative"
    @req internal_ordering[1] in [:lex, :deglex, :degrevlex] "ordering of the elementary variables must be one of :lex, :deglex, :degrevlex"
    @req internal_ordering[2] in [:top, :pot] "extension must be one of :top (term-over-position) or :pot (position-over-term)"
    upoly_ring = universal_polynomial_ring(R)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
    
    return new{T}(upoly_ring, elementary_symbols, ndiffs, internal_ordering, jet_to_var, var_to_jet, jet_to_upoly_idx)
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
  for i in 1:nelementary_symbols
    __add_new_jetvar!(dpr, i, fill(0, ndiffs))
  end
  return (dpr, map(i -> dpr.jet_to_var[(i, fill(0, ndiffs))], 1:nelementary_symbols))
end

@doc raw"""
    difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int) -> Tuple{DifferencePolyRing, Vector{DifferencePolyRingElem}}

Return a tuple consisting of the difference polynomial ring over the ring `R` with the specified elementary variables and number of commuting endomorphisms, and the vector of
these elementary variables. Note that the multiindex [0..0] of length 'ndiffs' is appended to the variable names provided.
"""
function difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int; internal_ordering = (:lex, :top)::Tuple{Symbol, Symbol})
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, elementary_symbols, ndiffs, internal_ordering)
  for i in 1:length(elementary_symbols)
    __add_new_jetvar!(dpr, i, fill(0, ndiffs))
  end
  return (dpr, map(i -> dpr.jet_to_var[(i, fill(0, ndiffs))], 1:length(elementary_symbols)))
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
    trailing_coefficient(p::ActionPolyRingElem)

Return the trailing coefficient of the polynomial `p`, i.e. the coefficient of the last nonzero term, or zero if the polynomial is zero.
"""
trailing_coefficient(apre::ActionPolyRingElem) = parent(apre)(trailing_coefficient(__poly(apre)))

@doc raw"""
    is_constant(p::ActionPolyRingElem)

Return `true` if `p` is a degree zero polynomial or the zero polynomial, i.e. a constant polynomial. 
"""
is_constant(apre::ActionPolyRingElem) = is_constant(__poly(apre))

@doc raw"""
    vars(p::ActionPolyRingElem)

Return the variables actually occuring in `p`.
"""
vars(apre::ActionPolyRingElem) = parent(apre).(vars(__poly(apre)))

#combine_like_terms!
#sort_terms!

@doc raw"""
    gen(apr::ActionPolyRing, i::Int, midx::Vector{Int})

Return the `i`-th elementary variable with multiindex `midx` in the action polynomial
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
function diff_action(apre::ActionPolyRingElem{T}, i::Int) where {T}
    apr = parent(apre)
    @req i in 1:ndiffs(apr) "index out of range"

    # Remove constant term: d_i(1) = 0
    apre -= trailing_coefficient(apre)

    if is_zero(apre)
        return apre
    end

    upre = __poly(apre)

    # Precompute mapping from existing variable -> shifted variable position
    upoly_var_to_shifted_pos = Dict{AbstractAlgebra.Generic.UnivPoly{T}, Int}()
    for var in vars(upre)
        jet_idx = apr.var_to_jet[apr(var)]
        new_jet_idx = copy(jet_idx[2])
        new_jet_idx[i] += 1

        if !haskey(apr.jet_to_upoly_idx, (jet_idx[1], new_jet_idx))
            __add_new_jetvar!(apr, jet_idx[1], new_jet_idx)
        end

        shifted_var_elem = apr.jet_to_var[(jet_idx[1], new_jet_idx)]
        shifted_var = __poly(shifted_var_elem)
        shifted_var_pos = findfirst(==(shifted_var), gens(apr.upoly_ring))
        upoly_var_to_shifted_pos[var] = shifted_var_pos
    end

    C = MPolyBuildCtx(apr.upoly_ring)

    for term in terms(upre)
        coeff_t = coeff(term, 1)
        vars_t = vars(term)
        ev = append!(exponent_vector(term, 1), fill(0, length(exponent_vector(term, 1))))

        for var in vars_t
            exp = ev[findfirst(==(var), gens(apr.upoly_ring))]

            jet_idx = apr.var_to_jet[apr(var)]
            var_pos = findfirst(==(var), gens(apr.upoly_ring))
            shifted_var_pos = upoly_var_to_shifted_pos[var]

            new_exp_vec = copy(ev)
            new_exp_vec[var_pos] -= 1
            new_exp_vec[shifted_var_pos] += 1

            push_term!(C, coeff_t * exp, new_exp_vec)
        end
    end

    return apr(finish(C))
end

@doc raw"""
    getindex(apr::ActionPolyRing, i::Int, midx::Vector{Int})

Alias for `gen(apr, i, midx)`.
"""
getindex(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = gen(apr, i, jet)

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
  apr.jet_to_upoly_idx[(i,jet)] = ngens(apr.upoly_ring)
  apr.var_to_jet[apr(upoly_new_var)] = (i, jet)
  return jtv[(i, jet)]
end

