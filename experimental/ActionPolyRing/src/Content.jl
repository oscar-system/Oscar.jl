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
these elementary variables.This methods allows for all the keywords of the `set_ranking!` method, as well as the keywords `cached` and `internal_ordering` from the MPoly Interface. See there for further information.
"""
function difference_polynomial_ring(R::Ring, nelementary_symbols::Int, ndiffs::Int;
    partition_name::Symbol = :default,
    index_ordering_name::Symbol = :default,
    partition::Vector{Vector{Int}} = Vector{Int}[],
    index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0),
    cached::Bool = true,
    internal_ordering::Symbol = :lex
  )
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, nelementary_symbols, ndiffs, cached, internal_ordering)
  set_ranking!(dpr;
      partition_name = partition_name,
      index_ordering_name = index_ordering_name,
      partition = partition,
      index_ordering_matrix = index_ordering_matrix
    )
  return (dpr, __add_new_jetvar!(dpr, map(i -> (i, fill(0, ndiffs)), 1:nelementary_symbols); provide_info = false))
end

@doc raw"""
    difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int) -> Tuple{DifferencePolyRing, Vector{DifferencePolyRingElem}}

Return a tuple consisting of the difference polynomial ring over the ring `R` with the specified elementary variables and number of commuting endomorphisms, and the vector of
these elementary variables. Note that the multiindex [0..0] of length 'ndiffs' is appended to the variable names provided. This methods allows for all the keywords of the `set_ranking!`
method, as well as the keywords `cached` and `internal_ordering` from the MPoly Interface. See there for further information.
"""
function difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int;
    partition_name::Symbol = :default,
    index_ordering_name::Symbol = :default,
    partition::Vector{Vector{Int}} = Vector{Int}[],
    index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0),
    cached::Bool = true,
    internal_ordering::Symbol = :lex
  )
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, elementary_symbols, ndiffs, cached, internal_ordering)
  set_ranking!(dpr;
      partition_name = partition_name,
      index_ordering_name = index_ordering_name,
      partition = partition,
      index_ordering_matrix = index_ordering_matrix
    )
  return (dpr, __add_new_jetvar!(dpr, map(i -> (i, fill(0, ndiffs)), 1:length(elementary_symbols)); provide_info = false))
end

### Differential ###

##### Elements #####

### Difference ###
(dpr::DifferencePolyRing)() = DifferencePolyRingElem{elem_type(base_ring(dpr))}(dpr)

(dpr::DifferencePolyRing)(upre::AbstractAlgebra.Generic.UniversalPolyRingElem) = DifferencePolyRingElem{elem_type(base_ring(dpr))}(dpr, upre)

function (dpr::DifferencePolyRing{T})(a::DifferencePolyRingElem{T}) where {T}
  @req parent(a) === dpr "Wrong parent"
  return a
end

(dpr::DifferencePolyRing)(a::T) where {T<:RingElement} = dpr(dpr.upoly_ring(a))

### Differential ###

#######################################
#
#  Basic ring functionality 
#
#######################################

elem_type(::Type{DifferencePolyRing{T}}) where {T} = DifferencePolyRingElem{T}

parent_type(::Type{DifferencePolyRingElem{T}}) where {T} = DifferencePolyRing{T}

base_ring_type(::Type{<:ActionPolyRing{T}}) where {T} = parent_type(T)

is_domain_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_domain_type(T)

is_exact_type(::Type{<:ActionPolyRingElem{T}}) where {T} = is_exact_type(T)

function Base.deepcopy_internal(dpre::DifferencePolyRingElem{T}, dict::IdDict) where {T}
    # Avoid deepcopying the parent as it may refer back to it in one of its dictionaries 
    pp = deepcopy_internal(data(dpre), dict)
    return DifferencePolyRingElem{T}(parent(dpre), pp)
end

zero(dpr::DifferencePolyRing) = dpr(0)

one(dpr::DifferencePolyRing) = dpr(1)

is_square(apre::ActionPolyRingElem) = is_square(data(apre))

Base.sqrt(apre::ActionPolyRingElem; check::Bool = true) = parent(apre)(sqrt(data(apre); check = check))

is_irreducible(apre::ActionPolyRingElem) = is_irreducible(data(apre))

is_unit(apre::ActionPolyRingElem) = is_unit(data(apre))

characteristic(apr::ActionPolyRing) = characteristic(__upr(apr))

### Factorization ###
function __wrap_factorization_apr(f::Fac{<:UniversalPolyRingElem}, apr::ActionPolyRing)
   res = Fac{elem_type(apr)}()
   res.unit = apr(f.unit)
   for (fact, expo) in f
      AbstractAlgebra.mulpow!(res, apr(fact), expo)
   end
   return res
end

factor_squarefree(apr::ActionPolyRingElem) = __wrap_factorization_apr(factor_squarefree(data(apr)), parent(apr))

factor(apr::ActionPolyRingElem) = __wrap_factorization_apr(factor(data(apr)), parent(apr))

#######################################
#
#  Field access 
#
#######################################

##### Algebras #####

### Difference ###
base_ring(dpr::DifferencePolyRing) = base_ring(__upr(dpr))

ndiffs(dpr::DifferencePolyRing) = dpr.ndiffs

elementary_symbols(dpr::DifferencePolyRing) = dpr.elementary_symbols

nelementary_symbols(dpr::DifferencePolyRing) = length(elementary_symbols(dpr))

##### Elements #####

parent(dpre::DifferencePolyRingElem) = dpre.parent

#######################################
#
#  Basic polynomial functionality 
#
#######################################

is_monomial(apre::ActionPolyRingElem) = is_monomial(data(apre))

is_univariate(apre::ActionPolyRingElem) = is_univariate(data(apre))

is_univariate(apr::ActionPolyRing) = false

length(apre::ActionPolyRingElem) = length(data(apre))

@doc raw"""
    total_degree(p::ActionPolyRingElem) -> Int

Return the total degree of `p`.
"""
total_degree(apre::ActionPolyRingElem) = total_degree(data(apre))

@doc raw"""
    degree(p::ActionPolyRingElem, i::Int, jet::Vector{Int}) -> Int

Return the degree of the polynomial `p` in the `i`-th elementary variable with
multiindex `jet`. If this jet variable is valid but still untracked, return $0$.
"""
function degree(apre::ActionPolyRingElem, i::Int, jet::Vector{Int})
  apr = parent(apre)
  upr = __upr(apr)
  @req __is_valid_jet(apr, i, jet) "Invalid jet variable"
  jtv = __jtv(apr)
  if haskey(jtv, (i,jet))
    idx = findfirst(var -> var == data(jtv[(i,jet)]), gens(upr))
    return degree(data(apre), gen(upr, idx))
  end
  return 0
end

degree(apre::ActionPolyRingElem, jet_idx::Tuple{Int, Vector{Int}}) = degree(apre, jet_idx...)

@doc raw"""
    degree(p::ActionPolyRingElem, i::Int)

Return the degree of the polynomial `p` in the, among the currently tracked variables, 'i'-th largest one.
"""
degree(apre::ActionPolyRingElem, i::Int) = degree(apre, __vtj(parent(apre))[gen(parent(apre), i)])

@doc raw"""
    degrees(p::ActionPolyRingElem)

Return an array of the degrees of the polynomial `p` in terms of each variable. The variables are sorted with respect to the ranking of the action polynomial ring containing `p`, leading with the largest variable.
"""
degrees(apre::ActionPolyRingElem) = map(i -> degree(apre, i), 1:ngens(parent(apre)))

#=
@doc raw"""
    trailing_coefficient(p::ActionPolyRingElem)

Return the trailing coefficient of the polynomial `p`, i.e. the coefficient of the last nonzero term, or zero if the polynomial is zero.
"""
trailing_coefficient(apre::ActionPolyRingElem) = parent(apre)(trailing_coefficient(data(apre)))
=#

@doc raw"""
    is_constant(p::ActionPolyRingElem)

Return `true` if `p` is a degree zero polynomial or the zero polynomial, i.e. a constant polynomial. 
"""
is_constant(apre::ActionPolyRingElem) = is_constant(data(apre))

@doc raw"""
    vars(p::ActionPolyRingElem)

Return the variables actually occuring in `p`. The variables are sorted with respect to the ranking of the action polynomial ring containing `p`, leading with the largest variable.
"""
vars(apre::ActionPolyRingElem) = sort(parent(apre).(vars(data(apre))); rev = true)

is_gen(apre::ActionPolyRingElem) = is_gen(data(apre))

@doc raw"""
    gen(apr::ActionPolyRing, i::Int, midx::Vector{Int})

Return the `i`-th elementary variable with multiindex `midx` in the action polynomial
`apr`. If this jet variable was untracked, it is tracked afterwards.
"""
function gen(apr::ActionPolyRing, i::Int, jet::Vector)
  @req __is_valid_jet(apr, i, jet) "invalid jet variable"
  jtv = __jtv(apr)
  if haskey(jtv, (i, jet))
    return jtv[(i, jet)]
  end
  __add_new_jetvar!(apr, i, jet)
  jtv[(i, jet)]
end

gen(apr::ActionPolyRing, jet_idx::Tuple{Int, Vector{Int}}) = gen(apr, jet_idx...)

@doc raw"""
    gen(apr::ActionPolyRing, i::Int)

Return the, among the currently tracked variables of `apr`, 'i'-th largest one.
"""
gen(apr::ActionPolyRing, i::Int) = gens(apr)[i]

@doc raw"""
    gens(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})

Return the jet variables of the action polynomial ring `apr` specified by the entries of `jet_idxs` as
a vector and track all new variables.
"""
function gens(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})
  jtv = __jtv(apr)
  new_jet_idxs = filter(jet_idx -> !haskey(jtv, jet_idx), jet_idxs)
  __add_new_jetvar!(apr, new_jet_idxs)
  return map(jet_idx -> jtv[jet_idx], jet_idxs)
end

@doc raw"""
    gens(apr::ActionPolyRing)

Return the currently tracked variables of the action polynomial ring `apr` as a vector. The variables are sorted with respect to the ranking of `apr`, leading with the largest variable.
"""
function gens(apr::ActionPolyRing)
  return sort(collect(values(__jtv(apr))); rev = true)
end

number_of_generators(apr::ActionPolyRing) = number_of_generators(__upr(apr))

number_of_variables(apr::ActionPolyRing) = number_of_variables(__upr(apr))

@doc raw"""
    getindex(apr::ActionPolyRing, i::Int, midx::Vector{Int})

Alias for `gen(apr, i, midx)`.
"""
getindex(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = gen(apr, i, jet)

internal_ordering(apr::ActionPolyRing) = internal_ordering(__upr(apr))

#######################################
#
#  Action polynomial functionality 
#
#######################################

@doc raw"""
    diff_action(p::ActionPolyRingElem, i::Int) -> ActionPolyRingElem

Apply the `i`-th diff-action to the polynomial `p`.
"""
function diff_action(apre::ActionPolyRingElem{T}, i::Int) where {T}
  apr = parent(apre)
  @req i in 1:ndiffs(apr) "index out of range"

  # Remove constant term: d_i(1) = 0
  apre -= trailing_coefficient(apre)

  if is_zero(apre)
    return apre
  end

  upre = data(apre)
  upr = __upr(apr)

  # Precompute mapping from existing variable -> shifted variable position
  upoly_var_to_shifted_pos = Dict{AbstractAlgebra.Generic.UnivPoly{T}, Int}()
  for var in vars(upre)
    jet_idx = __vtj(apr)[apr(var)]
    new_jet_idx = copy(jet_idx[2])
    new_jet_idx[i] += 1

    if !haskey(__jtu_idx(apr), (jet_idx[1], new_jet_idx))
      __add_new_jetvar!(apr, jet_idx[1], new_jet_idx; provide_info = false)
    end

    shifted_var_elem = __jtv(apr)[(jet_idx[1], new_jet_idx)]
    shifted_var = data(shifted_var_elem)
    shifted_var_pos = findfirst(==(shifted_var), gens(upr))
    upoly_var_to_shifted_pos[var] = shifted_var_pos
  end

  C = MPolyBuildCtx(upr)

  for term in terms(upre)
    coeff_t = coeff(term, 1)
    vars_t = vars(term)
    ev = append!(exponent_vector(term, 1), fill(0, length(exponent_vector(term, 1))))

    for var in vars_t
      exp = ev[findfirst(==(var), gens(upr))]

      jet_idx = __vtj(apr)[apr(var)]
      var_pos = findfirst(==(var), gens(upr))
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
    diff_action(p::ActionPolyRingElem, d::Vector{Int}) -> ActionPolyRingElem

Successively apply the `i`-th diff-action `d[i]`-times to the polynomial `p`, where $i = 1, \ldots, length(d)$. 
"""
function diff_action(apre::ActionPolyRingElem{T}, d::Vector{Int}) where {T}
  len = length(d)
  @req len == ndiffs(parent(apre)) && all(j -> j >= 0, d) "Invalid vector of diff multiplicities"
  res = apre
  for i in 1:len
    for j in 1:d[i]
      res = diff_action(res, i)
    end
  end
  return res
end

@doc raw"""
    leader(p::ActionPolyRingElem)

Return the leader of the polynomial `p`, that is the largest variable with respect to the ranking of `parent(p)`.
"""
leader(apre::ActionPolyRingElem) = max(vars(apre)...)

#######################################
#
#  Misc 
#
#######################################


#######################################
#
#  Aux Action polynomial rings 
#
#######################################

#Getters for internals
data(dpre::DifferencePolyRingElem) = dpre.p

__upr(dpr::DifferencePolyRing) = dpr.upoly_ring

__jtv(dpr::DifferencePolyRing) = dpr.jet_to_var

__vtj(dpr::DifferencePolyRing) = dpr.var_to_jet

__jtu_idx(dpr::DifferencePolyRing) = dpr.jet_to_upoly_idx

#Setters for internals
function __set_are_perms_up_to_date!(dpr::DifferencePolyRing, update::Bool)
  dpr.are_perms_up_to_date = update
end

function __set_perm_for_sort!(dpr::DifferencePolyRing)
  __set_are_perms_up_to_date!(dpr, true)
  dpr.permutation = sortperm(dpr.(gens(__upr(dpr))); rev = true)
end

__perm_for_sort(dpr::DifferencePolyRing) = dpr.permutation

#Check if the jet_to_var dictionary of apr could contain the key (i,jet).
__is_valid_jet(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = i in 1:nelementary_symbols(apr) && length(jet) == ndiffs(apr) && all(j -> j >= 0, jet)

function __add_new_jetvar!(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}}; provide_info::Bool = true)
  __set_are_perms_up_to_date!(apr, false)
  s_vec = map(jet_idx -> string(elementary_symbols(apr)[jet_idx[1]]) * "[" * join(jet_idx[2]) * "]", jet_idxs)
  if provide_info
    @info "New variables: " * join(s_vec, ", ")
  end
  upr = __upr(apr)
  new_vars = apr.(gens(upr, s_vec))
  for k in 1:length(new_vars)
    var = new_vars[k]
    jetvaridx = jet_idxs[k]
    __jtv(apr)[jetvaridx] = var
    __vtj(apr)[var] = jetvaridx
    __jtu_idx(apr)[jetvaridx] = ngens(upr)
  end
  return new_vars
end

__add_new_jetvar!(apr::ActionPolyRing, i::Int, jet::Vector{Int}; provide_info::Bool = true) = __add_new_jetvar!(apr, [(i, jet)]; provide_info = provide_info)

#######################################
#
#  Construction / Getter
#
#######################################

@doc raw"""
    ranking(dpr::DifferencePolyRing) -> DifferenceRanking

Return the ranking of the jet variables of the difference polynomial ring `dpr`
"""
ranking(dpr::DifferencePolyRing) = dpr.ranking
  
@doc raw"""
    set_ranking!(dpr::DifferencePolyRing;
                 partition_name::Symbol = :default,
                 index_ordering_name::Symbol = :default,
                 partition::Vector{Vector{Int}} = Vector{Int}[],
                 index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0)) 
                 -> DifferenceRanking

This method configures the ranking of the difference polynomial ring `dpr`, using an ordered partition of the elementary symbols and a monomial ordering on the indices. The ranking can be specified either by choosing predefined naming options or by explicitly providing a custom configuration.

# Keyword Arguments
- `partition_name`: Determines the partition of the elementary symbols of `dpr`. Supported values are:
  - `:top`: groups all variables into a single block,
  - `:pot`: separates each variable into its own block,
  - `:default`: uses `:top` unless a custom partition is specified.

- `index_ordering_name`: Specifies the ordering on the multiindices. Supported values are:
  - `:lex`: lexicographic ordering,
  - `:deglex`: degree lexicographic ordering,
  - `:invlex`: inverse lexicographic ordering,
  - `:deginvlex`: degree inverse lexicographic ordering,
  - `:degrevlex`: degree reverse lexicographic ordering,
  - `:default`: uses `:lex` unless a custom matrix is specified.

- `partition`: A custom partition of the elementary symbols, represented as a vector of characteristic vectors. The elementary symbols corresponding to the first characteristic vectors are considered largest and so on.

- `index_ordering_matrix`: A custom matrix representing a monomial ordering on the indices. Its number of columns must equal `ndiffs(dpr)`.
"""
function set_ranking!(dpr::DifferencePolyRing;
    partition_name::Symbol = :default,
    index_ordering_name::Symbol = :default,
    partition::Vector{Vector{Int}} = Vector{Int}[],
    index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0)
  )::DifferenceRanking{elem_type(base_ring(dpr))}
  
  __set_are_perms_up_to_date!(dpr, false)
  m = nelementary_symbols(dpr)
  @req partition_name in [:top, :pot, :default] "Invalid name of partition"
  if partition_name == :default
    if is_empty(partition)
      partition = [fill(1, m)] #Use :top by default
    else
      # Otherwise the input is used. Check its validity:
      @req __is_valid_partition(partition, m) "Not a partition of the number of elementary symbols"
    end
  else
    if is_empty(partition)
      if partition_name == :top
        partition = [fill(1, m)]
      else
        partition = [[i == j ? 1 : 0 for j in 1:m] for i in 1:m]
      end
    else # This case is only accessed if both a partition and a name are provided. Then a consistency check is required.
      @req __is_valid_partition(partition, m) "Not a partition of the number of elementary symbols"
      if partition_name == :top
        @req partition == [fill(1, m)] "The partition provided does not match its name"
      else
        @req partition == [[i == j ? 1 : 0 for j in 1:m] for i in 1:m] "The partition provided does not match its name"
      end
    end
  end

  n = ndiffs(dpr)
  @req index_ordering_name in [:default, :lex, :deglex, :invlex, :deginvlex, :degrevlex] "Invalid name of index ordering"
  R, _ = polynomial_ring(QQ, n; cached = false) # We use a dummy polynomial ring to extract canonical matrices of monomial orderings
  if index_ordering_name == :default
    if is_empty(index_ordering_matrix)
      index_ordering_matrix = canonical_matrix(lex(R)) #Use :lex by default
    else
      # Otherwise the input is used. Check its validity:
      @req ncols(index_ordering_matrix) == n "The number of columns of the matrix provided must equal $n" 
    end
  else
    @req is_empty(index_ordering_matrix) "Providing both a name and a matrix is not supported. Please just choose one."
    #if is_empty(index_ordering_matrix)
      if index_ordering_name == :lex
        index_ordering_matrix = canonical_matrix(lex(R))
      elseif index_ordering_name == :deglex
        index_ordering_matrix = canonical_matrix(deglex(R))
      elseif index_ordering_name == :invlex
        index_ordering_matrix = canonical_matrix(invlex(R))
      elseif index_ordering_name == :deginvlex
        index_ordering_matrix = canonical_matrix(deginvlex(R))
      elseif index_ordering_name == :degrevlex
        index_ordering_matrix = canonical_matrix(degrevlex(R))
      end
    #end
  end
  ran = DifferenceRanking{elem_type(base_ring(dpr))}(dpr, partition, index_ordering_matrix)
  dpr.ranking = ran
  return ran
end

#######################################
#
#  basic functionality 
#
#######################################

base_ring(ran::DifferenceRanking) = ran.ring

@doc raw"""
    partition(r::DifferenceRanking) -> Vector{Vector{Int}}

Return the partition of the elementary symbols defined by the ranking `r`
of the difference polynomial ring `dpr`, where `r = ranking(dpr)`.
"""
partition(ran::DifferenceRanking) = ran.partition

@doc raw"""
    index_ordering_matrix(r::DifferenceRanking) -> ZZMatrix

Return the matrix inducing the monomial ordering of the multiindices defined
by the ranking `r` of the difference polynomial ring `dpr`, where `r = ranking(dpr)`.
"""
index_ordering_matrix(ran::DifferenceRanking) = ran.index_ordering_matrix

function riquier_matrix(ran::DifferenceRanking)
  if !isdefined(ran, :riquier_matrix)
    par = partition(ran)
    dpr = base_ring(ran)
    upper_part = block_diagonal_matrix([matrix(ZZ, length(par)-1, nelementary_symbols(dpr), vcat(par[1:end-1]...)), index_ordering_matrix(ran)])
    lower_part = block_diagonal_matrix([__in_block_tie_breaking_matrix(par), zero_matrix(ZZ, 0, ndiffs(dpr))])
    ran.riquier_matrix = vcat(upper_part, lower_part)
  end
  return ran.riquier_matrix
end

#######################################
#
#  Aux Rankings 
#
#######################################

__is_valid_partition(par::Vector{Vector{Int}}, m::Int) = all(par_elt -> !is_zero(par_elt), par) && sum(par) == fill(1, m) && all(par_elt -> all(j -> par_elt[j] == 0 || par_elt[j] == 1, 1:m), par)

function __in_block_tie_breaking_matrix(partition::Vector{Vector{Int}})
    m = length(partition[1])
    max_ones = maximum(count(isone, vec) for vec in partition)
    result = zero_matrix(ZZ, max_ones - 1, m)

    for vec in partition
        pos = 1
        ones_seen = 0
        total_ones = count(isone, vec)

        for j in 1:m
            if isone(vec[j])
                ones_seen += 1
                if ones_seen == total_ones
                    break
                end
                result[pos, j] = ZZ(1)
                pos += 1
            end
        end
    end

    return result
end

