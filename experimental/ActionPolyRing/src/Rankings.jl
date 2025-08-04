# For a ranking of jet variables, e.g. u1[123], u3[041], ... one needs:
# - An ordering to compare the elementary symbols u_1, ..., u_n
# - An ordering of the multiindices (this is like a monomial_ordering)
# - A decision on which of the two has priority, e.g. position-over-term (pot)
#   or term-over-position (top). Here, the position of the jet variable corresponding
#   to (i, [a_1, ..., a_n]) is the first coordinate, i.e. i. More generally, this decision means
#   choosing an ordered partition of the set {1, ...,m}, pot and top being the choices of the two trivial
#   partitions.

export Ranking,
       DifferenceRanking,
       ranking,
       set_ranking!,
       partition,
       index_ordering_matrix,
       riquier_matrix

#######################################
#
#  Data types 
#
#######################################

abstract type Ranking end

mutable struct DifferenceRanking{T} <: Ranking where {T}
  ring::DifferencePolyRing{T}
  partition::Vector{Vector{Int}}
  index_ordering_matrix::ZZMatrix
  riquier_matrix::ZZMatrix

  function DifferenceRanking{T}(dpr::DifferencePolyRing{T}, partition::Vector{Vector{Int}}, index_ordering_matrix::ZZMatrix) where {T}  
    return new{T}(dpr, partition, index_ordering_matrix)
  end
  
end

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
  R, _ = polynomial_ring(QQ, n) # We use a dummy polynomial ring to extract canonical matrices of monomial orderings
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
#  I/O 
#
#######################################

function Base.show(io::IO, ::MIME"text/plain", ran::DifferenceRanking)
  io = pretty(io)
  print(io, "Ranking of $(base_ring(ran))\n")
  print(io, "with elementary symbols partitioned by\n")
  print(io, Indent())
  print(io, partition(ran))
  print(io, Dedent())
  print(io, "\nand ordering of the indices defined by\n")
  print(io, Indent())
  show(io, "text/plain", index_ordering_matrix(ran))
  print(io, Dedent())
end

function Base.show(io::IO, ran::DifferenceRanking)
  io = pretty(io)
  if is_terse(io)
    print(io, "Ranking of difference polynomial ring")
  end
  print(terse(io), "Ranking of difference polynomial ring")
end

#######################################
#
#  Aux 
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

