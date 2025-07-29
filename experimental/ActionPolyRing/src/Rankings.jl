# For a ranking of jet variables, e.g. u1[123], u3[041], ... one needs:
# - An ordering to compare the elementary symbols u_1, ..., u_n
# - An ordering of the multiindices (this is like a monomial_ordering)
# - A decision on which of the two has priority, e.g. position-over-term (pot)
#   or term-over-position (top). Here, the position of the jet variable corresponding
#   to (i, [a_1, ..., a_n]) is the first coordinate, i.e. i. More generally, this decision means
#   choosing a partition of the set {1, ...,m}, pot and top being the choices of the two trivial
#   partitions.
#   
#   We represent a ranking by a Riquier matrix M. Then two variables are compared as follows:
#   We embed (i, [a_1, ..., a_n]) -> [a_1, ..., a_n, e_i]^T where e_i is the i-th unit row.
#   Then multiply the columns from the right to M and compare the entries lexicographically.


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

function ranking(dpr::DifferencePolyRing)
  if !isdefined(dpr, :ranking)
    set_ranking!(dpr) #Default ranking
    return dpr.ranking
  end
  return dpr.ranking
end

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
      partition == [fill(1, m)] #Use :top by default
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
  R, _ = polynomial_ring(QQ, n) # We use a dummy polynomial ring to extract coefficient matrices of monomial orderings
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

partition(ran::DifferenceRanking) = ran.partition

index_ordering_matrix(ran::DifferenceRanking) = ran.index_ordering_matrix

function riquier_matrix(ran::DifferenceRanking)
  if !is_defined(ran, :riquier_matrix)
    #compute it
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

__is_valid_partition(par::Vector{Vector{Int}}, m::Int) = sum(par) == fill(1, m) && all(par_elt -> all(j -> par_elt[j] == 0 || par_elt[j] == 1, 1:m), par)

