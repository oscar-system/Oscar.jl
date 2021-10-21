module Multiindices

export OrderedMultiindex
export linear_index
export index_signature
export HomogMultiindex
export power
export getindex


abstract type Multiindex end

mutable struct OrderedMultiindex <: Multiindex 
  n::Int
  p::Int
  index::Vector{Int}
  function OrderedMultiindex(n::Int, p::Int)
    return new( n, p, collect(1:p) )
  end
end

function Base.show( io::Base.TTY, i::OrderedMultiindex ) 
  outstr= "0"
  for k in i.index
    outstr *= " < " * string(i.index[k])
  end
  outstr = outstr * " <= " * string( i.n ) * "\n"
  #Base.show( io, "n = " * string(i.n) * "; p = " * string(i.p) * "; index = " * string( i.index ) )
  Base.print( io, outstr )
end

function Base.iterate( i::OrderedMultiindex ) 
  for k in 1:i.p
    i.index[k] = k
  end
  return i, i.index
end

function Base.iterate( i::OrderedMultiindex, state ) 
  # determine the first entry that has not yet reached its 
  # limit 
  k = i.p
  while ( k>0 ) && ( i.index[k] == i.n - (i.p-k) )
    k = k-1
  end

  # whenever this statement is true, the iteration is exhausted
  if k == 0 
    return nothing
  end

  # otherwise iterate by one
  i.index[k] = i.index[k]+1
  # and adjust all the follow-up indices
  for j in k+1:i.p
    i.index[j] = i.index[j-1]+1
  end
  return i, i.index
end

# return the linear index in the enumeration pattern 
# implemented in the iterator 
function linear_index( i::OrderedMultiindex )
  N = binomial( i.n, i.p )
  for k in i.p:-1:1 
    N = N - binomial( i.n-i.index[k], i.p-k+1 )
  end
  return N
end

function ordered_multiindex( l::Int, n::Int, p::Int )
  i = OrderedMultiindex( n, p )
  if l > binomial( n, p ) 
    return nothing
  end
  for k in 1:p
    while l > binomial( n-i.index[k], p-k )
      l = l - binomial( n-i.index[k], p-k )
      i.index[k] = i.index[k]+1
    end
    for r = k+1:p
      i.index[r] = i.index[r-1]+1
    end
  end
  return i
end

function index_signature( alpha::OrderedMultiindex )
  sign = 1
  for k in alpha.p:-1:1
    if mod( k-alpha.index[k], 2 ) != 0
      sign = sign*(-1)
    end
  end
  return sign
end

mutable struct HomogMultiindex
  n::Int		# The number of variables
  d::Int		# The degree
  a::Vector{Int}	# The multi-exponent
  i::OrderedMultiindex	# The underlying multiindex


  function HomogMultiindex( n::Int, d::Int, a::Vector{Int}, i::OrderedMultiindex )
    n > 0 || error("number of variables must be positive")
    if d < 0
      error( "invalid degree." )
    end
    if sum(a) != d 
      error( "values of a don't match the requested degree." )
    end 
    return new( n, d, a, i )
  end
  
  function HomogMultiindex( n::Int, d::Int )
    if n<= 0
      error( "invalid number of variables." )
    end
    if d<0 
      error( "invalid degree." )
    end
    a = vcat( [ 0 for i in (2:n) ], [d] )
    i = OrderedMultiindex( d+n-1, n-1 ) 
    #@show i
    #@show i.index
    #@show i.n
    return HomogMultiindex( n, d, a, i )
  end
end

function extract_exponent( alpha::HomogMultiindex )
  #@show alpha.i.index, alpha.i.p
  if alpha.i.p > 0
    a = [ alpha.i.index[1] - 1 ]
    append!( a, [ alpha.i.index[k+1] - alpha.i.index[k] - 1 for k in 1:alpha.i.p-1 ] )
    a = vcat( a, [ alpha.i.n - alpha.i.index[alpha.i.p] ])
  else
    a = [ alpha.d ]
  end
  return a
end

function Base.show( io::Base.TTY, alpha::HomogMultiindex ) 
  outstr = "( "
  for k in (1:alpha.n-1) 
    outstr = outstr * string(alpha.a[k]) * ", "
  end
  #@show alpha.a
  #@show alpha.d
  outstr = outstr * "$(alpha.a[alpha.n]) )\n"
  Base.print( io, outstr )
end
  
function Base.iterate( alpha::HomogMultiindex )
  Base.iterate( alpha.i )
  alpha.a = extract_exponent( alpha )
  return alpha, alpha.a 
end

function Base.iterate( alpha::HomogMultiindex, state )
  #@show alpha.i
  x = Base.iterate( alpha.i, state )
  if x == nothing 
    return nothing	# Iteration is exhausted
  end
  #@show x
  #@show typeof( x )
  alpha.i = x[1]
  alpha.a = extract_exponent( alpha )
  return alpha, alpha.a
end

function power( x::Vector, alpha::HomogMultiindex )
  n = length( x )
  if n != alpha.n
    error( "number of variables not compatible." )
  end
  if n == 0 
    return Int(1)
  end
  return prod(k -> x[k]^(alpha.a[k]), 1:n)
end

function linear_index( alpha::HomogMultiindex ) 
  return linear_index( alpha.i )
end

# Use Ordered multiindices for indexing 2-dimensional arrays

function Base.getindex( v::Matrix, α::OrderedMultiindex )
  if α.p < 2 
    error( "Multiindex does not have sufficient length" )
  end
  return Base.getindex( v, α.index[1], α.index[2] )
end

function Base.setindex!( v::Array{Any,2}, cont::Any, α::OrderedMultiindex )
  if α.p < 2 
    error( "Multiindex does not have sufficient length" )
  end
  return Base.setindex!( v, cont, α.index[1], α.index[2] )
end

#Base.setindex!( v::Array{Any,2}, α::OrderedMultiindex ) = Base.setindex!( v, α.index[1], α.index[2] )

function Base.getindex( i::OrderedMultiindex, k::Int )
  return i.index[k]
end

end


