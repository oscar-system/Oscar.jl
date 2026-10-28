# This file implements types of complex reflection groups following the
# classification of Shephard & Todd (1954).
#
# The type of an irreducible reflection group G <= GL_n(C) is its conjugacy
# class in GL_n(C). In general, the type is determined by the dimension of the
# fixed space and the unordered product of the types of the nontrivial
# irreducible factors. This is finer than the isomorphism type of the underlying
# abstract group.
#
# References:
#
# * Lehrer, G. I., & Taylor, D. E. (2009). Unitary reflection groups (Vol. 20, p. viii). Cambridge University Press, Cambridge.
#
# * Thiel, U. (2014). On restricted rational Cherednik algebras. TU Kaiserslautern.
#
# * Geck, M., & Malle, G. (2006). Reflection groups. In Handbook of algebra. Vol. 4 (Vol. 4, pp. 337–383)
#
# Ulrich Thiel, 2023

const ComplexReflectionGroupComponent = Union{Int, Tuple{Int, Int, Int}}

"""
    ComplexReflectionGroupType(; fixed_space_dimension=0)
    ComplexReflectionGroupType(types::Vector; fixed_space_dimension=0)
    ComplexReflectionGroupType(i::Int; fixed_space_dimension=0)
    ComplexReflectionGroupType(m::Int, p::Int, n::Int;
                               fixed_space_dimension=0)

The type of a finite complex reflection group as a complex linear group.

For an irreducible group `G` in `GL(n, C)`, its type is its conjugacy class in
`GL(n, C)`. The Shephard--Todd symbols stored here are classification labels for
the nontrivial irreducible conjugacy classes. In general, the type is the pair
consisting of the dimension of the fixed space and the unordered product of the
nontrivial irreducible types. It therefore distinguishes ambient realizations
that differ by a trivial summand.

The keyword `fixed_space_dimension` supplies a trivial summand in addition to
any trivial one-dimensional factors occurring in `types`. This object stores no
concrete matrix representative, basis, coefficient field, or reference
marking.

The label `(1, 1, n)` for `n >= 2` denotes the essential reflection
representation of the symmetric group, of degree `n - 1`. Its permutation
representation has one additional fixed dimension. A label `(m, m, 1)` is the
trivial one-dimensional representation and is normalized by adding one to the
fixed-space dimension.
"""
struct ComplexReflectionGroupType
  type::Tuple{Vararg{ComplexReflectionGroupComponent}}
  fixed_space_dimension::Int

  # Argument checking Cannot be more specific than Vector for the constructor argument
  # since [ 4, (4,1,4) ] is of type Vector{Any}
  function ComplexReflectionGroupType(
      type::Vector; fixed_space_dimension::IntegerUnion=0
    )
    @req fixed_space_dimension >= 0 "The fixed-space dimension must be nonnegative"
    fixed_dimension = Int(fixed_space_dimension)

    # The normalized type with unique labels (see below)
    typenorm = ComplexReflectionGroupComponent[]

    for t in type
      @req isa(t, Int) || isa(t, Tuple{Int,Int,Int}) "Type must be a vector of integers (for exceptional groups) or of 3-element tuples (for the infinite series)"

      # Will be the normalized type
      tnorm = t

      if isa(t, Int)
        @req t >= 4 && t <= 37 "Exceptional types are between 4 and 37"
      end

      if isa(t, Tuple)
        (m,p,n) = t
        @req m >= 1 && p >= 1 && n >= 1 "m,p,n >= 1 required"
        @req is_divisible_by(m,p) "p must be a divisor of m"
        # Normalization of special overlap cases to get unique labeling
        # See # Lehrer & Taylor (2009). p 27
        if (m,p,n) == (2,2,2)
          # G(2, 2, 2) is the product of two groups of type A1.
          push!(typenorm, (1,1,2))
          push!(typenorm, (1,1,2))
          continue
        elseif (m,p,n) == (4,4,2)
          tnorm = (2,1,2)
        elseif (m,p,n) == (3,3,2)
          tnorm = (1,1,3)
        elseif (m,p,n) == (2,2,3)
          tnorm = (1,1,4)
        elseif n == 1
          m = div(m,p)
          tnorm = (m,1,1)
          if m == 1
            # The trivial group in dimension one is an irreducible trivial
            # representation, hence contributes one fixed dimension.
            fixed_dimension += 1
            continue
          elseif m == 2 # (2,1,1) = (1,1,2)
            tnorm = (1,1,2)
          end
        end

      end
      push!(typenorm, tnorm)
    end

    # The order of the irreducible factors carries no mathematical information.
    sort!(typenorm; by=t -> isa(t, Int) ? (1,t,0,0) : (0,t[1],t[2],t[3]))
    return new(Tuple(typenorm), fixed_dimension)
  end
end

Base.:(==)(G1::ComplexReflectionGroupType, G2::ComplexReflectionGroupType) =
  G1.type == G2.type &&
  fixed_space_dimension(G1) == fixed_space_dimension(G2)
Base.isequal(G1::ComplexReflectionGroupType, G2::ComplexReflectionGroupType) =
  isequal(G1.type, G2.type) &&
  isequal(fixed_space_dimension(G1), fixed_space_dimension(G2))
Base.hash(G::ComplexReflectionGroupType, h::UInt) = hash(
  (G.type, fixed_space_dimension(G)), hash(:ComplexReflectionGroupType, h)
)

# Printing
function Base.show(io::IO, ::MIME"text/plain", G::ComplexReflectionGroupType)
  print(io, "Complex reflection group type ")
  counter = 0
  N = length(G.type)
  if N == 0
    print(
      io,
      "trivial in dimension ",
      fixed_space_dimension(G),
    )
    return
  end
  for t in G.type
    print(io,"G")
    if isa(t, Int)
      print(io,t)
    else
      print(io, "(", t[1], ",", t[2], ",", t[3], ")")
    end
    counter += 1
    if counter < N
      print(io, " x ")
    end
  end
  if fixed_space_dimension(G) > 0
    print(
      io,
      " with ",
      fixed_space_dimension(G),
      "-dimensional fixed space",
    )
  end
end

# Convenience constructors
ComplexReflectionGroupType(; fixed_space_dimension::IntegerUnion=0) =
  ComplexReflectionGroupType(
    ComplexReflectionGroupComponent[]; fixed_space_dimension
  )
ComplexReflectionGroupType(i::Int; fixed_space_dimension::IntegerUnion=0) =
  ComplexReflectionGroupType(
    convert(Vector{Union{Int, Tuple{Int,Int,Int}}}, [i]);
    fixed_space_dimension,
  )
ComplexReflectionGroupType(
    m::Int, p::Int, n::Int; fixed_space_dimension::IntegerUnion=0
  ) = ComplexReflectionGroupType(
    convert(Vector{Union{Int, Tuple{Int,Int,Int}}}, [(m,p,n)]);
    fixed_space_dimension,
  )
ComplexReflectionGroupType(
    t::Tuple{Int,Int,Int}; fixed_space_dimension::IntegerUnion=0
  ) = ComplexReflectionGroupType(
    convert(Vector{Union{Int, Tuple{Int,Int,Int}}}, [t]);
    fixed_space_dimension,
  )

# Non-camel versions
complex_reflection_group_type(; fixed_space_dimension::IntegerUnion=0) =
  ComplexReflectionGroupType(; fixed_space_dimension)
complex_reflection_group_type(i::Int; fixed_space_dimension::IntegerUnion=0) =
  ComplexReflectionGroupType(i; fixed_space_dimension)
complex_reflection_group_type(
    m::Int, p::Int, n::Int; fixed_space_dimension::IntegerUnion=0
  ) = ComplexReflectionGroupType(m, p, n; fixed_space_dimension)
complex_reflection_group_type(
    t::Tuple{Int,Int,Int}; fixed_space_dimension::IntegerUnion=0
  ) = ComplexReflectionGroupType(t; fixed_space_dimension)
complex_reflection_group_type(
    X::Vector; fixed_space_dimension::IntegerUnion=0
  ) = ComplexReflectionGroupType(X; fixed_space_dimension)

# Triviality
function is_trivial(G::ComplexReflectionGroupType)
  return length(G.type) == 0
end

"""
    fixed_space_dimension(G::ComplexReflectionGroupType)

Return the dimension of the subspace fixed pointwise by the represented
reflection group.
"""
fixed_space_dimension(G::ComplexReflectionGroupType) = G.fixed_space_dimension

"""
    components(G::ComplexReflectionGroupType)

Return the nontrivial irreducible reflection-group factors of `G` in canonical
order. The trivial summand is recorded separately by
[`fixed_space_dimension`](@ref).
"""
function components(G::ComplexReflectionGroupType)
  return [ComplexReflectionGroupType(t) for t in G.type]
end

"""
    number_of_components(G::ComplexReflectionGroupType)

Return the number of nontrivial irreducible reflection-group factors of `G`.
"""
function number_of_components(G::ComplexReflectionGroupType)
  return length(G.type)
end

"""
    is_irreducible(G::ComplexReflectionGroupType)

Return whether the represented complex linear group acts irreducibly.

In particular, the trivial representation is irreducible precisely in ambient
dimension one.
"""
function is_irreducible(G::ComplexReflectionGroupType)
  return (length(G.type) == 1 && fixed_space_dimension(G) == 0) ||
         (isempty(G.type) && fixed_space_dimension(G) == 1)
end

"""
    is_equivalent(G1::ComplexReflectionGroupType, G2::ComplexReflectionGroupType)

Return whether `G1` and `G2` describe the same `GL(V)`-conjugacy type of complex
reflection group.
"""
function is_equivalent(G1::ComplexReflectionGroupType, G2::ComplexReflectionGroupType)
  return G1 == G2
end

# Group order
function order(G::ComplexReflectionGroupType)
  ord = ZZ(1)

  # Orders of the exceptional groups (grabbed from Magma)
  excords = [ 24, 72, 48, 144, 96, 192, 288, 576, 48, 96, 144, 288, 600, 1200, 1800, 3600,
  360, 720, 240, 120, 336, 648, 1296, 2160, 1152, 7680, 14400, 46080, 155520,
  51840, 39191040, 51840, 2903040, 696729600 ]

  for t in G.type
    if isa(t, Int)
      i = t
      ord *= ZZ(excords[i-3]) #index shift by 3 because exceptionals start with G4
    else
      m = ZZ(t[1])
      p = ZZ(t[2])
      n = ZZ(t[3])
      ord *= div( (m^n*factorial(n)), p )
    end
  end
  return ord
end

"""
    rank(G::ComplexReflectionGroupType)

Return the dimension of the support of `G`.
"""
function rank(G::ComplexReflectionGroupType)
  N = ZZ(0)

  # Ranks of the exceptional groups (grabbed from Magma)
  excranks = [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 6, 6, 7, 8 ]

  for t in G.type
    if isa(t, Int)
      i = t
      N += ZZ(excranks[i-3]) #index shift by 3 because exceptionals start with G4
    else
      m = ZZ(t[1])
      p = ZZ(t[2])
      n = ZZ(t[3])
      if m == 1 && p == 1 #symmetric group case is of rank n-1
        N += n-1
      else
        N += n
      end
    end
  end
  return N
end

"""
    degree(G::ComplexReflectionGroupType)

Return the dimension of the ambient reflection representation.
"""
degree(G::ComplexReflectionGroupType) =
  Int(rank(G)) + fixed_space_dimension(G)

# Imprimitivity
function is_imprimitive(G::ComplexReflectionGroupType)
  @req is_irreducible(G) "Imprimitivity is defined only for irreducible types"

  is_trivial(G) && return false

  t = G.type[1]
  if isa(t, Int) #exceptionals are primitive
    return false
  end

  (m,p,n) = t
  if rank(G) == 1
    return false
  end

  # S_n^ref is primitive for n >= 5. The apparent n = 4 exception in some
  # statements is excluded by S_4^ref = G(2, 2, 3), which is imprimitive.
  if m == 1 && p == 1 && n >= 5
    return false
  end

  return true
end

function is_primitive(G::ComplexReflectionGroupType)
  return !is_imprimitive(G)
end

"""
    direct_product(X::Vector{ComplexReflectionGroupType})
    direct_product(G1::ComplexReflectionGroupType,
                   G2::ComplexReflectionGroupType)

Return the type of the block-diagonal direct product. Both the nontrivial
irreducible factors and the fixed-space dimensions are combined.
"""
function direct_product(X::Vector{ComplexReflectionGroupType})
  type = ComplexReflectionGroupComponent[]
  fixed_dimension = 0
  for G in X
    append!(type, G.type)
    fixed_dimension += fixed_space_dimension(G)
  end

  return ComplexReflectionGroupType(type; fixed_space_dimension=fixed_dimension)

end

function direct_product(G1::ComplexReflectionGroupType, G2::ComplexReflectionGroupType)
  return direct_product([G1,G2])
end

function number_of_reflections(G::ComplexReflectionGroupType)

  N = ZZ(0)

  # See Lehrer & Taylor (2009). TABLE D.3 (p 275)
  excnum = [ 8, 16, 14, 22, 18, 30, 34, 46, 12, 18, 28, 34, 48, 78, 88, 118, 40, 70, 30, 15, 21, 24, 33, 45, 24, 40, 60, 60, 80, 45, 126, 36, 63, 120 ]

  for t in G.type
    if isa(t, Int)
      i = t
      N += ZZ(excnum[i-3]) #index shift by 3 because exceptionals start with G4
    else
      m = ZZ(t[1])
      p = ZZ(t[2])
      n = ZZ(t[3])
      N += div(m*(n^2-n),2) + n*(div(m,p)-1)
      # Formula follows from Lehrer & Taylor (2009), Lemma 2.8.
    end
  end

  return N
end

function is_well_generated(G::ComplexReflectionGroupType)

  # Grabbed from Magma
  excwell = [ 4, 5, 6, 8, 9, 10, 14, 16, 17, 18, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37 ]

  for t in G.type
    if isa(t, Int)
      if !(t in excwell)
        return false
      end
    else
      # See Lehrer & Taylor (2009), Section 2.7 (p 35)
      m = t[1]
      p = t[2]
      n = t[3]
      if !(p == 1 || p == m)
        return false
      end
    end
  end

  return true
end

function degrees(G::ComplexReflectionGroupType)

  degrees = fill(ZZ(1), fixed_space_dimension(G))

  # Lehrer & Taylor (2009), Table D.3 (p 275)
  excdegrees =
  [
    [ 4, 6 ],
    [ 6, 12 ],
    [ 4, 12 ],
    [ 12, 12 ],
    [ 8, 12 ],
    [ 8, 24 ],
    [ 12, 24 ],
    [ 24, 24 ],
    [ 6, 8 ],
    [ 8, 12 ],
    [ 6, 24 ],
    [ 12, 24 ],
    [ 20, 30 ],
    [ 20, 60 ],
    [ 30, 60 ],
    [ 60, 60 ],
    [ 12, 30 ],
    [ 12, 60 ],
    [ 12, 20 ],
    [ 2, 6, 10 ],
    [ 4, 6, 14 ],
    [ 6, 9, 12 ],
    [ 6, 12, 18 ],
    [ 6, 12, 30 ],
    [ 2, 6, 8, 12 ],
    [ 4, 8, 12, 20 ],
    [ 2, 12, 20, 30 ],
    [ 8, 12, 20, 24 ],
    [ 12, 18, 24, 30 ],
    [ 4, 6, 10, 12, 18 ],
    [ 6, 12, 18, 24, 30, 42 ],
    [ 2, 5, 6, 8, 9, 12 ],
    [ 2, 6, 8, 10, 12, 14, 18 ],
    [ 2, 8, 12, 14, 18, 20, 24, 30 ]
  ]

  for t in G.type
    if isa(t, Int)
      i = t
      append!(degrees, convert(Vector{ZZRingElem}, excdegrees[i-3]) )
    else
      # Lehrer & Taylor (2009), Appendix D.2 (p 274)
      (m,p,n) = t
      if m == 1 && p == 1 # symmetric group case
        append!(degrees, convert(Vector{ZZRingElem}, collect(2:n) ))
      elseif n == 1 # cyclic group case
        append!(degrees, convert(Vector{ZZRingElem}, [m]) )
      else
        append!(degrees, sort([[ZZ(m)*i for i=1:n-1] ; [div(ZZ(n)*ZZ(m),ZZ(p))]]))
      end
    end
  end

  return degrees

end

function codegrees(G::ComplexReflectionGroupType)

  codegrees = fill(ZZ(-1), fixed_space_dimension(G))

  # Lehrer & Taylor (2009), Table D.3 (p 275)
  exccodegrees =
  [
    [ 0, 2 ],  #4
    [ 0, 6 ],  #5
    [ 0, 8 ],  #6
    [ 0, 12 ], #7
    [ 0, 4 ],  #8
    [ 0, 16 ], #9
    [ 0, 12 ], #10
    [ 0, 24 ], #11
    [ 0, 10 ], #12
    [ 0, 16 ], #13
    [ 0, 18 ], #14
    [ 0, 24 ], #15
    [ 0, 10 ], #16
    [ 0, 40 ], #17
    [ 0, 30 ], #18
    [ 0, 60 ], #19
    [ 0, 18 ], #20
    [ 0, 48 ], #21
    [ 0, 28 ], #22
    [ 0, 4, 8 ],   #23
    [ 0, 8, 10 ],  #24
    [ 0, 3, 6 ],   #25
    [ 0, 6, 12 ],  #26
    [ 0, 18, 24 ], #27
    [ 0, 4, 6, 10 ],   #28
    [ 0, 8, 12, 16 ],  #29
    [ 0, 10, 18, 28 ], #30
    [ 0, 12, 16, 28 ], #31
    [ 0, 6, 12, 18 ],  #32
    [ 0, 6, 8, 12, 14 ], #33
    [ 0, 12, 18, 24, 30, 36 ], #34
    [ 0, 3, 4, 6, 7, 10 ], #35
    [ 0, 4, 6, 8, 10, 12, 16 ], #36
    [ 0, 6, 10, 12, 16, 18, 22, 28 ] #37
  ]

  for t in G.type
    if isa(t, Int)
      i = t
      append!(codegrees, convert(Vector{ZZRingElem}, exccodegrees[i-3]) )
    else
      # Lehrer & Taylor (2009), Appendix D.2 (p 274)
      (m,p,n) = t
      if m == 1 && p == 1 # symmetric group case
        append!(codegrees, convert(Vector{ZZRingElem}, collect(0:n-2) ))
      elseif n == 1 # cyclic group case
        append!(codegrees, convert(Vector{ZZRingElem}, [0]) )
      else
        if p != m
          append!(codegrees, [ZZ(m)*i for i=0:n-1])
        else
          append!(codegrees, sort([[ZZ(m)*i for i=0:n-2]; [ZZ(n-1)*ZZ(m)-ZZ(n)]]))
        end
      end
    end
  end

  return codegrees

end

function exponents(G::ComplexReflectionGroupType)

  # See Lehrer & Taylor (2009), p 51
  deg = degrees(G)
  return [ d-1 for d in deg ]

end

function coexponents(G::ComplexReflectionGroupType)

  # See Lehrer & Taylor (2009), p 205
  codeg = codegrees(G)
  return [ d+1 for d in codeg ]

end

function number_of_hyperplanes(G::ComplexReflectionGroupType)

  # See Geck & Malle (2003), Section 1.7

  coexp = coexponents(G)
  return sum(coexp)

end

function is_pseudo_real(G::ComplexReflectionGroupType)

  # Grabbed from Magma
  excpreal = [ 12, 13, 22, 23, 24, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37 ]

  for t in G.type
    if isa(t, Int)
      if !(t in excpreal)
        return false
      end
    else
      # See Lehrer & Taylor (2009). Exercise 2.15 (p 38)
      m = t[1]
      p = t[2]
      n = t[3]
      if !(m == p || ZZ(m) == 2*ZZ(p))
        return false
      end
    end
  end

  return true
end

function is_real(G::ComplexReflectionGroupType)

  excreal = [ 23, 28, 30, 35, 36, 37 ]

  for t in G.type
    if isa(t, Int)
      if !(t in excreal)
        return false
      end
    else
      m = t[1]
      p = t[2]
      n = t[3]
      # Only real groups are symmetric, dihedral, type B, and type D
      if !( (m == 1 && p == 1) || (m == p && n == 2) || (m == 2 && p == 1) || ( m == 2 && p == 2) )
        return false
      end
    end
  end

  return true

end

function is_coxeter_group(G::ComplexReflectionGroupType)
  return is_real(G)
end

function is_rational(G::ComplexReflectionGroupType)

  excrat = [ 28, 35, 36, 37 ]

  for t in G.type
    if isa(t, Int)
      if !(t in excrat)
        return false
      end
    else
      m = t[1]
      p = t[2]
      n = t[3]
      # Only rational groups are symmetric, dihedral order 12, type B, and type D
      if !( (m == 1 && p == 1) || (m == 6 && p == 6 && n == 2) || (m == 2 && p == 1) || ( m == 2 && p == 2) )
        return false
      end
    end
  end

  return true
end

function is_weyl_group(G::ComplexReflectionGroupType)
  return is_rational(G)
end

function coxeter_number(G::ComplexReflectionGroupType)

  @req is_irreducible(G) && !is_trivial(G) "The Coxeter number is defined only for irreducible nontrivial types"

  return div(number_of_reflections(G) + number_of_hyperplanes(G), rank(G))

end

# function is_shephard_group(G::ComplexReflectionGroupType)

#   # See Lehrer & Taylor (2009), Tables 6.1, 6.2, 6.3
#   excshep = [ 4, 5, 6, 8, 9, 10, 14, 16, 17, 18, 20, 21 ]

#   for t in G.type
#     if length(t) == 1
#       if !(t[1] in excshep)
#         return false
#       end
#     else
#       m = t[1]
#       p = t[2]
#       n = t[3]
#       if p != 1
#         return false
#       end
#     end
#   end

#   return true

# end

# function is_duality_group(G::ComplexReflectionGroupType)

# end

function is_spetsial(G::ComplexReflectionGroupType)

  # See Achar (2009)
  excspet = [ 4, 6, 8, 14, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37 ]

  for t in G.type
    if isa(t, Int)
      if !(t in excspet)
        return false
      end
    else
      m = t[1]
      p = t[2]
      n = t[3]
      if !( (p == 1) || (m == p) )
        return false
      end
    end
  end

  return true
end

function number_of_reflection_classes(G::ComplexReflectionGroupType)

  N = ZZ(0)

  # Grabbed from Magma
  excnum = [ 2, 4, 3, 5, 3, 4, 5, 6, 1, 2, 3, 4, 4, 5, 6, 7, 2, 3, 1, 1, 1, 2, 3, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1 ]

  for t in G.type
    if isa(t, Int)
      i = t
      N += ZZ(excnum[i-3]) #index shift by 3 because exceptionals start with G4
    else
      (m,p,n) = t

      # See Thiel (2014), Theorem 15.27
      if n == 1
        # Each nonidentity element of the cyclic group G(m, 1, 1) is a reflection
        # and forms its own conjugacy class.
        N += ZZ(m - 1)
      elseif n > 2 || (n == 2 && is_odd(p))
        N += div(ZZ(m),ZZ(p))
      else
        N += div(ZZ(m),ZZ(p)) + 1
      end
    end
  end

  return N

end

function center(G::ComplexReflectionGroupType)

  orders = ZZRingElem[]

  # The center of an irreducible group is always cyclic of order equal to the gcd
  # of the degrees.
  # See Thiel (2014), Proposition 13.32
  for C in components(G)
    push!(orders, gcd(degrees(C)))
  end

  return abelian_group(orders)

end

# function field_of_definition(G::ComplexReflectionGroupType)
#   # Reference is Geck & Malle (2006)

#   fields = Field[]

#   for t in G.type
#     if isa(t, Int)
#       i = t
#       if i == 4
#         K,z = cyclotomic_field(3)
#       elseif i == 5
#         K,z = cyclotomic_field(3)
#       elseif i == 6
#         K,z = cyclotomic_field(12)
#       elseif i == 7
#         K,z = cyclotomic_field(12)
#       elseif i == 8
#         K,z = cyclotomic_field(4)
#       elseif i == 9
#         K,z = cyclotomic_field(8)
#       elseif i == 10
#         K,z = cyclotomic_field(12)
#       elseif i == 11
#         K,z = cyclotomic_field(24)
#       elseif i == 12
#         R,X = polynomial_ring(QQ, "X")
#         f = X^2 + 2
#         K,z = number_field(f)
#       elseif i == 13
#         K,z = cyclotomic_field(8)
#       elseif i == 14
#         R,X = polynomial_ring(QQ, "X")
#         f1 = X^3 - 1
#         f2 = X^2 + 2
#         K,(a,b) = number_field([f1,f2])
#       else
#         nothing #Not finished yet!
#       end
#       push!(fields, K)
#     else
#       m = ZZ(t[1])
#       p = ZZ(t[2])
#       n = ZZ(t[3])
#     end
#   end

#   if length(fields) == 1
#     return fields[1]
#   else
#     println(fields)
#     # From Lars Göttgens
#     return reduce(compositum, fields; init=rationals_as_number_field()[1])[1]
#   end

# end
