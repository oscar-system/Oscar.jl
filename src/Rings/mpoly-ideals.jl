export saturation, quotient, elimination
export radical, primary_decomposition, minimal_primes, equidimensional_decomposition_weak,
          equidimensional_decomposition_radical, equidimensional_hull,
          equidimensional_hull_radical
export issubset, iscontained, isprime, isprimary
export ngens, gens

# elementary operations #######################################################

@doc Markdown.doc"""
    :^(I::MPolyIdeal, m::Int)

Returns the m-th power of `I`. 
"""
function Base.:^(I::MPolyIdeal, m::Int)
  singular_assure(I)
  return MPolyIdeal(I.gens.Ox, I.gens.S^m)
end

@doc Markdown.doc"""
    :+(I::MPolyIdeal, J::MPolyIdeal)

Returns the sum of `I` and `J`. 
"""
function Base.:+(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, I.gens.S + J.gens.S)
end
Base.:-(I::MPolyIdeal, J::MPolyIdeal) = I+J

@doc Markdown.doc"""
    :*(I::MPolyIdeal, J::MPolyIdeal)

Returns the product of `I` and `J`. 
"""
function Base.:*(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, I.gens.S * J.gens.S)
end

#######################################################

# ideal intersection #######################################################
@doc Markdown.doc"""
    intersect(I::MPolyIdeal, Js::MPolyIdeal...)

Returns the intersection of two or more ideals.
"""
function Base.intersect(I::MPolyIdeal, Js::MPolyIdeal...)
  singular_assure(I)
  si = I.gens.S
  for J in Js
    singular_assure(J)
    si = Singular.intersection(si, J.gens.S)
  end
  return MPolyIdeal(I.gens.Ox, si)
end

#######################################################

# ideal quotient #######################################################
@doc Markdown.doc"""
    quotient(I::MPolyIdeal, J::MPolyIdeal)
    
Returns the ideal quotient of `I` by `J`. Alternatively, use `I:J`. 
"""
function quotient(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, Singular.quotient(I.gens.S, J.gens.S))
end

(::Colon)(I::MPolyIdeal, J::MPolyIdeal) = quotient(I, J)

#######################################################

# saturation #######################################################
@doc Markdown.doc"""
    saturation(I::MPolyIdeal, J::MPolyIdeal)
    
Returns the saturation of `I` with respect to `J`.
"""
function saturation(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return MPolyIdeal(I.gens.Ox, Singular.saturation(I.gens.S, J.gens.S))
end
#######################################################

# elimination #######################################################
@doc Markdown.doc"""
    eliminate(I::MPolyIdeal, polys::Array{MPolyElem, 1})

Given a list of polynomials which are variables, these variables are eliminated from `I`. That is,
the function returns the ideal of all polynomials in `I` which only depend on the remaining variables.
"""
function eliminate(I::MPolyIdeal, l::Array{<:MPolyElem, 1})
  singular_assure(I)
  B = BiPolyArray(l)
  S = base_ring(I.gens.S)
  s = Singular.eliminate(I.gens.S, [S(x) for x = l]...)
  return MPolyIdeal(base_ring(I), s)
end

@doc Markdown.doc"""
    eliminate(I::MPolyIdeal, polys::AbstractArray{Int, 1})

Given a list of indices which specify variables, these variables are eliminated from `I`. That is,
the function returns the ideal of all polynomials in `I` which only depend on the remaining variables.
"""
function eliminate(I::MPolyIdeal, l::AbstractArray{Int, 1})
  R = base_ring(I)
  return eliminate(I, [gen(R, i) for i=l])
end

### todo: wenn schon GB bzgl. richtiger eliminationsordnung bekannt ...
### Frage: return MPolyIdeal(base_ring(I), s) ???

###################################################

# primary decomposition #######################################################

#######################################################
@doc Markdown.doc"""
    radical(I::MPolyIdeal)
    
Returns the radical of `I`. If the base ring of `I` is a polynomial
ring over a field, a combination of the algorithms of Krick and Logar 
(with modifications by Laplagne) and Kemper is used. For polynomial
rings over the integers, the algorithm proceeds as suggested by 
Pfister, Sadiq, and Steidel.
"""
function radical(I::MPolyIdeal)
  singular_assure(I)
  R = base_ring(I)
  if elem_type(base_ring(R)) <: FieldElement
  J = Singular.LibPrimdec.radical(I.gens.Sx, I.gens.S)
  elseif base_ring(I.gens.Sx) isa Singular.Integers
  J = Singular.LibPrimdecint.radicalZ(I.gens.Sx, I.gens.S)
  else
   error("not implemented for base ring")
  end
  return ideal(R, J)
end
#######################################################
@doc Markdown.doc"""
    primary_decomposition(I::MPolyIdeal)

Returns a primary decomposition of the ideal `I`. If `I` is the unit ideal, `[ideal(1)]` is returned.
If the base ring of `I` is a polynomial ring over a field, the algorithm of Gianni, Trager and Zacharias 
is used by default. Alternatively, the algorithm by Shimoyama and Yokoyama can be used by specifying 
`alg=:SY`.  For polynomial rings over the integers, the algorithm proceeds as suggested by 
Pfister, Sadiq, and Steidel.
"""
function primary_decomposition(I::MPolyIdeal; alg=:GTZ)
  R = base_ring(I)
  singular_assure(I)
  if elem_type(base_ring(R)) <: FieldElement
    if alg == :GTZ
      L = Singular.LibPrimdec.primdecGTZ(I.gens.Sx, I.gens.S)
    elseif alg == :SY
      L = Singular.LibPrimdec.primdecSY(I.gens.Sx, I.gens.S)
    else
      error("algorithm invalid")
    end
  elseif base_ring(I.gens.Sx) isa Singular.Integers
    L = Singular.LibPrimdecint.primdecZ(I.gens.Sx, I.gens.S)
  else
    error("base ring not implemented")
  end
  return [(ideal(R, q[1]), ideal(R, q[2])) for q in L]
end
#######################################################y
@doc Markdown.doc"""
    minimal_primes(I::MPolyIdeal; alg=:GTZ)

Returns an array of the minimal associated prime ideals of `I`.
If `I` is the unit ideal, `[ideal(1)]` is returned.
If the base ring of `I` is a polynomial ring over a field, the algorithm of
Gianni-Trager-Zacharias is used by default and characteristic sets may be
used by specifying `alg=:charSets`. For polynomial rings over the integers, 
the algorithm proceeds as suggested by Pfister, Sadiq, and Steidel.
"""
function minimal_primes(I::MPolyIdeal; alg = :GTZ)
  R = base_ring(I)
  singular_assure(I)
  if elem_type(base_ring(R)) <: FieldElement
    if alg == :GTZ
      l = Singular.LibPrimdec.minAssGTZ(I.gens.Sx, I.gens.S)
    elseif alg == :charSets
      l = Singular.LibPrimdec.minAssChar(I.gens.Sx, I.gens.S)
    else
      error("algorithm invalid")
    end
  elseif base_ring(I.gens.Sx) isa Singular.Integers
    l = Singular.LibPrimdecint.minAssZ(I.gens.Sx, I.gens.S)
  else
    error("base ring not implemented")
  end
  return [ideal(R, i) for i in l]
end
#######################################################
@doc Markdown.doc"""
    equidimensional_decomposition_weak(I::MPolyIdeal)

Returns an array of equidimensional ideals where the last element is the
equidimensional hull of `I`, that is, the intersection of the primary
components of `I` of maximal dimension. Each of the previous elements
is an ideal of lower dimension whose associated primes are exactly the associated
primes of `I` of that dimension. If `I` is the unit ideal, `[ideal(1)]` is returned.
Uses ideas of Eisenbud, Huneke, and Vasconcelos.
"""
function equidimensional_decomposition_weak(I::MPolyIdeal)
  R = base_ring(I)
  singular_assure(I)
  l = Singular.LibPrimdec.equidim(I.gens.Sx, I.gens.S)
  return [ideal(R, i) for i in l]
end

@doc Markdown.doc"""
    equidimensional_decomposition_radical(I::MPolyIdeal)

Returns an array of equidimensional radical ideals increasingly ordered by dimension.
For each dimension, the returned radical ideal is the intersection of the associated primes 
of `I` of that dimension. If `I` is the unit ideal, `[ideal(1)]` is returned.
Uses a combination of the algorithms of Krick and Logar (with modifications by Laplagne) and Kemper.
"""
function equidimensional_decomposition_radical(I::MPolyIdeal)
  R = base_ring(I)
  singular_assure(I)
  l = Singular.LibPrimdec.prepareAss(I.gens.Sx, I.gens.S)
  return [ideal(R, i) for i in l]
end
#######################################################
@doc Markdown.doc"""
    equidimensional_hull(I::MPolyIdeal)

If the base ring of `I` is a polynomial ring over a field, return the intersection
of the primary components of `I` of maximal dimension. In the case of polynomials
over the integers, return the intersection of the primary components of I of
minimal height.  If `I` is the unit ideal, `[ideal(1)]` is returned. 
For polynomial rings over a field, the algorithm relies on ideas as used by
Gianni, Trager, and Zacharias or Krick and Logar. For polynomial rings over the integers, 
the algorithm proceeds as suggested by Pfister, Sadiq, and Steidel.
"""
function equidimensional_hull(I::MPolyIdeal)
  R = base_ring(I)
  singular_assure(I)
  if elem_type(base_ring(R)) <: FieldElement
    i = Singular.LibPrimdec.equidimMax(I.gens.Sx, I.gens.S)
  elseif base_ring(I.gens.Sx) isa Singular.Integers
    i = Singular.LibPrimdecint.equidimZ(I.gens.Sx, I.gens.S)
  else
    error("base ring not implemented")
  end
  return ideal(R, i)
end
#######################################################
@doc Markdown.doc"""
    equidimensional_hull_radical(I::MPolyIdeal)

Returns the intersection of the associated primes of `I` of maximal dimension.
If `I` is the unit ideal, `[ideal(1)]` is returned. 
Uses a combination of the algorithms of Krick and Logar 
(with modifications by Laplagne) and Kemper. 
"""
function equidimensional_hull_radical(I::MPolyIdeal)
  R = base_ring(I)
  singular_assure(I)
  i = Singular.LibPrimdec.equiRadical(I.gens.Sx, I.gens.S)
  return ideal(R, i)
end

#######################################################
@doc Markdown.doc"""
    :(==)(I::MPolyIdeal, J::MPolyIdeal)

Returns `true` if `I=J`, `false` otherwise.
"""
function Base.:(==)(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return Singular.equal(I.gens.S, J.gens.S)
end

### todo: wenn schon GB's  bekannt ...

#######################################################
@doc Markdown.doc"""
    issubset(I::MPolyIdeal, J::MPolyIdeal)

Returns `true` if `I` is contained in `J`, `false` otherwise.
"""
function Base.issubset(I::MPolyIdeal, J::MPolyIdeal)
  singular_assure(I)
  singular_assure(J)
  return Singular.contains(J.gens.S, I.gens.S)
end

### todo: wenn schon GB's  bekannt ...

#######################################################
@doc Markdown.doc"""
    iscontained(f::MPolyElem, J::MPolyIdeal)

Returns `true` if `f` is contained in `J` and `false`, otherwise.
"""
function iscontained(f::MPolyElem, J::MPolyIdeal)
  return issubset(ideal([f]), J)
end

################################################################################
@doc Markdown.doc"""
    isprime(I::MPolyIdeal)

Returns `true` if the ideal `I` is prime, `false` otherwise. Proceeds by computing a primary decomposition.
"""
function isprime(I::MPolyIdeal)
  D = primary_decomposition(I)
  return length(D) == 1 && issubset(D[1][2], D[1][1])
end

################################################################################
@doc Markdown.doc"""
    isprimary(I::MPolyIdeal)

Return `true` if the ideal `I` is primary, `false` otherwise. Proceeds by computing a primary decomposition.
"""
function isprimary(I::MPolyIdeal)
  D = primary_decomposition(I)
  return length(D) == 1
end

#######################################################
@doc Markdown.doc"""
    ngens(I::MPolyIdeal)
Returns the number of generators of `I`.
"""
function ngens(I::MPolyIdeal)
  return length(I.gens)
end

#######################################################
@doc Markdown.doc"""
    gens(I::MPolyIdeal)

Returns the generators of `I` as an array of multivariate polynomials.
"""
function gens(I::MPolyIdeal)
  return [I.gens[Val(:O), i] for i=1:ngens(I)]
end

gen(I::MPolyIdeal, i::Int) = I.gens[Val(:O), i]
