export radical, primary_decomposition, minimal_primes, equidimensional_decomposition_weak,
          equidimensional_decomposition_radical, equidimensional_hull,
          equidimensional_hull_radical
# ideal quotient #######################################################
@doc Markdown.doc"""
    quotient(I::MPolyIdeal, J::MPolyIdeal)
    
Returns the ideal quotient of `I` by `J`.
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
###################################################

# primary decomposition #######################################################

#######################################################
@doc Markdown.doc"""
    radical(I::MPolyIdeal)
    
If the base ring of `I` is a polynomial ring over a field or over the integers, return
the radical of `I`.
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

Compute a primary decomposition of the ideal `I`. If the base ring of `I` is a polynomial
ring over a field, the algorithm of Gianni-Trager-Zacharias is used by default. Alternatively,
the Shimoyama-Yokoyama algorithm can be used by specifying `alg=:SY`.  For polynomial
rings over the integers, the algorithm proceeds as suggested by Pfister, Sadiq, and Steidel.
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

Return an array of the minimal associated prime ideals of `I`.
If `I` is the unit ideal, `[ideal(1)]` is returned.
If the base ring of `I` is a polynomial ring over a field, the algorithm of
Gianni-Trager-Zacharias is used by default and characteristic sets may be
used by specifying `alg=:charSets`.
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

Return an array of equidimensional ideals where the last element is the
equidimensional hull of `I`, that is, the intersection of the primary
components of `I` of maximal dimension, and each of the previous elements
is an ideal of lower dimension whose associated primes are exactly the associated
primes of `I` of that dimension.
If `I` is the unit ideal, `[ideal(1)]` is returned.
"""

function equidimensional_decomposition_weak(I::MPolyIdeal)
  R = base_ring(I)
  singular_assure(I)
  l = Singular.LibPrimdec.equidim(I.gens.Sx, I.gens.S)
  return [ideal(R, i) for i in l]
end

@doc Markdown.doc"""
    equidimensional_decomposition_radical(I::MPolyIdeal)

Return an array of equidimensional radical ideals increasingly ordered by dimension.
For each dimension, the returned radical ideal is the intersection of the associated primes 
of `I` of that dimension. If `I` is the unit ideal, `[ideal(1)]` is returned.
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
minimal height.
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
Return the intersection of associated primes of `I` of maximal dimension.
"""

function equidimensional_hull_radical(I::MPolyIdeal)
  R = base_ring(I)
  singular_assure(I)
  i = Singular.LibPrimdec.equiRadical(I.gens.Sx, I.gens.S)
  return ideal(R, i)
end

