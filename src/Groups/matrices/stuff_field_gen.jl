import Hecke: evaluate, field_extension, FinField, FinFieldElem, PolyElem

export primitive_element

function _change_type(f::PolyElem{T}) where T <: FinFieldElem
   e,p = ispower(order(base_ring(f)))
   F = GF(Int(p),Int(e))[1]
   t = PolynomialRing(F,"t")[2]
   return sum([t^i*F(lift(coeff(f,i))) for i in 0:degree(f)])
end

# return a generator for the unit group of F
function _centralizer(f::PolyElem{T}) where T <: FinFieldElem
  if typeof(f)!=fq_nmod_poly && typeof(f)!=fq_poly
     f = _change_type(f)
  end
  K = base_ring(f)
  d = degree(K)
  e = degree(f)
  q = order(K)
  L, mL = field_extension(f)
#  m = divexact(order(L)-1, order(K)-1)
#  U, mU = unit_group(L, n_quo = Int(m))
  U, mU = unit_group(L)
  g = mU(U[1])
  return mL\g
end


# return a primitive element of F, i.e. a group generator for F*
# TODO: are there faster procedures?
"""
    primitive_element(F::FinField)
Return a generator of the multiplicative group of `F`.
"""
function primitive_element(F::T) where T <: FinField
   z = gen(F)
   if isprime(order(F)) return z end
   f = _centralizer(defining_polynomial(F))
   return sum([z^i*F(coeff(f,i)) for i in 0:degree(f)])
end

# TODO very bold discrete log, waiting for a better one. Don't try with large fields!!
function _disc_log(a,b)
   done=false
   for g in 0:order(parent(a))
      if a^g==b
         done=true
         return g
      end
   end
   @assert done "Second element is not a power of the first one"
end
