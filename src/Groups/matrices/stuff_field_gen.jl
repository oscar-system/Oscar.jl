
# TODO : in this file, some functions for finite fields are defined just to make other files work,
# such as forms.jl, transform_form.jl, linear_conjugate.jl and linear_centralizer.jl
# TODO : functions in this file are only temporarily, and often inefficient.
# TODO: once similar working methods are defined in other files or packages (e.g. Hecke), 
# functions in this file are to be removed / moved / replaced
# TODO: when this happens, files mentioned above need to be modified too.

import Hecke: evaluate, field_extension, FinField, FinFieldElem, PolyElem

export primitive_element

# changes the base ring of a polynomial ring into fq_nmod
function _change_type(f::PolyElem{T}) where T <: FinFieldElem
   e,p = ispower(order(base_ring(f)))
   F = GF(Int(p),Int(e))[1]
   t = PolynomialRing(F,"t")[2]
   return sum([t^i*F(lift(coeff(f,i))) for i in 0:degree(f)])
end


# return a generator for the unit group of F = K[X] / (f), where K = base_ring(f)
function _centralizer(f::PolyElem{T}) where T <: FinFieldElem
  if typeof(f)!=fq_nmod_poly && typeof(f)!=fq_poly
     f = _change_type(f)
  end
  L, mL = field_extension(f)
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
function primitive_element(F::FinField)
   z = gen(F)
   isprime(order(F)) && return z
   f = _centralizer(defining_polynomial(F))
   return sum([z^i*F(coeff(f,i)) for i in 0:degree(f)])
end

# TODO very bold discrete log, waiting for a better one. Don't try with large fields!!
# return g such that a^g = b
function _disc_log(a,b)
   for g in 0:order(parent(a))
      if a^g==b
         return g
      end
   end
   error("Second element is not a power of the first one")
end
