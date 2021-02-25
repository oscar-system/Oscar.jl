<<<<<<< HEAD

# TODO : in this file, some functions for finite fields are defined just to make other files work,
# such as forms.jl, transform_form.jl, linear_conjugate.jl and linear_centralizer.jl
# TODO : functions in this file are only temporarily, and often inefficient.
# TODO: once similar working methods are defined in other files or packages (e.g. Hecke), 
# functions in this file are to be removed / moved / replaced
# TODO: when this happens, files mentioned above need to be modified too.


=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
import Hecke: evaluate, field_extension, FinField, FinFieldElem, PolyElem

export primitive_element

<<<<<<< HEAD
# changes the base ring of a polynomial ring into fq_nmod
=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
function _change_type(f::PolyElem{T}) where T <: FinFieldElem
   e,p = ispower(order(base_ring(f)))
   F = GF(Int(p),Int(e))[1]
   t = PolynomialRing(F,"t")[2]
   return sum([t^i*F(lift(coeff(f,i))) for i in 0:degree(f)])
end

<<<<<<< HEAD
# return a generator for the unit group of F = K[X] / (f), where K = base_ring(f)
=======
# return a generator for the unit group of F
>>>>>>> GAP: deal with matrices, vectors, finite field elements
function _centralizer(f::PolyElem{T}) where T <: FinFieldElem
  if typeof(f)!=fq_nmod_poly && typeof(f)!=fq_poly
     f = _change_type(f)
  end
<<<<<<< HEAD
  L, mL = field_extension(f)
=======
  K = base_ring(f)
  d = degree(K)
  e = degree(f)
  q = order(K)
  L, mL = field_extension(f)
#  m = divexact(order(L)-1, order(K)-1)
#  U, mU = unit_group(L, n_quo = Int(m))
>>>>>>> GAP: deal with matrices, vectors, finite field elements
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
<<<<<<< HEAD
function primitive_element(F::FinField)
   z = gen(F)
   isprime(order(F)) && return z
=======
function primitive_element(F::T) where T <: FinField
   z = gen(F)
   if isprime(order(F)) return z end
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   f = _centralizer(defining_polynomial(F))
   return sum([z^i*F(coeff(f,i)) for i in 0:degree(f)])
end

# TODO very bold discrete log, waiting for a better one. Don't try with large fields!!
<<<<<<< HEAD
# return g such that a^g = b
function _disc_log(a,b)
=======
function _disc_log(a,b)
   done=false
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   for g in 0:order(parent(a))
      if a^g==b
         done=true
         return g
      end
   end
<<<<<<< HEAD
   error("Second element is not a power of the first one")
=======
   @assert done "Second element is not a power of the first one"
>>>>>>> GAP: deal with matrices, vectors, finite field elements
end
