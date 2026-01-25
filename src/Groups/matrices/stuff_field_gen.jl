
# TODO : in this file, some functions for finite fields are defined just to make other files work,
# such as forms.jl, transform_form.jl, linear_conjugate.jl and linear_centralizer.jl
# TODO : functions in this file are only temporarily, and often inefficient.
# TODO: once similar working methods are defined in other files or packages (e.g. Hecke), 
# functions in this file are to be removed / moved / replaced
# TODO: when this happens, files mentioned above need to be modified too.

# changes the base ring of a polynomial ring into fqPolyRepFieldElem
function _change_type(f::PolyRingElem{T}) where T <: FinFieldElem
   e,p = is_perfect_power_with_data(order(base_ring(f)))
   F = GF(Int(p),Int(e))
   t = polynomial_ring(F,:t; cached=false)[2]
   return sum([t^i*F(lift(coeff(f,i))) for i in 0:degree(f)])
end

function _change_type(f::PolyRingElem{<: FqFieldElem})
   e,p = is_perfect_power_with_data(order(base_ring(f)))
   F = GF(Int(p),Int(e))
   t = polynomial_ring(F,:t; cached=false)[2]
   return sum([t^i*F(lift(ZZ, coeff(f,i))) for i in 0:degree(f)])
end

#TODO: _change_type not simpler???
# -> _change_type used only in the current file (and in a testfile ...),
#                      only in _centralizer

# if f in F[x] and z is a root of f in the splitting field of f over F,
# then return a polynomial g such that g(z) is a generator for the unit group of F(z)
# return a generator for the unit group of F = K[X] / (f), where K = base_ring(f)
function _centralizer(f::PolyRingElem{T}) where T <: FinFieldElem
  if !(f isa Union{fqPolyRepPolyRingElem, FqPolyRepPolyRingElem, FqPolyRingElem})
     _f = _change_type(f)
     return _centralizer(_f)
  end

  L, mL = field_extension(f)
  U, mU = unit_group(L)
  g = mU(U[1])
  return mL\g
end
