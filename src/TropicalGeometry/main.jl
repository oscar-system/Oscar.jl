include("TropicalNumbers/numbers.jl")

# Temporarily we will turn tropical polynomials into strings. This will be
# removed once Polymake.jl wraps its tropical polynomials and tropical numbers
#
# Warning: This function ignores all boundary cases!
function tropical_polynomial_to_polymake(f)
    convention = fun(base_ring(f))
    result = ""
    if convention == min
        result *= "min("
    else
        result *= "max("
    end
    td = total_degree(f)
    for i in 1:length(f)
        result *= repr(coeff(f,i).data)
        e = exponent_vector(f,i)
        if td - sum(e) != 0
            result *= "+"
            result *= repr(td-sum(e))
            result *= "x(0)"
        end
        if !iszero(e)
            for j in 1:length(e)
                if !iszero(e[j])
                    result *= "+"
                    result *= repr(e[j])
                    result *= "*x("
                    result *= repr(j)
                    result *= ")"
                end
            end
        end
        if i != length(f)
            result *= ","
        end
    end
    result *= ")"
end


function +(a::AbstractAlgebra.Generic.MPoly{T}, b::AbstractAlgebra.Generic.MPoly{T}) where {T <: RingElement}
   N = size(a.exps, 1)
   par = parent(a)
   r = par()
   fit!(r, length(a) + length(b))
   i = 1
   j = 1
   k = 1
   while i <= length(a) && j <= length(b)
      cmpexp = AbstractAlgebra.Generic.monomial_cmp(a.exps, i, b.exps, j, N, par, UInt(0))
      if cmpexp > 0
         r.coeffs[k] = a.coeffs[i]
         monomial_set!(r.exps, k, a.exps, i, N)
         i += 1
      elseif cmpexp == 0
         c = a.coeffs[i] + b.coeffs[j]
         if !iszero(c)
            r.coeffs[k] = c
            AbstractAlgebra.Generic.monomial_set!(r.exps, k, a.exps, i, N)
         else
            k -= 1
         end
         i += 1
         j += 1
      else
         r.coeffs[k] = b.coeffs[j]
         AbstractAlgebra.Generic.monomial_set!(r.exps, k, b.exps, j, N)
         j += 1
      end
      k += 1
   end
   while i <= length(a)
      r.coeffs[k] = a.coeffs[i]
      AbstractAlgebra.Generic.monomial_set!(r.exps, k, a.exps, i, N)
      i += 1
      k += 1
   end
   while j <= length(b)
      r.coeffs[k] = b.coeffs[j]
      AbstractAlgebra.Generic.monomial_set!(r.exps, k, b.exps, j, N)
      j += 1
      k += 1
   end
   r.length = k - 1
   return r
end


include("TropicalVariety/constructors.jl")
include("TropicalVariety/properties.jl")
