
include("semiring.jl")
include("valuation.jl")
include("poly.jl")
include("initial.jl")
include("groebner_basis.jl")
include("groebner_polyhedron.jl")
include("points.jl")

# Decompose a tropical polynomial into parts that Polymake can eat.
# First function deals with the coefficients,
# Second function then deals with the entire polynomial.
function homogenize_and_convert_to_pm(t::Oscar.TropicalSemiringElem{S}) where S<:Union{typeof(max), typeof(min)}
   Add = S == typeof(max) ? Polymake.Max : Polymake.Min
   if isinf(t)
      return Polymake.TropicalNumber{Add}()
   else
      return Polymake.TropicalNumber{Add}(Polymake.new_rational_from_fmpq(data(t)))
   end
end

function homogenize_and_convert_to_pm(f::Oscar.MPolyRingElem{Oscar.TropicalSemiringElem{S}}) where S<:Union{typeof(max), typeof(min)}
   Add = S == typeof(max) ? Polymake.Max : Polymake.Min
   coeffs = (Polymake.TropicalNumber{Add})[]
   td = total_degree(f)
   exps = (Vector{Int})[]
   for term in terms(f)
      push!(coeffs, homogenize_and_convert_to_pm(leading_coefficient(term)))
      exp = leading_exponent_vector(term)
      prepend!(exp, td-sum(exp))
      push!(exps, exp)
   end
   exps = matrix(ZZ, exps)
   coeffs = Polymake.Vector{Polymake.TropicalNumber{Add, Polymake.Rational}}(coeffs)
   return coeffs, exps
end



###
# Allow gcd of vectors of univariate rational polynomials
# to make their handling similar to that of integers
###
gcd(F::Vector{QQPolyRingElem}) = reduce(gcd, F)
gcd(F::QQPolyRingElem...) = reduce(gcd, F)



include("variety_supertype.jl")
include("variety.jl")
include("hypersurface.jl")
include("curve.jl")
include("linear_space.jl")
include("groebner_fan.jl")
