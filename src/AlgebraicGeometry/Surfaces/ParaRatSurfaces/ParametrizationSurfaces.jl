
##############################################################################
#
# Rational Parametrizations of rational surfaces
#
##############################################################################

@doc raw"""
    parametrization(X::AbsProjectiveVariety)

Given a smooth rational surface `X` which is linearly normal in the given embedding, return a rational parametrization of `X`.

!!! note
    The function does not check whether `X` is smooth. If you are uncertain, enter `is_smooth(X)` first.

!!! note
    The function does not check rationality. In fact, at current state, `OSCAR` does not offer a direct check for this.

!!! note
    The function makes use of the adjunction process. It returns an error message if the terminal object of the adjunction process is not the projective plane. See the OSCAR documentation for information on the adjunction process.

# Examples
```jldoctest
julia> X = bordiga()
Projective variety
  in projective 4-space over GF(31991) with coordinates [x, y, z, u, v]
defined by ideal with 4 generators

julia> dim(X)
2

julia> codim(X)
2

julia> phi = parametrization(X);

julia> domain(phi)
Multivariate polynomial ring in 5 variables over GF(31991) graded by
  x -> [1]
  y -> [1]
  z -> [1]
  u -> [1]
  v -> [1]

julia> codomain(phi)
Multivariate polynomial ring in 3 variables over GF(31991) graded by
  z[1] -> [1]
  z[2] -> [1]
  z[3] -> [1]

julia> [degree(phi(x)) for x in gens(ambient_coordinate_ring(X))]
5-element Vector{FinGenAbGroupElem}:
 [4]
 [4]
 [4]
 [4]
 [4]
```
"""
function parametrization(X::AbsProjectiveVariety)
   L = adjunction_process(X);
   @req is_zero(defining_ideal(L[4])) "The terminal object of the adjunction process is not the projective plane."
   RX = ambient_coordinate_ring(X)
   if L[2] == []
      return hom(RX, RX, gens(RX))
   end
   n = length(L[2])
   RXterminal = base_ring(L[2][n])
   H = gens(RXterminal)
   while n > 0    
     S = base_ring(L[2][n])
     phi = hom(S, RXterminal, H)
     M =  syz(transpose(map_entries(phi, L[2][n])))
     H = M[1, :]
     n = n-1  
   end
   return hom(RX, RXterminal, H)
end
