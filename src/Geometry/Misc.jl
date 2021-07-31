module Misc

using Oscar

export add_variables

function add_variables( R::MPolyRing, new_vars::Vector{String} )
  k = base_ring(R)
  old_vars = String.( symbols(R) )
  n = length( old_vars )
  vars = vcat( old_vars, new_vars )
  S, v = PolynomialRing( k, vars )
  @show S, v
  phi = AlgebraHomomorphism( R, S, gens(S)[1:n] )
  y = v[n+1:length(v)]
  return S, phi, y
end

end
