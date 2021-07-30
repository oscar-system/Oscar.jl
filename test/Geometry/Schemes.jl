@testset "affine schemes"
  An = affine_space( QQ, 4 )

  R = An.R
  x = gens(R)

  f = x[1]^2 + x[2]^2 -1

  Y, phi= localize( An, f )

  y = gens( Y.R )
  g = y[1]^3*y[3]+y[2]^4 -5

  Z, psi = localize( Y, g )

  mu = compose( psi, phi )
end
