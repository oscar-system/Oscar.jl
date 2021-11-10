@testset "mpolyquo-localizations.jl" begin
  R, v = QQ["x", "y", "u", "v"]
  x = v[1]
  y = v[2] 
  u = v[3]
  v = v[4]
  f = x*v-y*u
  I = ideal(R, f)
  Q, p = quo(R, I)
  S = MPolyComplementOfKPointIdeal(R, [QQ(1), QQ(0), QQ(1), QQ(0)])
  T = MPolyComplementOfKPointIdeal(R, [QQ(0), QQ(0), QQ(0), QQ(0)])
  L = Localization(Q, S)
  a = L(x)
  b = L(y)
  c = L(u)
  d = L(v)
  
  R, v = ZZ["x", "y"]
  x = v[1]
  y = v[2] 
  f = (x^2 + y^2)
  T = MPolyComplementOfKPointIdeal(R, [ZZ(0), ZZ(0)])
  I = ideal(R, f)
  V = MPolyQuoLocalizedRing(R, I, T)

  S = R
  U = MPolyPowersOfElement(S, [f-1])
  J = ideal(S, zero(S))
  W = MPolyQuoLocalizedRing(S, J, U)

  h = MPolyQuoLocalizedRingHom(W, V, localized_ring(V).([x//(y-1), y//(x-5)]))
  reduce_fraction(h(W(f//(f-1)^9)))
end 
