@testset "mpolyquo-localizations.jl" begin
  R, var = QQ["x", "y", "u", "v"]
  x = var[1]
  y = var[2] 
  u = var[3]
  v = var[4]
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
end 
