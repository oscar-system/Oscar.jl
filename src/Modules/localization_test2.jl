R, (x,y) = QQ["x", "y"]
U = MPolyComplementOfKPointIdeal(R, [1, -1])
L = Localization(U)

F = FreeMod(L, 2)
G = FreeMod(L, 1)

f = hom(F, G, [x//y^2*G[1], 1//y*G[1]])
K, inc_K = kernel(f)
for g in ambient_representatives_generators(K)
  @show f(g)
end

g = [x//y^2*F[1], 1//y*F[1]] 
q = [(x-1)*(y+1)*F[1]]
M = SubQuo(F, g, q)
coordinates(x*F[1], M)
