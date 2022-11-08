K = GF(2)
Kx,(x1,x2,x3,x4) = PolynomialRing(K,4);
inI = ideal([x2 + x3 + x4, x1 + x3])
println(inI) # same print as input for tropical_link
println(base_ring(inI)) # same print as input for tropical_link

x = gens(Kx)
inI1 = inI + ideal(Kx,[x[1]+1])

singular_assure(inI1)
singularIdeal = inI1.gens.S
singularRing = base_ring(singularIdeal)
singularIdeal = Singular.satstd(singularIdeal,Singular.MaximalIdeal(singularRing,1))
inI1 = ideal(Kx,singularIdeal)

L,t = RationalFunctionField(K,"t")
Lx,x = PolynomialRing(L,symbols(Kx))
inI1 = ideal(Lx,[change_base_ring(L,g) for g in gens(inI1)])

hyperplanes = [x[i]-t for i in [2,3,4]]
append!(hyperplanes,[t*x[i]-1 for i in [2,3,4]])
for hyperplane in hyperplanes
  display(inI1)
  inI0 = inI1+ideal(Lx,[hyperplane])
  display(inI0)
end














Kx,(x1,x2,x3,x4) = PolynomialRing(QQ,4);
p = 2;
I = ideal([x1-p*x2+(p+1)*x3,3*x2-p^2*x3+(p^2+1)*x4]);
w = Int[-1, 1, -1, 1]

###
# Step 1: compute a tropical Groebner basis via
#   groebner_basis(I,val,w), val=2-adic valuation
###

# Step 1.1: running simulate_valuation
G = gens(I)
Rtx,tx = PolynomialRing(ZZ,vcat([:t],symbols(parent(G[1]))))
vvG = [p-tx[1]]
for f in G
  fRtx = MPolyBuildCtx(Rtx)
  for (cK,expvKx) = zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    cR = numerator(cK)
    expvRtx = vcat([0],expvKx)
    push_term!(fRtx,cR,expvRtx)
  end
  push!(vvG,tighten_simulation(finish(fRtx),val))
end

# Step 1.2: compute standard basis in Singular
vvw = -w
pushfirst!(vvw,-1)
S,_ = Singular.PolynomialRing(singular_ring(base_ring(Rtx)), map(string, Nemo.symbols(Rtx)), ordering = Singular.ordering_a(vvw)*Singular.ordering_dp())
SI = Singular.Ideal(S, [S(g) for g in vvG])
vvGB = Singular.gens(Singular.satstd(SI,Singular.MaximalIdeal(S,1)))
vvGB = [Rtx(g) for g in vvGB]

# Step 1.3: desimulate_valuation
Rx = parent(vvGB[1])
x = copy(symbols(Rx))
popfirst!(x)

K = QQ
Kx,_ = PolynomialRing(K,x)

GB = []
for i = 2:3
  vvg = evaluate(vvGB[i],[1],[ZZ(p)])
  g = MPolyBuildCtx(Kx)
  for (c, expvRtx) = Base.Iterators.zip(AbstractAlgebra.coefficients(vvg), AbstractAlgebra.exponent_vectors(vvg))
    expvKx = expvRtx
    popfirst!(expvKx)
    push_term!(g,K(c),expvKx)
  end
  append!(GB,[finish(g)])
end


###
# Step 2: construct initial ideal from the tropical groebner basis
###
f = GB[1]
kx, x = PolynomialRing(GF(2),[repr(x) for x in gens(parent(f))])

initialf = MPolyBuildCtx(kx)
