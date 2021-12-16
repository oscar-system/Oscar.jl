R, (x,y) = QQ["x", "y"]
S, (t) = PolynomialRing(QQ, ["t"])
t = t[1]

A = Spec(R)
B = Spec(S)
f = x^2+x^3 -y^2
C = subscheme(A, f)

phi = maximal_extension(C, B, [x//y])
psi = maximal_extension(B, C, [(1-t^2)//t^2, (1-t^2)//t^3])
psi = restriction(psi, SpecOpen(B, [t*(1-t^2)]), domain(phi))
phi_res = restriction(phi, domain(phi), domain(psi))
#psi_res = restriction(psi, domain(psi), domain(phi))
p = compose(psi, phi)
q = compose(phi_res, psi)
@show restriction(p, domain(p), domain(p)) == identity_map(domain(p))
@show restriction(q, domain(q), domain(q)) == identity_map(domain(q))

