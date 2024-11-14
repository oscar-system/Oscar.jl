Qx, x = QQ["x"]
K, a = number_field(x^9 - 3*x^8 + x^6 + 15*x^5 - 13*x^4 -
                    3*x^3 + 4*x - 1, "a")
G, C = galois_group(K)
subsG = subgroups(G)
H = first(H for H in subsG if order(H) == 27)
k, = simplify(fixed_field(C, H))
# output
(Number field of degree 8 over QQ, Map: k -> number field)
