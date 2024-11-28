@testset "canonical divisor of a Kleinian singularity" begin
# We calculate the canonical divisor of the blow up of a Kleinian singularity in this document.
R, (x, y, z) = QQ[:x, :y, :z] 
I = ideal(R, [x^3 - y*z])  

IA3 = spec(R)
X, inc_X = sub(IA3, I) #X is the A2 Kleinian singularity
I_sing = ideal(R, gens(R))

pr1 = blow_up(IA3, I_sing) #blow up the affine space IA3
bl_IA3 = domain(pr1) 
E = exceptional_divisor(pr1) 
bl_X, inc_bl_X, pr2 = strict_transform(pr1, inc_X) #bl_X is the blow up of X

EX = pullback(inc_bl_X, E)
EX_weil = weil_divisor(EX)
EX_weil_decom = irreducible_decomposition(EX_weil)
C1 = components(WeilDivisor, EX_weil_decom)[1]
C2 = components(WeilDivisor, EX_weil_decom)[2]

bl_X_canonical_divisor = weil_divisor(bl_X, ZZ) 
Z11 = Oscar.self_intersection_via_adjunction(bl_X_canonical_divisor, C1, 0)
Z22 = Oscar.self_intersection_via_adjunction(bl_X_canonical_divisor, C2, 0)
Z12 = intersect(C1, C2)
Z21 = intersect(C2, C1)
Intersection_matrix = ZZ[Z11 Z12; Z21 Z22]
end






