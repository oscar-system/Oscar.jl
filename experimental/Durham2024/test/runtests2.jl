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

#the canonical divisor of bl_IA3
bl_IA3_canonical_divisor = pullback(pr1, Oscar.canonical_divisor(IA3)) + (dim(IA3)-1)*E

# the cartier divisor of IA3 that correspdonds to X
X_ideal = image_ideal(inc_X)
X_ideal_sheaf = ideal_sheaf(codomain(pr1),  IA3, X_ideal)
X_divisor = EffectiveCartierDivisor(X_ideal_sheaf)

#the cartier divisor of bl_IA3 that correspdonds to bl_X
bl_X_divisor = strict_transform(pr1, X_divisor) 

# the canonical divisor of bl_X
bl_X_canonical_divisor = pullback(inc_bl_X, bl_IA3_canonical_divisor + bl_X_divisor)
end






