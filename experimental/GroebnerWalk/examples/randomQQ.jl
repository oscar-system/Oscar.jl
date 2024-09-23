#random polynomial system in 4 variables over Q, computed in M2 
using Oscar
R, (a, b, c, d)  = polynomial_ring(QQ, [:a, :b, :c, :d])


f1 = 15013d^3 + 5213*c*d^2 + 15893*c^2*d - 2837*b*d^2 - 10238*c^3 + 13538*b*c*d - 4578*a*d^2 + 5035*b*c^2 +
-9975*b^2*d - 9777*a*c*d - 2177*b^2*c + 9963*a*c^2 - 15569*a*b*d + 2513*b^3 - 3108*a*b*c - 9437*a^2*d - 5291*a*b^2 -
13120*a^2*c + 714*a^2*b + 13507*a^3

f2 = 2733*d^3 + 4541*c*d^2 + 4333*c^2*d - 6525*b*d^2 + 633*c^3 + 5406*b*c*d + 1762*a*d^2 + 11272*b*c^2 + 14879*b^2*d -
12968*a*c*d + 12580*b^2*c + 15050*a*c^2 + 15360*a*b*d - 2276*b^3 + 4773*a*b*c - 10499*a^2*d - 15025*a*b^2 -
9910*a^2*c - 3196*a^2*b + 1874*a^3

f3 = 11293*d^3 + 6929*c*d^2 - 13730*c^2*d - 7748*b*d^2 + 7572*c^3 - 9254*b*c*d + 15477*a*d^2 + 8782*b*c^2 + 7914*b^2*d +
5188*a*c*d + 15934*b^2*c - 10006*a*c^2 + 8830*a*b*d + 11842*b^3 + 7559*a*b*c + 11735*a^2*d - 7116*a*b^2 +
8771*a^2*c + 9251*a^2*b + 2099*a^3

I = ideal([f1, f2, f3])


Gdegrevlex = groebner_basis(I, complete_reduction = true)
set_verbosity_level(:groebner_walk, 1)

t_b = @elapsed Gb = groebner_basis(I, ordering = lex(R), complete_reduction = true) #8.6s

t_s = @elapsed Gs = groebner_walk(I, lex(R), algorithm =:standard) #11.3, one conversion

t_g = @elapsed Gg = groebner_walk(I, lex(R), algorithm =:generic) #14.4

t_p = @elapsed Gp = groebner_walk(I, lex(R), algorithm =:perturbed) #11.49

