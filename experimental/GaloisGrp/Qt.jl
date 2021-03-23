#= musings about Qt - before I forget
 
 f in Q(t)[x] - or better Z[t][x]

 for all t (outside some bad points), the roots are power series over C
 
   R_i(z) =: R(z) = sum a_n (z-t)^n

 By Taylor, Cauchy, Dan and his partner, 
 
 a_n = 1/2/pi/i int_{|z| = r} R(z)/(z-t)^(n+1) dz

 Now f(z, R(z)) = 0

 f = sum f_i(t) x^i

 By standart bounds on polynomials (assume f monic)

 roots of f(z)(x) are bounded in abs. value by |f_i(z)| + 1 
 (or even 2*|f_i(z)|^(1/(n-i)) or so, 
 https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots
 Donald Knuth)

 if |f|_infty = B (largest coeff in Z), then
   |f_i(z)| <= degree(f_i) B |z|^degree(f_i)

 so |R(z)| has an explicit upper bound, growing no worse than |z|^degree(f)

 so 

 |a_n| <= 1/2/pi (2 pi r) deg(f) B (|t|+r)^deg(f)/r^(n+1)
        = deg(f) B (|t|+r)^deg(f)/ r^n

 This is valid for all r where R is analytic.

 Let d = disc(f), then R should be analytic if r < min(|t-s| for d(s) = 0)
 
 if this is > 1, the a_n -> 0    

=#

