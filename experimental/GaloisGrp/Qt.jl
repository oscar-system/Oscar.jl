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

=======================================================

OK, s is a root, so the coeffs s = sum s_n (x-t)^n 
satisfy, by above, |s_n| <= B/r^n for some explicit r > 1

Lemma:
  let s, t be power series with
    |s_n| <= B (n+1)^k/r^n for some B, k
    |t_n| <= C (n+1)^l/r^n
  then
    st = sum d_i (x-t)^n with
    |d_n| <= BC (n+1)^(k+l+1)/r^n

pf:
  Cauchy product:

  |d_n| =  |sum_{i+j = n} s_i t_j|
        <= sum |s_i t_j|
        <= sum Bi^k/r^i C j^l/r^j
        <= BC sum (n+1)^k (n+1)^l/r^(i+j) = BC (n+1)^(k+l) (n+1) 1/r^n

Clear: 
  |(s+t)_n| <= ...


Lemma:
  if |s_n| <= B (n+1)^k/r^n, then the largest coeff is either at
    floor(x) or ceil(x) for x = (k/log r) - 1

pf:
 as a function of n:
   (B (n+1)^k/r^n)' = B k (n+1)^(k-1)/r^n - B (n+1)^k log r /r^n
 which is zero iff
   k - (n+1) log r = 0 or n+1 = k/log r


Possibly we can get better estimates for sum_{i+j=n} i^k j^l
on particular for powering:
  sum_{|alpha| = k} prod alpha^k

===============================================================
Also see 
https://math.uni-paderborn.de/fileadmin/mathematik/AG-Computeralgebra/Publications-klueners/function.pdf

for analysis of the denominator and the infinite valuations  
=#

