using Oscar
R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
_F = free_module(R, 3) 


kQ = monoid_algebra([[1,0],[0,1]],QQ)
x,y = gens(kQ)
M1 = quotient_ring_as_module(ideal(kQ,[x^3,y^3]))
M2 = quotient_ring_as_module(ideal(kQ,[x^3,y^3]))
_M = direct_sum(M1,M2)[1]
M,_ = sub(_M,[x*_M[1],y*_M[2]])
N,_ = quo(M,sub(M,[y^2*M[1]-x*y*M[2],x*y*M[1]-x^2*M[2]])[1])

is_graded(N)
gens(N)
relations(N)

#injective resolution of N 
I = injective_resolution(N,3)
I.cochain_maps[1]

irr_res = I.Q_graded_part
i = irr_res.cochain_maps[1]
W0 = irr_res.irr_sums[1].Q_graded_part
_i = hom(N,W0,matrix(i))

#free resolution of N 
F = free_resolution(N)
p = map(F,0)



# computation with 4 generators
kQ = monoid_algebra([[1,0],[0,1]],QQ)
x,y = gens(kQ)
I1 = quotient_ring_as_module(ideal(kQ,[x^3]))
I2 = quotient_ring_as_module(ideal(kQ,[y^3]))
I3 = quotient_ring_as_module(ideal(kQ,[]))
I4 = quotient_ring_as_module(ideal(kQ,[]))

_M = direct_sum([I1,I2,I3,I4])[1]
M,_ = sub(_M,[x*_M[1],y*_M[2],x*y*_M[3],x*y*_M[4]])
N,_ = quo(M,sub(M,[y*M[1] - M[3],x*M[2]-M[4]])[1])

gens(N)
relations(N)

res = injective_resolution(N,3)
irr_res = res.Q_graded_part

i = res.cochain_maps[1]

I = irr_res.cochain_complex
i = map(I,0)
# why is this map not from N to I^0???

f = free_resolution(N)
C = f.C

e = matrix(map(C,0))