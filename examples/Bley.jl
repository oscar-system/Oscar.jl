#=

Qx, x = QQ["x"]
k, _ = number_field([x^2-65, x^2-185])
k, _ = absolute_simple_field(k)
Ik = idele_class_gmodule(k)

C = cyclotomic_field(ClassField, 65)
C = subfields(C, degree = 4)
c = number_field(AbsSimpleNumField, C[1])

K = compositum(k, c)

IK = idele_class_gmodule(K[1])

Ic = idele_class_gmodule(c)
Ic = add_prime(Ic, [37])
Ic = change_precision(Ic, [(37, 1)])

Psi = induce_hom(Ic, IK, K[3])
Phi = induce_hom(Ik, IK, K[2])



Zg = free_res(group_algebra(ZZ, domain(Ik.mG)); side = :right)
ZG = free_res(group_algebra(ZZ, domain(IK.mG)); side = :right)
Zc = free_res(group_algebra(ZZ, domain(Ic.mG)); side = :right)
qKk = Oscar.GaloisCohomology_Mod.fixed_group(IK.mG, Ik.mG, K[2])
qKc = Oscar.GaloisCohomology_Mod.fixed_group(IK.mG, Ic.mG, K[3])

inf_cK = change_group(ZG, Zc, qKc[2])

hc = hom(Zc, Ic.data[1])
hK = hom(ZG, IK.data[1])

kc = kernel(map(hc, 2))
qc = quo(kc[1], image(map(hc, 1))[1])
sc = snf(qc[1])
fc = kc[2](preimage(qc[2], sc[2](sc[1][1])))
Oscar.GaloisCohomology_Mod.magic_inf(Ic.data[1], fc, inf_cK[2], Psi, hK[2])
fcK = ans;
kK = kernel(map(hK, 2))
qK = quo(kK[1], image(map(hK, 1))[1])
sK = snf(qK[1])
preimage(kK[2], fcK)
map(hK, 2)(fcK)

