P3 = projective_space(QQ,[:x0,:x1,:x2,:x3]);
(x0,x1,x2,x3) = homogeneous_coordinates(P3);
R=homogeneous_coordinate_ring(P3);
########################################
I = ideal(R,[x0,x3]);
#L,inc_L = sub(P3,I);
I1 = ideal(R,[x0-x1,x2]);
#L1,inc_L1 = sub(P3, I1);
I2 = ideal(R,[x1,x2]);
#L2,inc_L2 = sub(P3,ideal(R,[x1,x2]));
I3 = ideal(R,[x1,x2-x3]);
#L3,inc_L3 = sub(P3,ideal(R,[x1,x2-x3]));
# (X2,E2,inc_E2) = blowup(inc_L2); This does not work.
IH = ideal(x0-x2);
##########################################
Y = covered_scheme(P3);
II1=ideal_sheaf(P3,I1);
II2=ideal_sheaf(P3,I2);
II3=ideal_sheaf(P3,I3);
II=ideal_sheaf(P3,I);
IIH=ideal_sheaf(P3,IH);
H=effective_cartier_divisor(IIH);
###############################################
#We do the first blowup, by blowing up L
bl = blow_up(Y,II);
X= domain(bl);
II11=strict_transform(bl,II1);
II12=strict_transform(bl,II2);
II13=strict_transform(bl,II3);
H1=strict_transform(bl,H);
E=exceptional_divisor(bl);
##################################################
#We do the second blowup, by blowing up L2
bl2 = blow_up(X,II12);
X2= domain(bl2);
II21=strict_transform(bl2,II11);
II23=strict_transform(bl2,II13);
H2=strict_transform(bl2,H1);
EE=strict_transform(bl2,E);
E2=exceptional_divisor(bl2);
###################################################
#We do the third blowup, by blowing up L1
bl1 = blow_up(X2,II21);
X1= domain(bl1);
II33=strict_transform(bl1,II23);
H3=strict_transform(bl1,H2);
EEE=strict_transform(bl1,EE);
EE2=strict_transform(bl1,E2);
E1=exceptional_divisor(bl1);
####################################################
#We do the fourth blowup, by blowing up L3
bl3 = blow_up(X1,II33);
X_final= domain(bl3);
H_final = strict_transform(bl3,H3);
E_final=strict_transform(bl3,EEE);
E2_final=strict_transform(bl3,EE2);
E1_final=strict_transform(bl3,E1);
E3_final=exceptional_divisor(bl3);
###################################################
Int1=intersect(E_final,E1_final);
Int23=intersect(E2_final,E3_final);
Int21=intersect(E2_final,E1_final);
E2E3E1=intersect(Int23,E1_final);
integral(E2E3E1);
##################################################
#Now we want to do intersections inside the surface. Lets start with E2.
inc_E2 = embedding(E2_final);
l21=pullback(inc_E2,E1_final);
l23=pullback(inc_E2,E3_final);
Int123=intersect(l21,l23);
integral(Int123);
#####################################################
inc_E1 = embedding(E1_final);
l21=pullback(inc_E2,E1_final);
l23=pullback(inc_E2,E3_final);
Int123=intersect(l21,l23);
integral(Int123);
kX=-4*H_final+E1_final+E2_final+E3_final+E_final;
###################################################
inc_H = embedding(H_final);
l21=pullback(inc_E2,E1_final);
l23=pullback(inc_E2,E3_final);
Int123=intersect(l21,l23);
integral(Int123);

E2_sub, inc_E2 = sub(ideal_sheaf(E2_final))
pb_E1 = pullback(inc_E2, E1_final)
pb_E1_weil = weil_divisor(pb_E1)
int1 = integral(intersect(pb_E1, pb_E1_weil))
move = Oscar.move_divisor(pb_E1_weil)
int2 = integral(intersect(pb_E1, move))
@test int1 == int2




