# These are three examples that produce the correct result.
#
# To be put in a proper test suite.
R, (x,y,z,w) = QQ["x", "y", "z", "w"]
I = ideal(R, [z^3, y^3*z*w^2, x^2*y^4*w^2, x*y*z^2])
@show Oscar._hilbert_numerator_from_leading_ideal(I)

W = [1 2 3 4; 0 0 5 8]
P, (X, Y, Z, W) = grade(R, W)

I = ideal(P, [X^2, Y, Z^3])

@show Oscar._hilbert_numerator_from_leading_ideal(I)

R, (x,y,z) = QQ["x", "y", "z"]
W = [1 1 1; 0 0 -1]
Q, (X, Y, Z) = grade(R, W)
I = ideal(Q, [X^3*Y, Y*Z^2, y^2*Z, Z^4])
Oscar._hilbert_numerator_from_leading_ideal(I)
