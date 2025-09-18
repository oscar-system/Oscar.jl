using Oscar

function cyclic5(k::Field=QQ)
    R, (a,b,c,d,x) = k[:a,:b,:c,:d,:x]
    F = [
        a + b + c + d + x,
        a*b + b*c + c*d + d*x + x*a,
        a*b*c + b*c*d + c*d*x + d*x*a + x*a*b,
        a*b*c*d + b*c*d*x + c*d*x*a + d*x*a*b + x*a*b*c,
        a*b*c*d*x - 1
    ]

    return ideal(F), lex(R), default_ordering(R)
end

function cyclic6(k::Field=QQ)
    R, (z0,z1,z2,z3,z4,z5) = k[:z0,:z1,:z2,:z3,:z4,:z5]
    F = [
      z0 + z1 + z2 + z3 + z4 + z5,
      z0*z1 + z1*z2 + z2*z3 + z3*z4 + z4*z5 + z5*z0,
      z0*z1*z2 + z1*z2*z3 + z2*z3*z4 + z3*z4*z5 + z4*z5*z0 + z5*z0*z1,
      z0*z1*z2*z3 + z1*z2*z3*z4 + z2*z3*z4*z5 + z3*z4*z5*z0 + z4*z5*z0*z1 
      + z5*z0*z1*z2,
      z0*z1*z2*z3*z4 + z1*z2*z3*z4*z5 + z2*z3*z4*z5*z0 + z3*z4*z5*z0*z1 
      + z4*z5*z0*z1*z2 + z5*z0*z1*z2*z3,
      z0*z1*z2*z3*z4*z5 - 1
    ]

    return ideal(F), lex(R), default_ordering(R)
end

function cyclic7(k::Field=QQ)
    R, (z0, z1, z2, z3, z4, z5, z6) = QQ[:z0, :z1, :z2, :z3, :z4, :z5, :z6]
    F = [
        z0 + z1 + z2 + z3 + z4 + z5 + z6,
        z0 * z1 + z1 * z2 + z2 * z3 + z3 * z4 + z4 * z5 + z5 * z6 + z6 * z0,
        z0 * z1 * z2 + z1 * z2 * z3 + z2 * z3 * z4 + z3 * z4 * z5 + z4 * z5 * z6 + z5 * z6 * z0 + z6 * z0 * z1,
        z0 * z1 * z2 * z3 + z1 * z2 * z3 * z4 + z2 * z3 * z4 * z5 + z3 * z4 * z5 * z6 + z4 * z5 * z6 * z0
        + z5 * z6 * z0 * z1 + z6 * z0 * z1 * z2,
        z0 * z1 * z2 * z3 * z4 + z1 * z2 * z3 * z4 * z5 + z2 * z3 * z4 * z5 * z6 + z3 * z4 * z5 * z6 * z0
        + z4 * z5 * z6 * z0 * z1 + z5 * z6 * z0 * z1 * z2 + z6 * z0 * z1 * z2 * z3,
        z0 * z1 * z2 * z3 * z4 * z5 + z1 * z2 * z3 * z4 * z5 * z6 + z2 * z3 * z4 * z5 * z6 * z0 + z3 * z4 * z5 * z6 * z0 * z1
        + z4 * z5 * z6 * z0 * z1 * z2 + z5 * z6 * z0 * z1 * z2 * z3 + z6 * z0 * z1 * z2 * z3 * z4,
        z0 * z1 * z2 * z3 * z4 * z5 * z6 - 1
    ]

    return ideal(F), lex(R), default_ordering(R)
end

function cyclic8(k::Field=QQ)
    R, (z0, z1, z2, z3, z4, z5, z6, z7) = k[:z0, :z1, :z2, :z3, :z4, :z5, :z6, :z7]

    F = [
        z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7,

        z0*z1 + z1*z2 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z7 + z7*z0,

        z0*z1*z2 + z1*z2*z3 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z7
        + z6*z7*z0 + z7*z0*z1,

        z0*z1*z2*z3 + z1*z2*z3*z4 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z7
        + z5*z6*z7*z0 + z6*z7*z0*z1 + z7*z0*z1*z2,

        z0*z1*z2*z3*z4 + z1*z2*z3*z4*z5 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z7
        + z4*z5*z6*z7*z0 + z5*z6*z7*z0*z1 + z6*z7*z0*z1*z2 + z7*z0*z1*z2*z3,

        z0*z1*z2*z3*z4*z5 + z1*z2*z3*z4*z5*z6 + z2*z3*z4*z5*z6*z7 + z3*z4*z5*z6*z7*z0
        + z4*z5*z6*z7*z0*z1 + z5*z6*z7*z0*z1*z2 + z6*z7*z0*z1*z2*z3 + z7*z0*z1*z2*z3*z4,
    
        z0*z1*z2*z3*z4*z5*z6 + z1*z2*z3*z4*z5*z6*z7 + z2*z3*z4*z5*z6*z7*z0
        + z3*z4*z5*z6*z7*z0*z1 + z4*z5*z6*z7*z0*z1*z2 + z5*z6*z7*z0*z1*z2*z3
        + z6*z7*z0*z1*z2*z3*z4 + z7*z0*z1*z2*z3*z4*z5,

        z0*z1*z2*z3*z4*z5*z6*z7 - 1
    ]

    return ideal(F), lex(R), default_ordering(R)
end

