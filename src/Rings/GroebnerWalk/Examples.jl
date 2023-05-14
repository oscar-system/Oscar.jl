function katsura4()
    dim = 5
    ve = [1, 1, 1, 1, 1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x0, x1, x2, x3, x4) = polynomial_ring(
        GF(32003),
        ["x0", "x1", "x2", "x3", "x4"],
        ordering= :degrevlex,
    )

    f1 = 2 * x4^2 + 2 * x3^2 + 2 * x2^2 + 2 * x1^2 + x0^2 - x0
    f2 = 2 * x3 * x4 + 2 * x3 * x2 + 2 * x2 * x1 + 2 * x1 * x0 - x1
    f3 = 2 * x4 * x2 + 2 * x3 * x1 + 2 * x1^2 + 2 * x2 * x0 - x2
    f4 = 2 * x4 * x1 + 2 * x2 * x1 + 2 * x3 * x0 - x3
    f5 = 2 * x4 + 2 * x3 + 2 * x2 + 2 * x1 + x0 - 1
    return ideal(R, [f1, f2, f3, f4, f5])
end


#Katsura6
function katsura5()
    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x0, x1, x2, x3, x4, x5) = polynomial_ring(
        QQ,
        ["x0", "x1", "x2", "x3", "x4", "x5"],
        ordering= :degrevlex,
    )
    f1 = 2 * x5^2 + 2 * x4^2 + 2 * x3^2 + 2 * x2^2 + 2 * x1^2 + x0^2 - x0
    f2 = 2 * x5 * x4 + 2 * x4 * x3 + 2 * x3 * x2 + 2 * x2 * x1 + 2 * x1 * x0 - x1
    f3 = 2 * x5 * x3 + 2 * x4 * x2 + 2 * x3 * x1 + 2 * x1^2 + 2 * x2 * x0 - x2
    f4 = 2 * x5 * x2 + 2 * x4 * x1 + 2 * x2 * x1 + 2 * x3 * x0 - x3
    f5 = x2^2 + 2 * x5 * x0 + 2 * x4 * x0 + 2 * x3 * x0 - x4
    f6 = 2 * x5 + 2 * x4 + 2 * x3 + 2 * x2 + 2 * x1 + x0 - 1

    return ideal(R, [f1, f2, f3, f4, f5, f6])
end


function katsura6()
    dim = 7
    ve = [1, 1, 1, 1, 1, 1,1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x0,x1,x2,x3,x4,x5,x6) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x0","x1","x2","x3","x4","x5","x6"],
        ordering = Singular.ordering_M(StartOrd),
    )
    f1 = x0+2*x1+2*x2+2*x3+2*x4+2*x5+2*x6-1
    f2 = 2*x2*x3+2*x1*x4+2*x0*x5+2*x1*x6-x5
    f3 = x2^2+2*x1*x3+2*x0*x4+2*x1*x5+2*x2*x6-x4
    f4 = 2*x1*x2+2*x0*x3+2*x1*x4+2*x2*x5+2*x3*x6-x3
    f5 = x1^2+2*x0*x2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6-x2
    f6 = 2*x0*x1+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6-x1
    f7 = x0^2+2*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2-x0

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7])
end



function katsura7()
    dim = 8
    ve = [1, 1, 1, 1, 1, 1, 1,1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (u0, u1, u2, u3, u4, u5, u6, u7) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["u0","u1","u2","u3", "u4","u5","u6", "u7"],
        ordering = Singular.ordering_M(StartOrd),
    )
f1=    u0+2*u1+2*u2+2*u3+2*u4+2*u5+2*u6+2*u7-1
f2=    u3^2+2*u2*u4+2*u1*u5+2*u0*u6+2*u1*u7-u6
f3=    2*u2*u3+2*u1*u4+2*u0*u5+2*u1*u6+2*u2*u7-u5
f4=    u2^2+2*u1*u3+2*u0*u4+2*u1*u5+2*u2*u6+2*u3*u7-u4
f5=    2*u1*u2+2*u0*u3+2*u1*u4+2*u2*u5+2*u3*u6+2*u4*u7-u3
f6=    u1^2+2*u0*u2+2*u1*u3+2*u2*u4+2*u3*u5+2*u4*u6+2*u5*u7-u2
f7=    2*u0*u1+2*u1*u2+2*u2*u3+2*u3*u4+2*u4*u5+2*u5*u6+2*u6*u7-u1
f8=    u0^2+2*u1^2+2*u2^2+2*u3^2+2*u4^2+2*u5^2+2*u6^2+2*u7^2-u0

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7, f8])
end

function katsura8()
    dim = 9
    ve = [1, 1, 1, 1, 1, 1,1,1,1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (u0,u1,u2,u3,u4,u5,u6,u7,u8) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["u0","u1","u2","u3","u4","u5","u6","u7","u8"],
        ordering = Singular.ordering_M(StartOrd),
    )
    f1 = u0+2*u1+2*u2+2*u3+2*u4+2*u5+2*u6+2*u7+2*u8-1
    f2 = 2*u3*u4+2*u2*u5+2*u1*u6+2*u0*u7+2*u1*u8-u7
    f3 = u3^2+2*u2*u4+2*u1*u5+2*u0*u6+2*u1*u7+2*u2*u8-u6
    f4 = 2*u2*u3+2*u1*u4+2*u0*u5+2*u1*u6+2*u2*u7+2*u3*u8-u5
    f5 = u2^2+2*u1*u3+2*u0*u4+2*u1*u5+2*u2*u6+2*u3*u7+2*u4*u8-u4
    f6 = 2*u1*u2+2*u0*u3+2*u1*u4+2*u2*u5+2*u3*u6+2*u4*u7+2*u5*u8-u3
    f7 = u1^2+2*u0*u2+2*u1*u3+2*u2*u4+2*u3*u5+2*u4*u6+2*u5*u7+2*u6*u8-u2
    f8 = 2*u0*u1+2*u1*u2+2*u2*u3+2*u3*u4+2*u4*u5+2*u5*u6+2*u6*u7+2*u7*u8-u1
    f9 = u0^2+2*u1^2+2*u2^2+2*u3^2+2*u4^2+2*u5^2+2*u6^2+2*u7^2+2*u8^2-u0

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7, f8, f9])
end


function cyclic7()
    dim = 7
    ve = [1, 1, 1, 1, 1, 1, 1]
    example = "Cyclic7"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5, x6, x7) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7"],
        ordering = Singular.ordering_M(StartOrd),
    )
    f1 = x1 + x2 + x3 + x4 + x5 + x6 + x7
    f2 = x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x5 * x6 + x1 * x7 + x6 * x7
    f3 =
        x1 * x2 * x3 +
        x2 * x3 * x4 +
        x3 * x4 * x5 +
        x4 * x5 * x6 +
        x1 * x2 * x7 +
        x1 * x6 * x7 +
        x5 * x6 * x7
    f4 =
        x1 * x2 * x3 * x4 +
        x2 * x3 * x4 * x5 +
        x3 * x4 * x5 * x6 +
        x1 * x2 * x3 * x7 +
        x1 * x2 * x6 * x7 +
        x1 * x5 * x6 * x7 +
        x4 * x5 * x6 * x7
    f5 =
        x1 * x2 * x3 * x4 * x5 +
        x2 * x3 * x4 * x5 * x6 +
        x1 * x2 * x3 * x4 * x7 +
        x1 * x2 * x3 * x6 * x7 +
        x1 * x2 * x5 * x6 * x7 +
        x1 * x4 * x5 * x6 * x7 +
        x3 * x4 * x5 * x6 * x7
    f6 =
        x1 * x2 * x3 * x4 * x5 * x6 +
        x1 * x2 * x3 * x4 * x5 * x7 +
        x1 * x2 * x3 * x4 * x6 * x7 +
        x1 * x2 * x3 * x5 * x6 * x7 +
        x1 * x2 * x4 * x5 * x6 * x7 +
        x1 * x3 * x4 * x5 * x6 * x7 +
        x2 * x3 * x4 * x5 * x6 * x7
    f7 = x1 * x2 * x3 * x4 * x5 * x6 * x7 - 1

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7])
end


function cyclic6()
    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "Cyclic6"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = x1 + x2 + x3 + x4 + x5 + x6
    f2 = x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x1 * x6 + x5 * x6
    f3 =
        x1 * x2 * x3 +
        x2 * x3 * x4 +
        x3 * x4 * x5 +
        x1 * x2 * x6 +
        x1 * x5 * x6 +
        x4 * x5 * x6
    f4 =
        x1 * x2 * x3 * x4 +
        x2 * x3 * x4 * x5 +
        x1 * x2 * x3 * x6 +
        x1 * x2 * x5 * x6 +
        x1 * x4 * x5 * x6 +
        x3 * x4 * x5 * x6
    f5 =
        x1 * x2 * x3 * x4 * x5 +
        x1 * x2 * x3 * x4 * x6 +
        x1 * x2 * x3 * x5 * x6 +
        x1 * x2 * x4 * x5 * x6 +
        x1 * x3 * x4 * x5 * x6 +
        x2 * x3 * x4 * x5 * x6
    f6 = x1 * x2 * x3 * x4 * x5 * x6 - 1

    return ideal(R, [f1, f2, f3, f4, f5, f6])
end

function cyclic5()
    dim = 5
    ve = [1, 1, 1, 1, 1]
    example = "Cyclic5"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (v, w, x, y, z) = polynomial_ring(
        Singular.n_ZpField(32003),
        ["v", "w", "x", "y", "z"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = v + w + x + y + z
    f2 = v * w + w * x + x * y + v * z + y * z
    f3 = v * w * x + w * x * y + v * w * z + v * y * z + x * y * z
    f4 =
        v * w * x * y +
        v * w * x * z +
        v * w * y * z +
        v * x * y * z +
        w * x * y * z
    f5 = v * w * x * y * z - 1

    return ideal(R, [f1, f2, f3, f4, f5])
end

function cyclic4()
    dim = 4
    ve = [1, 1, 1, 1]
    example = "Cyclic4"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (w, x, y, z) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["w", "x", "y", "z"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = w + x + y + z
    f2 = w * x + x * y + w * z + y * z
    f3 = w * x * y + w * x * z + w * y * z + x * y * z
    f4 = w * x * y * z - 1
    return ideal(R, [f1, f2, f3, f4])
end

function cyclic8()
    dim = 8
    ve = [1, 1, 1, 1, 1 ,1 ,1 ,1]
    example = "Cyclic8"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1,x2,x3,x4,x5,x6,x7,x8) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1","x2", "x3","x4","x5","x6","x7","x8"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = x1+x2+x3+x4+x5+x6+x7+x8
    f2 = x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7+x1*x8+x7*x8
    f3 = x1*x2*x3+x2*x3*x4+x3*x4*x5+x4*x5*x6+x5*x6*x7+x1*x2*x8+x1*x7*x8+x6*x7*x8
    f4 = x1*x2*x3*x4+x2*x3*x4*x5+x3*x4*x5*x6+x4*x5*x6*x7+x1*x2*x3*x8+x1*x2*x7*x8+x1*x6*x7*x8+x5*x6*x7*x8
    f5 = x1*x2*x3*x4*x5+x2*x3*x4*x5*x6+x3*x4*x5*x6*x7+x1*x2*x3*x4*x8+x1*x2*x3*x7*x8+x1*x2*x6*x7*x8+x1*x5*x6*x7*x8+x4*x5*x6*x7*x8
    f6 = x1*x2*x3*x4*x5*x6+x2*x3*x4*x5*x6*x7+x1*x2*x3*x4*x5*x8+x1*x2*x3*x4*x7*x8+x1*x2*x3*x6*x7*x8+x1*x2*x5*x6*x7*x8+x1*x4*x5*x6*x7*x8+x3*x4*x5*x6*x7*x8
    f7 = x1*x2*x3*x4*x5*x6*x7+x1*x2*x3*x4*x5*x6*x8+x1*x2*x3*x4*x5*x7*x8+x1*x2*x3*x4*x6*x7*x8+x1*x2*x3*x5*x6*x7*x8+x1*x2*x4*x5*x6*x7*x8+x1*x3*x4*x5*x6*x7*x8+x2*x3*x4*x5*x6*x7*x8
    f8 = x1*x2*x3*x4*x5*x6*x7*x8-1
    return ideal(R, [f1, f2, f3, f4, f5, f6, f7, f8])
end


function eco6()
    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "eco6"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )


    f1 = x1 + x2 + x3 + x4 + x5 + 1
    f2 = x5 * x6 - 5
    f3 = x1 * x5 * x6 + x4 * x6 - 4
    f4 = x1 * x4 * x6 + x2 * x5 * x6 + x3 * x6 - 3
    f5 = x1 * x3 * x6 + x2 * x4 * x6 + x3 * x5 * x6 + x2 * x6 - 2
    f6 = x1 * x2 * x6 + x2 * x3 * x6 + x3 * x4 * x6 + x4 * x5 * x6 + x1 * x6 - 1

    return ideal(R, [f1, f2, f3, f4, f5, f6])
end

function eco7()
    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "eco7"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = x1 + x2 + x3 + x4 + x5 + x6 + 1
    f2 = x6 * x7 - 6
    f3 = x1 * x6 * x7 + x5 * x7 - 5
    f4 = x1 * x5 * x7 + x2 * x6 * x7 + x4 * x7 - 4
    f5 = x1 * x4 * x7 + x2 * x5 * x7 + x3 * x6 * x7 + x3 * x7 - 3
    f6 = x1 * x3 * x7 + x2 * x4 * x7 + x3 * x5 * x7 + x4 * x6 * x7 + x2 * x7 - 2
    f7 =
        x1 * x2 * x7 +
        x2 * x3 * x7 +
        x3 * x4 * x7 +
        x4 * x5 * x7 +
        x5 * x6 * x7 +
        x1 * x7 - 1


    return ideal(R, [f1, f2, f3, f4, f5, f6, f7])
end

function noon5()
    dim = 5
    ve = [1, 1, 1, 1, 1]
    example = "noon5"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 =
        10 * x1^2 * x5 + 10 * x2^2 * x5 + 10 * x3^2 * x5 + 10 * x4^2 * x5 -
        11 * x5 + 10
    f2 =
        10 * x1^2 * x4 + 10 * x2^2 * x4 + 10 * x3^2 * x4 + 10 * x4 * x5^2 -
        11 * x4 + 10
    f3 =
        10 * x1^2 * x3 + 10 * x2^2 * x3 + 10 * x3 * x4^2 + 10 * x3 * x5^2 -
        11 * x3 + 10
    f4 =
        10 * x1 * x2^2 + 10 * x1 * x3^2 + 10 * x1 * x4^2 + 10 * x1 * x5^2 -
        11 * x1 + 10
    f5 =
        10 * x1^2 * x2 + 10 * x2 * x3^2 + 10 * x2 * x4^2 + 10 * x2 * x5^2 -
        11 * x2 + 10

    return ideal(R, [f1, f2, f3, f4, f5])
end

function noon6()
    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "noon6"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5, x6) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 =
        10 * x1^2 * x6 +
        10 * x2^2 * x6 +
        10 * x3^2 * x6 +
        10 * x4^2 * x6 +
        10 * x5^2 * x6 - 11 * x6 + 10
    f2 =
        10 * x1^2 * x5 +
        10 * x2^2 * x5 +
        10 * x3^2 * x5 +
        10 * x4^2 * x5 +
        10 * x5 * x6^2 - 11 * x5 + 10
    f3 =
        10 * x1^2 * x4 +
        10 * x2^2 * x4 +
        10 * x3^2 * x4 +
        10 * x4 * x5^2 +
        10 * x4 * x6^2 - 11 * x4 + 10
    f4 =
        10 * x1^2 * x3 +
        10 * x2^2 * x3 +
        10 * x3 * x4^2 +
        10 * x3 * x5^2 +
        10 * x3 * x6^2 - 11 * x3 + 10
    f5 =
        10 * x1 * x2^2 +
        10 * x1 * x3^2 +
        10 * x1 * x4^2 +
        10 * x1 * x5^2 +
        10 * x1 * x6^2 - 11 * x1 + 10
    f6 =
        10 * x1^2 * x2 +
        10 * x2 * x3^2 +
        10 * x2 * x4^2 +
        10 * x2 * x5^2 +
        10 * x2 * x6^2 - 11 * x2 + 10

    return ideal(R, [f1, f2, f3, f4, f5, f6])
end

function noon7()
    dim = 7
    ve = [1, 1, 1, 1, 1, 1, 1]
    example = "noon7"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5, x6, x7) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 =
        10 * x1^2 * x7 +
        10 * x2^2 * x7 +
        10 * x3^2 * x7 +
        10 * x4^2 * x7 +
        10 * x5^2 * x7 +
        10 * x6^2 * x7 - 11 * x7 + 10
    f2 =
        10 * x1^2 * x6 +
        10 * x2^2 * x6 +
        10 * x3^2 * x6 +
        10 * x4^2 * x6 +
        10 * x5^2 * x6 +
        10 * x6 * x7^2 - 11 * x6 + 10
    f3 =
        10 * x1^2 * x5 +
        10 * x2^2 * x5 +
        10 * x3^2 * x5 +
        10 * x4^2 * x5 +
        10 * x5 * x6^2 +
        10 * x5 * x7^2 - 11 * x5 + 10
    f4 =
        10 * x1^2 * x4 +
        10 * x2^2 * x4 +
        10 * x3^2 * x4 +
        10 * x4 * x5^2 +
        10 * x4 * x6^2 +
        10 * x4 * x7^2 - 11 * x4 + 10
    f5 =
        10 * x1^2 * x3 +
        10 * x2^2 * x3 +
        10 * x3 * x4^2 +
        10 * x3 * x5^2 +
        10 * x3 * x6^2 +
        10 * x3 * x7^2 - 11 * x3 + 10
    f6 =
        10 * x1 * x2^2 +
        10 * x1 * x3^2 +
        10 * x1 * x4^2 +
        10 * x1 * x5^2 +
        10 * x1 * x6^2 +
        10 * x1 * x7^2 - 11 * x1 + 10
    f7 =
        10 * x1^2 * x2 +
        10 * x2 * x3^2 +
        10 * x2 * x4^2 +
        10 * x2 * x5^2 +
        10 * x2 * x6^2 +
        10 * x2 * x7^2 - 11 * x2 + 10

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7])
end

function noon8()
    dim = 8
    ve = [1, 1, 1, 1, 1, 1, 1, 1]
    example = "noon8"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, x5, x6, x7, x8) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 =
        10 * x1^2 * x8 +
        10 * x2^2 * x8 +
        10 * x3^2 * x8 +
        10 * x4^2 * x8 +
        10 * x5^2 * x8 +
        10 * x6^2 * x8 +
        10 * x7^2 * x8 - 11 * x8 + 10
    f2 =
        10 * x1^2 * x7 +
        10 * x2^2 * x7 +
        10 * x3^2 * x7 +
        10 * x4^2 * x7 +
        10 * x5^2 * x7 +
        10 * x6^2 * x7 +
        10 * x7 * x8^2 - 11 * x7 + 10
    f3 =
        10 * x1^2 * x6 +
        10 * x2^2 * x6 +
        10 * x3^2 * x6 +
        10 * x4^2 * x6 +
        10 * x5^2 * x6 +
        10 * x6 * x7^2 +
        10 * x6 * x8^2 - 11 * x6 + 10
    f4 =
        10 * x1^2 * x5 +
        10 * x2^2 * x5 +
        10 * x3^2 * x5 +
        10 * x4^2 * x5 +
        10 * x5 * x6^2 +
        10 * x5 * x7^2 +
        10 * x5 * x8^2 - 11 * x5 + 10
    f5 =
        10 * x1^2 * x4 +
        10 * x2^2 * x4 +
        10 * x3^2 * x4 +
        10 * x4 * x5^2 +
        10 * x4 * x6^2 +
        10 * x4 * x7^2 +
        10 * x4 * x8^2 - 11 * x4 + 10
    f6 =
        10 * x1^2 * x3 +
        10 * x2^2 * x3 +
        10 * x3 * x4^2 +
        10 * x3 * x5^2 +
        10 * x3 * x6^2 +
        10 * x3 * x7^2 +
        10 * x3 * x8^2 - 11 * x3 + 10
    f7 =
        10 * x1 * x2^2 +
        10 * x1 * x3^2 +
        10 * x1 * x4^2 +
        10 * x1 * x5^2 +
        10 * x1 * x6^2 +
        10 * x1 * x7^2 +
        10 * x1 * x8^2 - 11 * x1 + 10
    f8 =
        10 * x1^2 * x2 +
        10 * x2 * x3^2 +
        10 * x2 * x4^2 +
        10 * x2 * x5^2 +
        10 * x2 * x6^2 +
        10 * x2 * x7^2 +
        10 * x2 * x8^2 - 11 * x2 + 10

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7, f8])
end

function redeco7()
    dim = 7
    ve = [1, 1, 1, 1, 1, 1, 1]
    example = "redeco7"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, u7, x5, x6) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "u7", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = -6 * u7 + x6
    f2 = x1 + x2 + x3 + x4 + x5 + x6 + 1
    f3 = x1 * x6 - 5 * u7 + x5
    f4 = x1 * x5 + x2 * x6 + x4 - 4 * u7
    f5 = x1 * x4 + x2 * x5 + x3 * x6 + x3 - 3 * u7
    f6 = x1 * x3 + x2 * x4 + x3 * x5 + x4 * x6 + x2 - 2 * u7
    f7 = x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x5 * x6 + x1 - u7

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7])
end

function redeco8()
    dim = 8
    ve = [1, 1, 1, 1, 1, 1, 1, 1]
    example = "redeco8"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3, x4, u8, x5, x6, x7) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3", "x4", "u8", "x5", "x6", "x7"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = -7 * u8 + x7
    f2 = x1 + x2 + x3 + x4 + x5 + x6 + x7 + 1
    f3 = x1 * x7 - 6 * u8 + x6
    f4 = x1 * x6 + x2 * x7 + x5 - 5 * u8
    f5 = x1 * x5 + x2 * x6 + x3 * x7 + x4 - 4 * u8
    f6 = x1 * x4 + x2 * x5 + x3 * x6 + x4 * x7 + x3 - 3 * u8
    f7 = x1 * x3 + x2 * x4 + x3 * x5 + x4 * x6 + x5 * x7 + x2 - 2 * u8
    f8 = x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x5 * x6 + x6 * x7 + x1 - u8

    return ideal(R, [f1, f2, f3, f4, f5, f6, f7, f8])
end

function wang91()
    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "Wang-91"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x3, x2, x1, x0, b, a) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x3", "x2", "x1", "x0", "b", "a"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = 3 * x2 * x1 * a + 3 * x0^2
    f2 = 3 * x2 * x1 * b + 3 * x3^2
    f3 = 3 * x3 * x1 * b + 3 * x1 * x0 * a + 3 * x2^2
    f4 = 3 * x3 * x2 * b + 3 * x2 * x0 * a + 3 * x1^2

    return ideal(R, [f1, f2, f3, f4])
end

function cohn4()
    dim = 4
    ve = [1, 1, 1, 1]
    example = "cohn4"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x, y, z, t) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x", "y", "z", "t"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 =
        -x^3 * y^2 + 2 * x^2 * y^2 * z - x^2 * y * z^2 - 144 * x^2 * y^2 -
        207 * x^2 * y * z +
        288 * x * y^2 * z +
        78 * x * y * z^2 +
        x * z^3 - 3456 * x^2 * y - 5184 * x * y^2 - 9504 * x * y * z -
        432 * x * z^2 - 248832 * x * y + 62208 * x * z - 2985984 * x
    f2 =
        y^3 * t^3 - y^2 * z * t^3 + 4 * y^3 * t^2 - 2 * y^2 * z * t^2 +
        72 * y^2 * t^3 +
        71 * y * z * t^3 +
        288 * y^2 * t^2 +
        360 * y * z * t^2 +
        6 * z^2 * t^2 +
        1728 * y * t^3 - 464 * z * t^3 +
        432 * y * z * t +
        8 * z^2 * t +
        6912 * y * t^2 - 4320 * z * t^2 +
        13824 * t^3 +
        z^2 - 13824 * z * t + 55296 * t^2 - 13824 * z
    f3 =
        x^2 * y * t^3 - 2 * x * y^2 * t^3 + y^3 * t^3 + 8 * x^2 * y * t^2 -
        12 * x * y^2 * t^2 + 4 * y^3 * t^2 - 24 * x * y * t^3 +
        24 * y^2 * t^3 +
        20 * x^2 * y * t - 20 * x * y^2 * t - 160 * x * y * t^2 +
        96 * y^2 * t^2 +
        128 * x * t^3 +
        16 * x^2 * y +
        96 * x * y * t +
        2304 * x * t^2 +
        1152 * x * y +
        13824 * x * t +
        27648 * x
    f4 =
        -x^3 * z * t^2 + x^2 * z^2 * t^2 - 6 * x^3 * z * t +
        4 * x^2 * z^2 * t +
        32 * x^3 * t^2 - 72 * x^2 * z * t^2 - 87 * x * z^2 * t^2 - z^3 * t^2 -
        8 * x^3 * z - 432 * x^2 * z * t - 414 * x * z^2 * t +
        2592 * x * z * t^2 +
        864 * z^2 * t^2 - 1728 * x^2 * z - 20736 * x * z * t + 3456 * z^2 * t -
        186624 * z * t^2 - 124416 * x * z - 1492992 * z * t - 2985984 * z
    return ideal(R, [f1, f2, f3, f4])
end

function oberfr()
    dim = 3
    ve = [9, 9, 8]
    example = "Oberfranz"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x1, x2, x3) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x1", "x2", "x3"],
        ordering = Singular.ordering_M(StartOrd),
    )
    f1 = x2^3 + x1 * x2 * x3 + x2^2 * x3 + x1 * x3^3
    f2 = 3 + x1 * x2 + x1^2 * x2 + x2^2 * x3

    return ideal(R, [f1, f2])
end

function trinks1()
    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "trinks1"
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    R, (x, y, z, t, u, v) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["x", "y", "z", "t", "u", "v"],
        ordering = Singular.ordering_M(StartOrd),
    )

    f1 = 45 * y + 35 * u - 165 * v - 36
    f2 = 36 * y + 25 * z + 40 * t - 27 * u
    f3 = 25 * y * u - 165 * v^2 + 15 * x - 18 * z + 30t
    f4 = 15 * y * z + 20 * t * u - 9 * x
    f5 = -11 * v^3 + x * y + 2 * z * t
    f6 = -11 * u * v + 3 * v^2 + 99 * x


    return ideal(R, [f1, f2, f3, f4, f5, f6])
end

function ex1()
    dim = 4
    ve = [1, 1, 1, 1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (a, b, c, d) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["a", "b", "c", "d"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    return ideal(
        R,
        [
            d^2 +
            3 * a^2 * d +
            3 * c * d^2 +
            5 * a^3 * c +
            5 * a^2 * c^2 +
            4 * a^2 * c * d,
            3 + 2 * b + 3 * a * b + 3 * c^2 * d + 5 * a^4 + a * b * c^2,
        ],
    )
end

function ex2()
    dim = 4
    ve = [1, 1, 1, 1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (a, b, c, d) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["a", "b", "c", "d"],
        ordering = Singular.ordering_M(StartOrd),
    )

    S = change_order(R, TarOrd)
    return ideal(
        R,
        [
            2 * a^2 * b +
            3 * a * b^2 +
            3 * b^3 +
            4 * c^3 +
            4 * a * b * d +
            c^2 * d +
            2 * b * d^2 +
            2 * d^3 +
            4 * c^2 +
            2 * c * d +
            2 * c,
            2 * a^2 * b +
            5 * a^2 * c +
            2 * b * c^2 +
            c^2 * d +
            a * c +
            2 * b * d +
            2 * c * d +
            d^2 +
            a +
            4 * d,
        ],
    )
end

function v3g7d50()

    dim = 3
    ve = [1, 1, 1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (a, b, c) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["a", "b", "c"],
        ordering = Singular.ordering_M(StartOrd),
    )

    S = change_order(R, TarOrd)
    return ideal(R,[5*a^7+4*a^5*b^2+2*a^3*b^4+5*a^2*b^5+4*b^7+5*a^3*b^3*c+3*a^2*b^4*c+4*a^5*c^2+5*a^2*b^3*c^2+5*a^4*c^3+a^3*b*c^3+2*a^2*b^2*c^3+3*a*b^3*c^3+4*a^2*c^5+5*a*b*c^5+b^2*c^5+a^6+3*a^4*b^2+4*a*b^5+3*a^5*c+3*a*b^4*c+5*a^3*b*c^2+2*a^3*c^3+5*a^2*b*c^3+2*a^2*c^4+2*b^2*c^4+a*c^5+2*b*c^5+2*c^6+3*a^3*b^2+a*b^4+5*a*b^3*c+4*a^2*b*c^2+a*b^2*c^2+a^2*c^3+5*a*c^4+4*c^5+2*a^3*b+2*b^4+4*a*b^2*c+3*a*b*c^2+5*b*c^3+2*a*b^2+5*b^2*c+5*a*c^2+5*c^3+3*a^2+a*c+b*c+2*c^2+a+2*c,4*a^6*b+4*a^5*b^2+3*a^4*b^3+4*a^3*b^4+5*a*b^6+4*a^6*c+4*a^4*b^2*c+5*a^4*b*c^2+a^3*b^2*c^2+4*a^2*b^3*c^2+a^4*c^3+5*a^3*b*c^3+4*b^4*c^3+3*a^2*b*c^4+a*b^2*c^4+3*b^3*c^4+4*a^2*c^5+5*a*b*c^5+3*a*c^6+c^7+3*a^5*b+2*a^3*b^3+5*a^2*b^4+3*b^6+a^5*c+5*a^4*b*c+4*a^2*b^3*c+5*b^5*c+a^4*c^2+4*a^3*b*c^2+4*a^2*b^2*c^2+4*a^3*c^3+4*a^2*b*c^3+4*b^3*c^3+3*a^2*c^4+a^5+5*a^2*b^3+4*a*b^4+5*b^5+5*a^4*c+5*a^2*b^2*c+4*b^4*c+4*a^3*c^2+3*a^2*b*c^2+5*b^3*c^2+3*a^2*c^3+4*b^2*c^3+3*c^5+5*a^4+2*a^3*b+b^4+2*a^3*c+3*a^2*b*c+4*a*b^2*c+2*b^3*c+4*a^2*c^2+3*b^2*c^2+2*c^4+3*a^3+3*b^3+a*b*c+2*b^2*c+2*a*c^2+c^3+2*b^2+4*a*c+5*a+5])
end

function v4g3()

    dim = 4
    ve = [1, 1, 1, 1]
    StartOrd = ordering_as_matrix(:degrevlex, dim)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (a, b, c, d) = polynomial_ring(
        Singular.N_ZpField(32003),
        ["a", "b", "c", "d"],
        ordering = Singular.ordering_M(StartOrd),
    )

    S = change_order(R, TarOrd)
    return ideal(R,[3*a^3+3*a^2*b+5*a*b^2+2*b^3+a^2*c+2*b^2*c+3*a*c^2+3*a^2*d+5*b^2*d+3*a*c*d+a*d^2+b*d^2+3*c*d^2+5*d^3+4*b^2+3*a*c+5*b*c+3*c^2+4*a*d+2*b*d+c*d+d^2+4*a+4*b+c,3*a^2*b+5*a*b^2+2*a^2*c+5*a*b*c+4*b^2*c+4*a*c^2+5*c^3+5*a^2*d+a*b*d+4*b^2*d+2*a*c*d+3*b*c*d+4*c^2*d+5*b*d^2+4*a*b+3*b^2+2*a*c+5*c^2+5*a*d+2*d^2+4*a+5*c
])
end
