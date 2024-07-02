#Newell's teapot, patch 1
#a 2 dimensional ideal from tran 2004
using Oscar

function newellp1(k::Field)
    R, (x,y,z,u,v) = polynomial_ring(QQ, ["x","y","z","u","v"])
    
    F = [
        -x + 7//5 - 231//125 * v^2 + 39//80 * u^2 - 1//5 * u^3 + 99//400 * u * v^2 - 1287//2000 * u^2 * v^2 + 33//125 * u^3 * v^2 - 3//16 * u + 56//125 * v^3 - 3//50 * u * v^3 + 39//250 * u^2 * v^3 - 8//125 * u^3 * v^3,
        -y + 63//125 * v^2 - 294//125 * v + 56//125 * v^3 - 819//1000 * u^2 * v + 42//125 * u^3 * v - 3//50 * u * v^3 + 351//2000 * u^2 * v^2 + 39//250 * u^2 * v^3 - 9//125 * u^3 * v^2 - 8//125 * u^3 * v^3,
        -z + 12//5 - 63//160 * u^2 + 63//160 * u
    ]

    if k == QQ
        o1 = matrix_ordering(R, [1 1 1 0 0; 0 0 0 1 1; 0 0 0 1 0; 1 1 0 0 0; 1 0 0 0 0])
        o2 = matrix_ordering(R, [0 0 0 1 1; 1 1 1 0 0; 1 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0])

        return ideal(F), o2, o1
    end

    lcm_denom = [
        lcm(numerator.(coefficients(f))) for f in F
    ]
    integral_F = [
        change_coefficient_ring(
            k, l*f
        ) for (l, f) in zip(lcm_denom, F)
    ]

    R = parent(first(F))
    o1 = matrix_ordering(R, [1 1 1 0 0; 0 0 0 1 1; 0 0 0 1 0; 1 1 0 0 0; 1 0 0 0 0])
    o2 = matrix_ordering(R, [0 0 0 1 1; 1 1 1 0 0; 1 1 0 0 0; 1 0 0 0 0; 0 0 0 1 0])
    return ideal(integral_F), o2, o1
end
