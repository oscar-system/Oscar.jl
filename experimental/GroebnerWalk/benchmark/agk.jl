using Oscar

function agk4(k::Field)
    R, (x,y,z,u,v) = polynomial_ring(k, [:x,:y,:z,:u,:v])

    o1 = weight_ordering([1,1,1,0,0], degrevlex(R))
    o2 = weight_ordering([0,0,0,1,1], degrevlex(R))

    F = [
        u + u^2 - 2*v - 2*u^2*v + 2*u*v^2 - x,
        -6*u + 2*v + v^2 - 5*v^3 + 2*u*v^2 - 4*u^2*v^2 - y,
        -2 + 2*u^2 + 6*v - 3*u^2*v^2 - z
    ]

    return ideal(F), o2, o1
end

