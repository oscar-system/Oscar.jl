using Oscar

function katsura6(k::Field=QQ)
    R, (x1, x2, x3, x4, x5, x6, x7) = k[:x1, :x2, :x3, :x4, :x5, :x6, :x7]
    F = [
        1*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7-1,
        2*x4*x3+2*x5*x2+2*x6*x1+2*x7*x2-1*x6,
        1*x3^2+2*x4*x2+2*x5*x1+2*x6*x2+2*x7*x3-1*x5,
        2*x3*x2+2*x4*x1+2*x5*x2+2*x6*x3+2*x7*x4-1*x4,
        1*x2^2+2*x3*x1+2*x4*x2+2*x5*x3+2*x6*x4+2*x7*x5-1*x3,
        2*x2*x1+2*x3*x2+2*x4*x3+2*x5*x4+2*x6*x5+2*x7*x6-1*x2,
        1*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2-1*x1
    ];

    return ideal(F), lex(R), default_ordering(R)
end