
function simple(io)
    R,(x,y) = QQ[:x,:y]
    I = ideal([y^4+ x^3-x^2+x,x^4])
    
    print(io, "simple,")
    t = @belapsed groebner_walk($I; algorithm=:standard)
    print(io, t, ",")
    t = @belapsed groebner_walk($I; algorithm=:generic)
    print(io, t, ",")
    t = @belapsed groebner_basis($I; ordering=lex($R))
    println(io, t)
end