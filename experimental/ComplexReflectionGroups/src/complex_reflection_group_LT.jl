
###########################################################################################
# Models from Lehrer & Taylor (2009)
###########################################################################################
function complex_reflection_group_LT(n::Int)

    if n == 4
        # See Lehrer & Taylor (2009), page 85-86
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r1prime = ω//2 * matspace(K[-1-i -1+i ; 1+i -1+i])

        push!(gens, transpose(r1))
        push!(gens, transpose(r1prime))

        W = matrix_group(gens)

    elseif n == 5
         # See Lehrer & Taylor (2009), page 85-86
         R,x = polynomial_ring(QQ)
         K,i = number_field(x^2+1, "i")
         R,x = polynomial_ring(K)
         K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
         i = K(i)
         matspace = matrix_space(K, 2, 2)
         gens = elem_type(matspace)[]

         r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
         r2prime = ω//2 * matspace(K[-1+i 1-i ; -1-i -1-i])

         push!(gens, transpose(r1))
         push!(gens, transpose(r2prime))

         W = matrix_group(gens)

    elseif n == 6
        # See Lehrer & Taylor (2009), page 85-86
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r = matspace(K[1 0 ; 0 -1])
        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
    
        push!(gens, transpose(r))
        push!(gens, transpose(r1))

        W = matrix_group(gens)

    elseif n == 7
        # See Lehrer & Taylor (2009), page 85-86
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r = matspace(K[1 0 ; 0 -1])
        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r2 = ω//2 * matspace(K[-1+i -1+i ; 1+i -1-i])
    
        push!(gens, transpose(r))
        push!(gens, transpose(r1))
        push!(gens, transpose(r2))

        W = matrix_group(gens)

    elseif n == 8
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r4 = matspace(K[1 0 ; 0 i])
        r4prime = 1//2 * matspace(K[1+i -1+i ; -1+i 1+i])
    
        push!(gens, transpose(r4))
        push!(gens, transpose(r4prime))

        W = matrix_group(gens)

    elseif n == 9
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,sqrt2 = number_field(x^2-2, "√2")
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
        r4 = matspace(K[1 0 ; 0 i])
    
        push!(gens, transpose(r3))
        push!(gens, transpose(r4))

        W = matrix_group(gens)

    elseif n == 10
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r4prime = 1//2 * matspace(K[1+i -1+i ; -1+i 1+i])
    
        push!(gens, transpose(r1))
        push!(gens, transpose(r4prime))

        W = matrix_group(gens)

    elseif n == 11
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
        R,x = polynomial_ring(K)
        K,sqrt2 = number_field(x^2-2, "√2")
        i = K(i)
        ω = K(ω)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
        r4 = matspace(K[1 0 ; 0 i])

        push!(gens, transpose(r1))
        push!(gens, transpose(r3))
        push!(gens, transpose(r4))

        W = matrix_group(gens)

    elseif n == 12
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,sqrt2 = number_field(x^2-2, "√2")
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
        r3prime = 1//sqrt2 * matspace(K[1 1 ; 1 -1])
        r3primeprime = 1//sqrt2 * matspace(K[0 1+i ; 1-i 0])

        push!(gens, transpose(r3))
        push!(gens, transpose(r3prime))
        push!(gens, transpose(r3primeprime))

        W = matrix_group(gens)

    elseif n == 13
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,sqrt2 = number_field(x^2-2, "√2")
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r = matspace(K[1 0 ; 0 -1])
        r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
        r3primeprime = 1//sqrt2 * matspace(K[0 1+i ; 1-i 0])

        push!(gens, transpose(r))
        push!(gens, transpose(r3))
        push!(gens, transpose(r3primeprime))

        W = matrix_group(gens)

    elseif n == 14
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
        R,x = polynomial_ring(K)
        K,sqrt2 = number_field(x^2-2, "√2")
        i = K(i)
        ω = K(ω)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r3prime = 1//sqrt2 * matspace(K[1 1 ; 1 -1])

        push!(gens, transpose(r1))
        push!(gens, transpose(r3prime))

        W = matrix_group(gens)

    elseif n == 15
        # See Lehrer & Taylor (2009), page 87-88
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is primitive 3rd root of unity
        R,x = polynomial_ring(K)
        K,sqrt2 = number_field(x^2-2, "√2")
        i = K(i)
        ω = K(ω)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r = matspace(K[1 0 ; 0 -1])
        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])

        push!(gens, transpose(r))
        push!(gens, transpose(r1))
        push!(gens, transpose(r3))

        W = matrix_group(gens)

    elseif n == 16
        # See Lehrer & Taylor (2009), page 89-90
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,τ = number_field(x^2-x-1, "τ") #τ is golden ratio
        i = K(i)
        R,x = polynomial_ring(K)
        K,ζ = number_field(x^2 + (-τ + 1)*x + 1, "ζ") #ζ is prim 5th root of uni
        i = K(i)
        τ = K(τ)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])
        r5prime = -ζ^2//2 * matspace(K[τ-i -τ+1 ; τ-1 τ+i])

        push!(gens, transpose(r5))
        push!(gens, transpose(r5prime))

        W = matrix_group(gens)

    elseif n == 17
        # See Lehrer & Taylor (2009), page 89-90
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,τ = number_field(x^2-x-1, "τ") #τ is golden ratio
        i = K(i)
        R,x = polynomial_ring(K)
        K,ζ = number_field(x^2 + (-τ + 1)*x + 1, "ζ") #ζ is prim 5th root of uni
        i = K(i)
        τ = K(τ)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])
        r = matspace(K[1 0 ; 0 -1])

        push!(gens, transpose(r5))
        push!(gens, transpose(r))

        W = matrix_group(gens)

    elseif n == 18
        # See Lehrer & Taylor (2009), page 89-90
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,τ = number_field(x^2-x-1, "τ") #τ is golden ratio
        i = K(i)
        R,x = polynomial_ring(K)
        K,ζ = number_field(x^2 + (-τ + 1)*x + 1, "ζ") #ζ is prim 5th root of uni
        i = K(i)
        τ = K(τ)
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is prim 3rd root of unity
        i = K(i)
        τ = K(τ)
        ζ = K(ζ)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])

        push!(gens, transpose(r1^2))
        push!(gens, transpose(r5))

        W = matrix_group(gens)

    elseif n == 19
        # See Lehrer & Taylor (2009), page 89-90
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,τ = number_field(x^2-x-1, "τ") #τ is golden ratio
        i = K(i)
        R,x = polynomial_ring(K)
        K,ζ = number_field(x^2 + (-τ + 1)*x + 1, "ζ") #ζ is prim 5th root of uni
        i = K(i)
        τ = K(τ)
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is prim 3rd root of unity
        i = K(i)
        τ = K(τ)
        ζ = K(ζ)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r = matspace(K[1 0 ; 0 -1])
        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])

        push!(gens, transpose(r))
        push!(gens, transpose(r1))
        push!(gens, transpose(r5))

        W = matrix_group(gens)

    elseif n == 20
        # See Lehrer & Taylor (2009), page 89-90
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,τ = number_field(x^2-x-1, "τ") #τ is golden ratio
        i = K(i)
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is prim 3rd root of unity
        i = K(i)
        τ = K(τ)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
        r1primeprime = ω//2 * matspace(K[-τ*i-1 (-τ+1)*i ; (-τ+1)*i τ*i-1])

        push!(gens, transpose(r1))
        push!(gens, transpose(r1primeprime))

        W = matrix_group(gens)

    elseif n == 21
        # See Lehrer & Taylor (2009), page 89-90
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,τ = number_field(x^2-x-1, "τ") #τ is golden ratio
        i = K(i)
        R,x = polynomial_ring(K)
        K,ω = number_field(x^2+x+1, "ω") #ω is prim 3rd root of unity
        i = K(i)
        τ = K(τ)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r = matspace(K[1 0 ; 0 -1])
        r1primeprime = ω//2 * matspace(K[-τ*i-1 (-τ+1)*i ; (-τ+1)*i τ*i-1])

        push!(gens, transpose(r))
        push!(gens, transpose(r1primeprime))

        W = matrix_group(gens)

    elseif n == 22
        # See Lehrer & Taylor (2009), page 89-90
        R,x = polynomial_ring(QQ)
        K,i = number_field(x^2+1, "i")
        R,x = polynomial_ring(K)
        K,τ = number_field(x^2-x-1, "τ") #τ is golden ratio
        i = K(i)
        matspace = matrix_space(K, 2, 2)
        gens = elem_type(matspace)[]

        r = matspace(K[1 0 ; 0 -1])
        rprime = matspace(K[0 1 ; 1 0])
        rprimeprime = 1//2 * matspace(K[τ (τ-1)*i+1 ; (-τ+1)*i+1 -τ])

        push!(gens, transpose(r))
        push!(gens, transpose(rprime))
        push!(gens, transpose(rprimeprime))

        W = matrix_group(gens)

    elseif n == 23
        # See Lehrer & Taylor (2009), page 110
        # This group is realized by the line system H3
        R,x = polynomial_ring(QQ)
        K,τ = number_field(x^2-x-1, "τ")
        V = vector_space(K,3)
        matspace = matrix_space(K, 3, 3)

        # The line system H3
        cycshift = elem_type(matspace)[]
        push!(cycshift, identity_matrix(K, 3))
        push!(cycshift, matrix(K,3,3,[0 1 0; 0 0 1; 1 0 0])) #(2 3 1)
        push!(cycshift, matrix(K,3,3,[0 0 1; 1 0 0; 0 1 0])) #(3 1 2) 

        lines = [V([1,0,0])*A for A in cycshift]
        lines = vcat(lines, [1//2*V([K(1),τ,τ^-1])*A for A in cycshift])
        lines = vcat(lines, [1//2*V([K(1),τ,-τ^-1])*A for A in cycshift])
        lines = vcat(lines, [1//2*V([K(1),-τ,τ^-1])*A for A in cycshift])
        lines = vcat(lines, [1//2*V([K(1),-τ,-τ^-1])*A for A in cycshift])

        refls = Set([unitary_reflection(l) for l in lines])

        # I guessed the roots of generating reflections from the proof of 
        # Theorem 8.10 (page 144). 
        # The line system H3 is the star-closure of the following lines.
        a = V([0,1,0])
        b = 1//2*V([K(1),τ,τ^-1])
        p = V([1,0,0])  

        # Now, create the corresponding generators
        gens = elem_type(matspace)[]

        push!(gens, matrix(unitary_reflection(a)))
        push!(gens, matrix(unitary_reflection(b)))
        push!(gens, matrix(unitary_reflection(p)))

        W = matrix_group(gens)

        set_attribute!(W, :complex_reflections, refls)

    else

        error("LT model not yet implemented")

    end

    return W

end

function complex_reflection_group_LT(t::Tuple)
    (m,p,n) = t

    # See Lehrer & Taylor (2009), p 35-36
    
    # For the symmetric group case we do not want permutation matrices
    # since this is not an irreducible reflection representation
    if m == 1 && p == 1
        error("LT model not Specified for symmetric group case")
    end

    # See Lehrer & Taylor (2009), p 35-36
    K, z = cyclotomic_field(m)
    matspace = matrix_space(K, n, n)

    #the transpositions
    transp = elem_type(matspace)[]
    for i = 1:n-1
        r = matspace(1)
        r[i,i] = 0
        r[i+1,i+1] = 0
        r[i,i+1] = 1
        r[i+1,i] = 1
        push!(transp, r)
    end

    #the matrix t
    t = matspace(1)
    t[1,1] = z
                                                    
    #the matrix s
    if n > 1 
        s = t^-1*transp[1]*t
    end

    #now, create the list of generators
    if m == p
        gens = vcat([s], transp)
    elseif p == 1
        gens = vcat([t], transp)
    else
        gens = vcat([s], [t^p], transp)
    end

    W = matrix_group(gens)

    return W

end
