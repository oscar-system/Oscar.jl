# This file implements explicit models of complex reflection groups.
#
# References:
#
# * Lehrer, G. I., & Taylor, D. E. (2009). Unitary reflection groups (Vol. 20, p. viii). Cambridge University Press, Cambridge.
# 
# * Thiel, U. (2014). On restricted rational Cherednik algebras. TU Kaiserslautern.
#
# * Marin, I., & Michel, J. (2010). Automorphisms of complex reflection groups. Represent. Theory, 14, 747–788.
#
# Ulrich Thiel, 2023 

export complex_reflection_group
export is_complex_reflection_group
export complex_reflection_group_type
export complex_reflection_group_model
export complex_reflection_group_dual

function complex_reflection_group(G::ComplexReflectionGroupType; model=:CHEVIE)
    
    # this will be the list of matrix groups corresponding to the components of G
    component_groups = MatrixGroup[]

    # list of models
    modellist = []

    for C in components(G)

        t = C.type[1]
        gens = nothing
        Cmodel = nothing

         # Set default model (CHEVIE)
         if model === nothing
            Cmodel = :CHEVIE
        else
            Cmodel = model
        end

        # we now create the list "gens" of matrix generators for C in the selected model
        if isa(t, Int)

            # Exceptional groups

            # G4
            if t == 4

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 85-86
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    i = K(i)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r1prime = ω//2 * matspace(K[-1-i -1+i ; 1+i -1+i])

                    push!(gens, r1)
                    push!(gens, r1prime)
                
                    
                elseif Cmodel == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 ; -zeta_3 - 1 1])
                    s2 = matspace(K[1 zeta_3 + 1 ; 0 zeta_3])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif Cmodel == :CHEVIE
                    K,z = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 z])
                    s2 = 1//3 * matspace(K[(2*z + 1) (2*z - 2) ; (z - 1) (z + 2)])

                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G5
            elseif t == 5

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 85-86
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    i = K(i)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r2prime = ω//2 * matspace(K[-1+i 1-i ; -1-i -1-i])

                    push!(gens, r1)
                    push!(gens, r2prime)

                elseif Cmodel == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 ; -2*zeta_3 - 2 1])
                    s2 = matspace(K[1 zeta_3 + 1 ; 0 zeta_3])

                    push!(gens, s1)
                    push!(gens, s2)


                elseif Cmodel == :CHEVIE
                    K,z = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 z])
                    s2 = 1//3 * matspace(K[(z + 2) (-2*z + 2) ; (-z + 1) (2*z + 1)])

                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G6
            elseif t == 6

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 85-86
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    i = K(i)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r = matspace(K[1 0 ; 0 -1])
                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                
                    push!(gens, r)
                    push!(gens, r1)

                elseif Cmodel == :Magma
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; -zeta_12^3 - 2*zeta_12^2 + 1 1])
                    s2 = matspace(K[1 zeta_12^2 ; 0 zeta_12^2 - 1])

                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G7
            elseif t == 7

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 85-86
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    i = K(i)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r = matspace(K[1 0 ; 0 -1])
                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r2 = ω//2 * matspace(K[-1+i -1+i ; 1+i -1-i])
                
                    push!(gens, r)
                    push!(gens, r1)
                    push!(gens, r2)

                elseif Cmodel == :Magma
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-zeta_12^2 - zeta_12 zeta_12^2 ; -2*zeta_12^2 - 2*zeta_12 zeta_12^2 + zeta_12])
                    s2 = matspace(K[1 0 ; -zeta_12^3 - zeta_12^2 + zeta_12 + 2 zeta_12^2 - 1])
                    s3 = matspace(K[zeta_12^2 - 1 0 ; zeta_12^3 + zeta_12^2 - zeta_12 - 2 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end   
                
            # G8
            elseif t == 8

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 87-88
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r4 = matspace(K[1 0 ; 0 i])
                    r4prime = 1//2 * matspace(K[1+i -1+i ; -1+i 1+i])
                
                    push!(gens, r4)
                    push!(gens, r4prime)

                elseif Cmodel == :Magma
                    K,zeta_4 = cyclotomic_field(4)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_4 0 ; -zeta_4 1])
                    s2 = matspace(K[1 1 ; 0 zeta_4])

                    push!(gens, s1)
                    push!(gens, s2)

                end   

            # G9
            elseif t == 9

                if Cmodel == :LT
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
                
                    push!(gens, r3)
                    push!(gens, r4)

                elseif Cmodel == :Magma
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_8^3 + zeta_8^2 + zeta_8 1])
                    s2 = matspace(K[1 -zeta_8^3 ; 0 -zeta_8^2])

                    push!(gens, s1)
                    push!(gens, s2)

                end  

            # G10
            elseif t == 10

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 87-88
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    i = K(i)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r4prime = 1//2 * matspace(K[1+i -1+i ; -1+i 1+i])
                
                    push!(gens, r1)
                    push!(gens, r4prime)

                elseif Cmodel == :Magma
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-zeta_12^3 0 ; -zeta_12^2 + zeta_12 1])
                    s2 = matspace(K[1 zeta_12^2 ; 0 zeta_12^2 - 1])

                    push!(gens, s1)
                    push!(gens, s2)

                end  

            # G11
            elseif t == 11

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 87-88
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    R,x = polynomial_ring(K)
                    K,sqrt2 = number_field(x^2-2, "√2")
                    i = K(i)
                    ω = K(ω)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])
                    r4 = matspace(K[1 0 ; 0 i])

                    push!(gens, r1)
                    push!(gens, r3)
                    push!(gens, r4)

                elseif Cmodel == :Magma
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-zeta_24^6 0 ; -zeta_24^6 - zeta_24^3 - 1 1])
                    s2 = matspace(K[-zeta_24^6 - zeta_24^3 zeta_24^6 ; -2*zeta_24^6 - 2*zeta_24^3 - 1 zeta_24^6 + zeta_24^3])
                    s3 = matspace(K[-zeta_24^7 zeta_24^7 ; -zeta_24^7 - zeta_24^4 - zeta_24 zeta_24^7 + zeta_24^4])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G12
            elseif t == 12

                if Cmodel == :LT
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

                    push!(gens, r3)
                    push!(gens, r3prime)
                    push!(gens, r3primeprime)

                elseif Cmodel == :Magma
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; -zeta_8^3 - zeta_8 + 1 1])
                    s2 = matspace(K[1 zeta_8^3 + zeta_8 + 1 ; 0 -1])
                    s3 = matspace(K[zeta_8^3 + zeta_8 - 1 -2 ; -zeta_8^3 - zeta_8 - 1 -zeta_8^3 - zeta_8 + 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G13
            elseif t == 13

                if Cmodel == :LT
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

                    push!(gens, r)
                    push!(gens, r3)
                    push!(gens, r3primeprime)

                elseif Cmodel == :Magma
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_8^3 - zeta_8 - 2 1])
                    s2 = matspace(K[zeta_8^3 + zeta_8^2 -zeta_8^2 ; 2*zeta_8^3 + 2*zeta_8^2 - 1 -zeta_8^3 - zeta_8^2])
                    s3 = matspace(K[-zeta_8^3 + zeta_8 + 1 -1 ; -2*zeta_8^3 + 2*zeta_8 + 2 zeta_8^3 - zeta_8 - 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G14
            elseif t == 14

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 87-88
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    R,x = polynomial_ring(K)
                    K,sqrt2 = number_field(x^2-2, "√2")
                    i = K(i)
                    ω = K(ω)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r3prime = 1//sqrt2 * matspace(K[1 1 ; 1 -1])

                    push!(gens, r1)
                    push!(gens, r3prime)

                elseif Cmodel == :Magma
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_24^4 - 1 0 ; -zeta_24^4 1])
                    s2 = matspace(K[1 -zeta_24^5 + 2*zeta_24^4 - zeta_24^3 + zeta_24 - 1 ; 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)

                end 

            # G15
            elseif t == 15

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 87-88
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    R,x = polynomial_ring(K)
                    K,sqrt2 = number_field(x^2-2, "√2")
                    i = K(i)
                    ω = K(ω)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r = matspace(K[1 0 ; 0 -1])
                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r3 = 1//sqrt2 * matspace(K[1 -1 ; -1 -1])

                    push!(gens, r)
                    push!(gens, r1)
                    push!(gens, r3)

                elseif Cmodel == :Magma
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_24^5 - zeta_24^3 - zeta_24 - 2 1])
                    s2 = matspace(K[-zeta_24^6 + zeta_24^4 + zeta_24^2 + zeta_24 -zeta_24 ; -2*zeta_24^6 + zeta_24^4 + 2*zeta_24^2 + 2*zeta_24 zeta_24^6 - zeta_24^2 - zeta_24])
                    s3 = matspace(K[-zeta_24^5 + zeta_24^3 + zeta_24 + 1 -1 ; -2*zeta_24^5 + 2*zeta_24^3 + 2*zeta_24 + 2 zeta_24^5 - zeta_24^3 - zeta_24 - 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G16
            elseif t == 16

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 89-90
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ζ = number_field(x^4+x^3+x^2+x+1, "ζ")
                    i = K(i)
                    τ = ζ + ζ^-1 + 1
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])
                    r5prime = -ζ^2//2 * matspace(K[τ-i -τ+1 ; τ-1 τ+i])

                    push!(gens, r5)
                    push!(gens, r5prime)
                
                elseif Cmodel == :Magma
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_5 0 ; 1 1])
                    s2 = matspace(K[1 -zeta_5 ; 0 zeta_5])

                    push!(gens, s1)
                    push!(gens, s2)

                end 

            # G17
            elseif t == 17

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 89-90
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ζ = number_field(x^4+x^3+x^2+x+1, "ζ")
                    i = K(i)
                    τ = ζ + ζ^-1 + 1
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])
                    r = matspace(K[1 0 ; 0 -1])

                    push!(gens, r5)
                    push!(gens, r)

                elseif Cmodel == :Magma
                    K,zeta_20 = cyclotomic_field(20)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; -zeta_20^6 + zeta_20^5 + zeta_20^4 - 2*zeta_20^2 + 1 1])
                    s2 = matspace(K[1 zeta_20^2 ; 0 zeta_20^4])

                    push!(gens, s1)
                    push!(gens, s2)

                end 

            # G18
            elseif t == 18

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 89-90
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    R,x = polynomial_ring(K)
                    K,ζ = number_field(x^4+x^3+x^2+x+1, "ζ")
                    i = K(i)
                    ω = K(ω)
                    τ = ζ + ζ^-1 + 1
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])

                    push!(gens, r1^2)
                    push!(gens, r5)

                elseif Cmodel == :Magma
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_15^5 0 ; zeta_15^7 - 2*zeta_15^6 - zeta_15^3 + zeta_15^2 - zeta_15 - 1 1])
                    s2 = matspace(K[1 -zeta_15^7 + zeta_15^6 + zeta_15^3 - zeta_15^2 + 1 ; 0 zeta_15^3])

                    push!(gens, s1)
                    push!(gens, s2)

                end
                
            # G19
            elseif t == 19

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 89-90
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    R,x = polynomial_ring(K)
                    K,ζ = number_field(x^4+x^3+x^2+x+1, "ζ")
                    i = K(i)
                    ω = K(ω)
                    τ = ζ + ζ^-1 + 1
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r = matspace(K[1 0 ; 0 -1])
                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r5 = ζ^2//2 * matspace(K[-τ+i -τ+1 ; τ-1 -τ-i])

                    push!(gens, r)
                    push!(gens, r1)
                    push!(gens, r5)

                elseif Cmodel == :Magma
                    K,zeta_60 = cyclotomic_field(60)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_60^12 0 ; zeta_60^12 + zeta_60^10 - 1 1])
                    s2 = matspace(K[1 -1 ; 0 zeta_60^10 - 1])
                    s3 = matspace(K[-zeta_60^13 - zeta_60^11 + zeta_60^3 + zeta_60 zeta_60^13 - zeta_60^3 ; -zeta_60^15 - 2*zeta_60^13 - 2*zeta_60^11 + zeta_60^5 + zeta_60^3 + zeta_60 zeta_60^13 + zeta_60^11 - zeta_60^3 - zeta_60])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end
                
            # G20
            elseif t == 20

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 89-90
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    R,x = polynomial_ring(K)
                    K,ζ = number_field(x^4+x^3+x^2+x+1, "ζ")
                    i = K(i)
                    ω = K(ω)
                    τ = ζ + ζ^-1 + 1
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r1 = ω//2 * matspace(K[-1-i 1-i ; -1-i -1+i])
                    r1primeprime = ω//2 * matspace(K[-τ*i-1 (-τ+1)*i ; (-τ+1)*i τ*i-1])

                    push!(gens, r1)
                    push!(gens, r1primeprime)

                elseif Cmodel == :Magma
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_15^5 0 ; -2*zeta_15^5 + zeta_15^4 + zeta_15 - 2 1])
                    s2 = matspace(K[1 zeta_15^5 + 1 ; 0 zeta_15^5])

                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G21
            elseif t == 21

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 89-90
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ω = number_field(x^2+x+1, "ω")
                    R,x = polynomial_ring(K)
                    K,ζ = number_field(x^4+x^3+x^2+x+1, "ζ")
                    i = K(i)
                    ω = K(ω)
                    τ = ζ + ζ^-1 + 1
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r = matspace(K[1 0 ; 0 -1])
                    r1primeprime = ω//2 * matspace(K[-τ*i-1 (-τ+1)*i ; (-τ+1)*i τ*i-1])

                    push!(gens, r)
                    push!(gens, r1primeprime)

                elseif Cmodel == :Magma
                    K,zeta_60 = cyclotomic_field(60)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_60^11 - 2*zeta_60^10 + zeta_60^9 - zeta_60 + 1 1])
                    s2 = matspace(K[1 zeta_60^10 ; 0 zeta_60^10 - 1])

                    push!(gens, s1)
                    push!(gens, s2)

                end
            
            # G22
            elseif t == 22

                if Cmodel == :LT
                    # See Lehrer & Taylor (2009), page 89-90
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,ζ = number_field(x^4+x^3+x^2+x+1, "ζ")
                    i = K(i)
                    τ = ζ + ζ^-1 + 1
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r = matspace(K[1 0 ; 0 -1])
                    rprime = matspace(K[0 1 ; 1 0])
                    rprimeprime = 1//2 * matspace(K[τ (τ-1)*i+1 ; (-τ+1)*i+1 -τ])

                    push!(gens, r)
                    push!(gens, rprime)
                    push!(gens, rprimeprime)

                elseif Cmodel == :Magma
                    K,zeta_20 = cyclotomic_field(20)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_20^6 + zeta_20^5 - zeta_20^4 1])
                    s2 = matspace(K[1 zeta_20^6 - zeta_20^5 - zeta_20^4 ; 0 -1])
                    s3 = matspace(K[zeta_20^7 - zeta_20^6 - zeta_20^5 + zeta_20^4 + zeta_20^3 - 1 -2*zeta_20^6 + 2*zeta_20^4 ; -zeta_20^7 + 2*zeta_20^5 - zeta_20^3 + 1 -zeta_20^7 + zeta_20^6 + zeta_20^5 - zeta_20^4 - zeta_20^3 + 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end

            # G23
            elseif t == 23
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 ; -zeta_5^3 - zeta_5^2 1 0 ; 0 0 1])
                    s2 = matspace(K[1 -zeta_5^3 - zeta_5^2 0 ; 0 -1 0 ; 0 1 1])
                    s3 = matspace(K[1 0 0 ; 0 1 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                end
                
           
            # G24
            elseif t == 24
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_7 = cyclotomic_field(7)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 ; -1 1 0 ; -zeta_7^4 - zeta_7^2 - zeta_7 - 1 0 1])
                    s2 = matspace(K[1 -1 0 ; 0 -1 0 ; 0 -zeta_7^4 - zeta_7^2 - zeta_7 1])
                    s3 = matspace(K[1 0 zeta_7^4 + zeta_7^2 + zeta_7 ; 0 1 zeta_7^4 + zeta_7^2 + zeta_7 + 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                end

            # G25
            elseif t == 25
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 0 ; zeta_3 + 1 1 0 ; 0 0 1])
                    s2 = matspace(K[1 -zeta_3 - 1 0 ; 0 zeta_3 0 ; 0 -zeta_3 - 1 1])
                    s3 = matspace(K[1 0 0 ; 0 1 zeta_3 + 1 ; 0 0 zeta_3])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                end

            # G26
            elseif t == 26
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 0 ; zeta_3 + 1 1 0 ; 0 0 1])
                    s2 = matspace(K[1 -zeta_3 - 1 0 ; 0 zeta_3 0 ; 0 -zeta_3 + 1 1])
                    s3 = matspace(K[1 0 0 ; 0 1 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                end

            # G27
            elseif t == 27
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 ; -1 1 0 ; zeta_15^7 + zeta_15^4 - zeta_15^3 + zeta_15^2 + zeta_15 - 1 0 1])
                    s2 = matspace(K[1 -1 0 ; 0 -1 0 ; 0 zeta_15^7 + zeta_15^4 - zeta_15^3 + zeta_15^2 + zeta_15 1])
                    s3 = matspace(K[1 0 -zeta_15^4 - zeta_15 ; 0 1 -zeta_15^4 - zeta_15 + 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                end

            # G28
            elseif t == 28
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K = QQ
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 ; 1 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[1 1 0 0 ; 0 -1 0 0 ; 0 1 1 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 ; 0 1 2 0 ; 0 0 -1 0 ; 0 0 1 1])
                    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 1 ; 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)

                end

            # G29
            elseif t == 29
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_4 = cyclotomic_field(4)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 ; 1 1 0 0 ; zeta_4 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[1 1 0 0 ; 0 -1 0 0 ; 0 -zeta_4 + 1 1 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 -zeta_4 0 ; 0 1 zeta_4 + 1 0 ; 0 0 -1 0 ; 0 0 1 1])
                    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 1 ; 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                end

            # G30
            elseif t == 30
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 ; -zeta_5^3 - zeta_5^2 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[1 -zeta_5^3 - zeta_5^2 0 0 ; 0 -1 0 0 ; 0 1 1 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 ; 0 1 1 0 ; 0 0 -1 0 ; 0 0 1 1])
                    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 1 ; 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                end

            # G31
            elseif t == 31
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_4 = cyclotomic_field(4)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 ; 1 1 0 0 ; zeta_4 - 1 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[1 1 0 0 ; 0 -1 0 0 ; 0 -zeta_4 1 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 -zeta_4 - 1 0 ; 0 1 zeta_4 0 ; 0 0 -1 0 ; 0 0 1 1])
                    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 1 ; 0 0 0 -1])
                    s5 = matspace(K[-zeta_4 0 zeta_4 - 1 0 ; 0 1 0 0 ; -zeta_4 - 1 0 zeta_4 0 ; zeta_4 0 1 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)
                end

            # G32
            elseif t == 32
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 0 0 ; -zeta_3 - 1 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[1 zeta_3 + 1 0 0 ; 0 zeta_3 0 0 ; 0 zeta_3 + 1 1 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 ; 0 1 -zeta_3 - 1 0 ; 0 0 zeta_3 0 ; 0 0 -zeta_3 - 1 1])
                    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 zeta_3 + 1 ; 0 0 0 zeta_3])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                end

            # G3
            elseif t == 33
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 5, 5)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 0 ; 1 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ; 0 0 0 0 1])
                    s2 = matspace(K[1 1 0 0 0 ; 0 -1 0 0 0 ; 0 1 1 0 0 ; 0 -zeta_3 - 1 0 1 0 ; 0 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 0 ; 0 1 1 0 0 ; 0 0 -1 0 0 ; 0 0 1 1 0 ; 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 ; 0 1 0 zeta_3 0 ; 0 0 1 1 0 ; 0 0 0 -1 0 ; 0 0 0 1 1])
                    s5 = matspace(K[1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 1 ; 0 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)

                end

            # G34
            elseif t == 34
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 6, 6)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 0 0 ; 1 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s2 = matspace(K[1 1 0 0 0 0 ; 0 -1 0 0 0 0 ; 0 1 1 0 0 0 ; 0 -zeta_3 - 1 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 0 0 ; 0 1 1 0 0 0 ; 0 0 -1 0 0 0 ; 0 0 1 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 0 ; 0 1 0 zeta_3 0 0 ; 0 0 1 1 0 0 ; 0 0 0 -1 0 0 ; 0 0 0 1 1 0 ; 0 0 0 0 0 1])
                    s5 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 1 0 ; 0 0 0 0 -1 0 ; 0 0 0 0 1 1])
                    s6 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 1 ; 0 0 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)
                    push!(gens, s6)
                end

            # G35
            elseif t == 35
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K = QQ
                    matspace = matrix_space(K, 6, 6)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 1 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 0 ; 0 -1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 1 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s3 = matspace(K[1 0 1 0 0 0 ; 0 1 0 0 0 0 ; 0 0 -1 0 0 0 ; 0 0 1 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 0 ; 0 1 0 1 0 0 ; 0 0 1 1 0 0 ; 0 0 0 -1 0 0 ; 0 0 0 1 1 0 ; 0 0 0 0 0 1])
                    s5 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 1 0 ; 0 0 0 0 -1 0 ; 0 0 0 0 1 1])
                    s6 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 1 ; 0 0 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)
                    push!(gens, s6)
                end

            # G36
            elseif t == 36
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K = QQ
                    matspace = matrix_space(K, 7, 7)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 1 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 1 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s3 = matspace(K[1 0 1 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 -1 0 0 0 0 ; 0 0 1 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 1 0 0 0 ; 0 0 1 1 0 0 0 ; 0 0 0 -1 0 0 0 ; 0 0 0 1 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s5 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 1 0 0 ; 0 0 0 0 -1 0 0 ; 0 0 0 0 1 1 0 ; 0 0 0 0 0 0 1])
                    s6 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 1 0 ; 0 0 0 0 0 -1 0 ; 0 0 0 0 0 1 1])
                    s7 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 1 ; 0 0 0 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)
                    push!(gens, s6)
                    push!(gens, s7)

                end

            # G37
            elseif t == 37
                
                if Cmodel == :LT
                    nothing

                elseif Cmodel == :Magma
                    K = QQ
                    matspace = matrix_space(K, 8, 8)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 1 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 0 0 0 ; 0 -1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 1 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s3 = matspace(K[1 0 1 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 -1 0 0 0 0 0 ; 0 0 1 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 1 0 0 0 0 ; 0 0 1 1 0 0 0 0 ; 0 0 0 -1 0 0 0 0 ; 0 0 0 1 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s5 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 1 0 0 0 ; 0 0 0 0 -1 0 0 0 ; 0 0 0 0 1 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s6 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 1 0 0 ; 0 0 0 0 0 -1 0 0 ; 0 0 0 0 0 1 1 0 ; 0 0 0 0 0 0 0 1])
                    s7 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 1 0 ; 0 0 0 0 0 0 -1 0 ; 0 0 0 0 0 0 1 1])
                    s8 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 1 ; 0 0 0 0 0 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)
                    push!(gens, s6)
                    push!(gens, s7)
                    push!(gens, s8)
                end
            
            end


        else
            # The infinite series
            
            (m,p,n) = t

            # Symmetric group case needs special care
            if m == 1 && p == 1
                # Set default model for this case
                if model === nothing
                    Cmodel = :CHEVIE
                else
                    Cmodel = model
                end

                if Cmodel === :CHEVIE
                    # the generators correspond to the transpositions r_i=(i,i+1)
                    # but note that we consider the irreducible representation of 
                    # dimension n-1, so it's not simply the transposition matrix
                    matspace = matrix_space(QQ, n-1, n-1)
                    gens = elem_type(matspace)[]
                    for i = 1:n-1
                        ri = matspace(1)
                        ri[i,i] = -1
                        if i > 1
                            ri[i,i-1] = 1
                        end
                        if i < n-1
                            ri[i,i+1] = 1
                        end
                        push!(gens, ri)
                    end
                end
            else
                # Set default model for this case
                if model === nothing
                    Cmodel = :CHEVIE
                else
                    Cmodel = model
                end

                if Cmodel === :Magma
                    # See Lehrer & Taylor (2009), p 35-36
                    # What's implemented in Magma deviates slightly from Lehrer & Taylor, 
                    # see the remarks below.
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
                    t[n,n] = z #LT would take t[1,1]=z here

                    #the matrix s
                    if n > 1
                        s = t^-1*transp[n-1]*t #LT would take transp[1] here
                    end

                    #now, create the list of generators
                    #LT would order the transpositions to the end
                    if m == p
                        gens = vcat(transp, [s])
                    elseif p == 1
                        gens = vcat(transp, [t])
                    else
                        gens = vcat(transp, [s], [t^p])
                    end
                    
                elseif Cmodel === :CHEVIE
                    # I extracted this from Marin & Michel (2010), p 751. 
                    # This should match what is in CHEVIE.
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

                    #the matrix t'
                    tprime = t^p

                    #the matrix s_1'
                    s1prime = t^-1*transp[1]*t

                    #now, create the list of generators
                    if m == p
                        gens = vcat([sprime], transp)
                    elseif p == 1
                        gens = vcat([tprime], transp)
                    else
                        gens = vcat([tprime], [s1prime], transp)
                    end
                end
            end
        end

        # create matrix group from the generators
        if gens === nothing
            throw(ArgumentError("Specified model not found"))
        end

        matgrp = matrix_group(gens)

        # set attributes that are already known from type
        set_attribute!(matgrp, :order, order(C))
        set_attribute!(matgrp, :is_complex_reflection_group, true)
        set_attribute!(matgrp, :complex_reflection_group_type, C)
        set_attribute!(matgrp, :complex_reflection_group_model, [Cmodel])
        set_attribute!(matgrp, :is_irreducible, true)

        # add to list
        push!(component_groups, matgrp)
        push!(modellist, Cmodel)
    end

    if length(component_groups) == 1
        # treat this as a special case because direct_product([G]) will have type direct
        # product which looks weird for a single group.
        return component_groups[1]
    else
        matgrp = direct_product(component_groups)

        set_attribute!(matgrp, :order, order(G))
        set_attribute!(matgrp, :is_complex_reflection_group, true)
        set_attribute!(matgrp, :complex_reflection_group_type, G)
        set_attribute!(matgrp, :complex_reflection_group_model, modellist)
        set_attribute!(matgrp, :is_irreducible, false)
    end

end

# Convenience constructors
complex_reflection_group(i::Int; model=nothing) = complex_reflection_group(ComplexReflectionGroupType(i); model=model)

complex_reflection_group(m::Int, p::Int, n::Int; model=nothing) = complex_reflection_group(ComplexReflectionGroupType(m,p,n); model=model)

complex_reflection_group(X::Vector; model=nothing) = complex_reflection_group(ComplexReflectionGroupType(X); model=model)

function is_complex_reflection_group(G::MatrixGroup)
    if has_attribute(G, :is_complex_reflection_group)
        return get_attribute(G, :is_complex_reflection_group)
    end
    return false 
    #this should be upgraded later to work with a general matrix group
end

function complex_reflection_group_type(G::MatrixGroup)
    if has_attribute(G, :complex_reflection_group_type)
        return get_attribute(G, :complex_reflection_group_type)
    end
    return nothing
    # this should be upgraded later to work with a general matrix group (identifying the
    # type from scratch is not so easy though)
end

function complex_reflection_group_model(G::MatrixGroup)
    if has_attribute(G, :complex_reflection_group_model)
        return get_attribute(G, :complex_reflection_group_model)
    end
    return nothing
end

function complex_reflection_group_dual(W::MatrixGroup)

    WD = matrix_group([transpose(matrix(w^-1)) for w in gens(W)])

    set_attribute!(WD, :order, order(W))
    set_attribute!(WD, :is_complex_reflection_group, true)
    set_attribute!(WD, :complex_reflection_group_type, complex_reflection_group_type(W))

    return WD

end