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
export complex_reflection_group_type
export complex_reflection_group_model
export complex_reflection_group_dual
export complex_reflection_group_cartan_matrix

###########################################################################################
# Important remark: In OSCAR, matrices act by default from the right on vectors, so x*A. 
# For example the kernel of a matrix A is the space of all vectors such that x*A=0 
# (confusingly, this is usually referred to a s the *left kernel*....). The same holds for
# Magma and also CHEVIE.
# 
# This means that when we take matrices for the models from the literature in which 
# matrices act from the left (the usual convention for non-computer stuff), like the book
# by Lehrer & Taylor, we need to transpose these matrices! This is why down there in the
# code there are several transpose operation for the final generators.
###########################################################################################

function complex_reflection_group(G::ComplexReflectionGroupType, model::Symbol=:Magma)
    
    # this will be the list of matrix groups corresponding to the components of G
    component_groups = MatrixGroup[]

    # list of models of the components
    modellist = []

    for C in components(G)

        t = C.type[1]
        gens = nothing #this will be the list of generators
        refls = nothing #this will be the set of reflections

        # we now create the list "gens" of matrix generators for C in the selected model
        if isa(t, Int)

            # Exceptional groups

            # G4
            if t == 4

                if model == :LT
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
                
                    
                elseif model == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 ; -zeta_3 - 1 1])
                    s2 = matspace(K[1 zeta_3 + 1 ; 0 zeta_3])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 zeta_3])
                    s2 = matspace(K[1//3*(2*zeta_3 + 1) 1//3*(2*zeta_3 - 2) ; 1//3*(zeta_3 - 1) 1//3*(zeta_3 + 2)])

                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G5
            elseif t == 5

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 ; -2*zeta_3 - 2 1])
                    s2 = matspace(K[1 zeta_3 + 1 ; 0 zeta_3])

                    push!(gens, s1)
                    push!(gens, s2)


                elseif model == :CHEVIE
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 zeta_3])
                    s2 = matspace(K[1//3*(zeta_3 + 2) 1//3*(-2*zeta_3 + 2) ; 1//3*(-zeta_3 + 1) 1//3*(2*zeta_3 + 1)])
                    
                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G6
            elseif t == 6

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; -zeta_12^3 - 2*zeta_12^2 + 1 1])
                    s2 = matspace(K[1 zeta_12^2 ; 0 zeta_12^2 - 1])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//3*(-zeta_12^3 + 2*zeta_12) 1//3*(-2*zeta_12^3 + 4*zeta_12) ; 1//3*(-zeta_12^3 + 2*zeta_12) 1//3*(zeta_12^3 - 2*zeta_12)])
                    s2 = matspace(K[1 0 ; 0 zeta_12^2 - 1])
                    
                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G7
            elseif t == 7

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-zeta_12^2 - zeta_12 zeta_12^2 ; -2*zeta_12^2 - 2*zeta_12 zeta_12^2 + zeta_12])
                    s2 = matspace(K[1 0 ; -zeta_12^3 - zeta_12^2 + zeta_12 + 2 zeta_12^2 - 1])
                    s3 = matspace(K[zeta_12^2 - 1 0 ; zeta_12^3 + zeta_12^2 - zeta_12 - 2 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 -1])
                    s2 = matspace(K[1//2*(-zeta_12^3 + zeta_12^2 + zeta_12) 1//2*(-zeta_12^3 - zeta_12^2 + zeta_12) ; 1//2*(-zeta_12^3 + zeta_12^2 + zeta_12) 1//2*(zeta_12^3 + zeta_12^2 - zeta_12)])
                    s3 = matspace(K[1//2*(-zeta_12^3 + zeta_12^2 + zeta_12) 1//2*(-zeta_12^3 + zeta_12^2 + zeta_12) ; 1//2*(-zeta_12^3 - zeta_12^2 + zeta_12) 1//2*(zeta_12^3 + zeta_12^2 - zeta_12)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end   
                
            # G8
            elseif t == 8

                if model == :LT
                    # See Lehrer & Taylor (2009), page 87-88
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    r4 = matspace(K[1 0 ; 0 i])
                    r4prime = 1//2 * matspace(K[1+i -1+i ; -1+i 1+i])
                
                    push!(gens, transpose(r4))
                    push!(gens, transpose(r4prime))

                elseif model == :Magma
                    K,zeta_4 = cyclotomic_field(4)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_4 0 ; -zeta_4 1])
                    s2 = matspace(K[1 1 ; 0 zeta_4])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_4 = cyclotomic_field(4)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 zeta_4])
                    s2 = matspace(K[1//2*(zeta_4 + 1) 1//2*(zeta_4 - 1) ; 1//2*(zeta_4 - 1) 1//2*(zeta_4 + 1)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end   

            # G9
            elseif t == 9

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_8^3 + zeta_8^2 + zeta_8 1])
                    s2 = matspace(K[1 -zeta_8^3 ; 0 -zeta_8^2])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//2*(-zeta_8^3 + zeta_8) 1//2 ; 1 1//2*(zeta_8^3 - zeta_8)])
                    s2 = matspace(K[1 0 ; 0 zeta_8^2])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end  

            # G10
            elseif t == 10

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-zeta_12^3 0 ; -zeta_12^2 + zeta_12 1])
                    s2 = matspace(K[1 zeta_12^2 ; 0 zeta_12^2 - 1])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_12 = cyclotomic_field(12)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 zeta_12^2 - 1])
                    s2 = matspace(K[1//3*(zeta_12^3 - zeta_12^2 + zeta_12 + 2) 1//3*(zeta_12^3 + 2*zeta_12^2 - 2*zeta_12 - 1) ; 1//6*(zeta_12^3 + 2*zeta_12^2 - 2*zeta_12 - 1) 1//3*(2*zeta_12^3 + zeta_12^2 - zeta_12 + 1)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end  

            # G11
            elseif t == 11

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-zeta_24^6 0 ; -zeta_24^6 - zeta_24^3 - 1 1])
                    s2 = matspace(K[-zeta_24^6 - zeta_24^3 zeta_24^6 ; -2*zeta_24^6 - 2*zeta_24^3 - 1 zeta_24^6 + zeta_24^3])
                    s3 = matspace(K[-zeta_24^7 zeta_24^7 ; -zeta_24^7 - zeta_24^4 - zeta_24 zeta_24^7 + zeta_24^4])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//3*(2*zeta_24^7 - zeta_24^5 - zeta_24^3 - zeta_24) 1//3*(-zeta_24^7 + 2*zeta_24^5 + 2*zeta_24^3 - zeta_24) ; 1//6*(-zeta_24^7 - zeta_24^5 - zeta_24^3 + 2*zeta_24) 1//3*(-2*zeta_24^7 + zeta_24^5 + zeta_24^3 + zeta_24)])
                    s2 = matspace(K[1 0 ; 0 zeta_24^4 - 1])
                    s3 = matspace(K[1//3*(zeta_24^6 - zeta_24^4 + zeta_24^2 + 2) 1//3*(zeta_24^6 + 2*zeta_24^4 - 2*zeta_24^2 - 1) ; 1//6*(zeta_24^6 + 2*zeta_24^4 - 2*zeta_24^2 - 1) 1//3*(2*zeta_24^6 + zeta_24^4 - zeta_24^2 + 1)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G12
            elseif t == 12

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; -zeta_8^3 - zeta_8 + 1 1])
                    s2 = matspace(K[1 zeta_8^3 + zeta_8 + 1 ; 0 -1])
                    s3 = matspace(K[zeta_8^3 + zeta_8 - 1 -2 ; -zeta_8^3 - zeta_8 - 1 -zeta_8^3 - zeta_8 + 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//2 1//2*(zeta_8^3 + zeta_8 + 2) ; 1//4*(-zeta_8^3 - zeta_8 + 2) -1//2])
                    s2 = matspace(K[1//2 1//2*(-zeta_8^3 - zeta_8 + 2) ; 1//4*(zeta_8^3 + zeta_8 + 2) -1//2])
                    s3 = matspace(K[1 0 ; 0 -1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G13
            elseif t == 13

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_8^3 - zeta_8 - 2 1])
                    s2 = matspace(K[zeta_8^3 + zeta_8^2 -zeta_8^2 ; 2*zeta_8^3 + 2*zeta_8^2 - 1 -zeta_8^3 - zeta_8^2])
                    s3 = matspace(K[-zeta_8^3 + zeta_8 + 1 -1 ; -2*zeta_8^3 + 2*zeta_8 + 2 zeta_8^3 - zeta_8 - 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_8 = cyclotomic_field(8)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 -1])
                    s2 = matspace(K[1//2*(-zeta_8^3 + zeta_8) 1//2*(zeta_8^3 - zeta_8) ; 1//2*(zeta_8^3 - zeta_8) 1//2*(zeta_8^3 - zeta_8)])
                    s3 = matspace(K[1//2*(-zeta_8^3 + zeta_8) 1//2*(zeta_8^3 + zeta_8) ; 1//2*(-zeta_8^3 - zeta_8) 1//2*(zeta_8^3 - zeta_8)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G14
            elseif t == 14

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_24^4 - 1 0 ; -zeta_24^4 1])
                    s2 = matspace(K[1 -zeta_24^5 + 2*zeta_24^4 - zeta_24^3 + zeta_24 - 1 ; 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 -1])
                    s2 = matspace(K[1//2*(-zeta_24^7 + zeta_24^4 + zeta_24) 1//2*zeta_24^4 ; -1//2*zeta_24^4 1//2*(zeta_24^7 + zeta_24^4 - zeta_24)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end 

            # G15
            elseif t == 15

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_24^5 - zeta_24^3 - zeta_24 - 2 1])
                    s2 = matspace(K[-zeta_24^6 + zeta_24^4 + zeta_24^2 + zeta_24 -zeta_24 ; -2*zeta_24^6 + zeta_24^4 + 2*zeta_24^2 + 2*zeta_24 zeta_24^6 - zeta_24^2 - zeta_24])
                    s3 = matspace(K[-zeta_24^5 + zeta_24^3 + zeta_24 + 1 -1 ; -2*zeta_24^5 + 2*zeta_24^3 + 2*zeta_24 + 2 zeta_24^5 - zeta_24^3 - zeta_24 - 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_24 = cyclotomic_field(24)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//3*(2*zeta_24^7 - zeta_24^5 - zeta_24^3 - zeta_24) 1//3*(-zeta_24^7 + 2*zeta_24^5 + 2*zeta_24^3 - zeta_24) ; 1//6*(-zeta_24^7 - zeta_24^5 - zeta_24^3 + 2*zeta_24) 1//3*(-2*zeta_24^7 + zeta_24^5 + zeta_24^3 + zeta_24)])
                    s2 = matspace(K[1 0 ; 0 zeta_24^4 - 1])
                    s3 = matspace(K[1//3*(-zeta_24^6 + 2*zeta_24^2) 1//3*(2*zeta_24^6 - 4*zeta_24^2) ; 1//3*(zeta_24^6 - 2*zeta_24^2) 1//3*(zeta_24^6 - 2*zeta_24^2)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end 

            # G16
            elseif t == 16

                print(model)
                if model == :LT
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
                
                elseif model == :Magma
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_5 0 ; 1 1])
                    s2 = matspace(K[1 -zeta_5 ; 0 zeta_5])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 zeta_5])
                    s2 = matspace(K[1//5*(-zeta_5^3 - 2*zeta_5^2 + 2*zeta_5 + 1) 1//5*(2*zeta_5^3 - zeta_5^2 + zeta_5 - 2) ; 1//5*(-3*zeta_5^3 - zeta_5^2 + zeta_5 - 2) 1//5*(zeta_5^3 + 2*zeta_5^2 + 3*zeta_5 + 4)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end 

            # G17
            elseif t == 17

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_20 = cyclotomic_field(20)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; -zeta_20^6 + zeta_20^5 + zeta_20^4 - 2*zeta_20^2 + 1 1])
                    s2 = matspace(K[1 zeta_20^2 ; 0 zeta_20^4])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_20 = cyclotomic_field(20)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//5*(-3*zeta_20^7 + zeta_20^5 + zeta_20^3 + 2*zeta_20) 1//5*(zeta_20^7 - 2*zeta_20^5 + 3*zeta_20^3 - 4*zeta_20) ; 1//5*(zeta_20^7 - 2*zeta_20^5 + 3*zeta_20^3 - 4*zeta_20) 1//5*(3*zeta_20^7 - zeta_20^5 - zeta_20^3 - 2*zeta_20)])
                    s2 = matspace(K[1 0 ; 0 zeta_20^4])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end 

            # G18
            elseif t == 18

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_15^5 0 ; zeta_15^7 - 2*zeta_15^6 - zeta_15^3 + zeta_15^2 - zeta_15 - 1 1])
                    s2 = matspace(K[1 -zeta_15^7 + zeta_15^6 + zeta_15^3 - zeta_15^2 + 1 ; 0 zeta_15^3])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//5*(-3*zeta_15^7 + 4*zeta_15^5 - 2*zeta_15^4 - zeta_15 + 4) 1//5*(zeta_15^7 + 2*zeta_15^5 - zeta_15^4 - 3*zeta_15 + 2) ; 1//5*(zeta_15^7 - 3*zeta_15^5 + 4*zeta_15^4 + 2*zeta_15 - 3) 1//5*(3*zeta_15^7 + zeta_15^5 + 2*zeta_15^4 + zeta_15 + 1)])
                    s2 = matspace(K[1 0 ; 0 zeta_15^3])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end
                
            # G19
            elseif t == 19

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_60 = cyclotomic_field(60)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_60^12 0 ; zeta_60^12 + zeta_60^10 - 1 1])
                    s2 = matspace(K[1 -1 ; 0 zeta_60^10 - 1])
                    s3 = matspace(K[-zeta_60^13 - zeta_60^11 + zeta_60^3 + zeta_60 zeta_60^13 - zeta_60^3 ; -zeta_60^15 - 2*zeta_60^13 - 2*zeta_60^11 + zeta_60^5 + zeta_60^3 + zeta_60 zeta_60^13 + zeta_60^11 - zeta_60^3 - zeta_60])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_60 = cyclotomic_field(60)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//5*(zeta_60^15 - 3*zeta_60^11 + zeta_60^9 + 2*zeta_60^3 + 3*zeta_60) 1//5*(-2*zeta_60^15 + zeta_60^11 + 3*zeta_60^9 - 4*zeta_60^3 - zeta_60) ; 1//5*(-2*zeta_60^15 + zeta_60^11 + 3*zeta_60^9 - 4*zeta_60^3 - zeta_60) 1//5*(-zeta_60^15 + 3*zeta_60^11 - zeta_60^9 - 2*zeta_60^3 - 3*zeta_60)])
                    s2 = matspace(K[1//5*(-zeta_60^14 - 3*zeta_60^12 + 2*zeta_60^10 + zeta_60^8 + zeta_60^6 + 2*zeta_60^4 + 2*zeta_60^2 - 1) 1//5*(2*zeta_60^14 + zeta_60^12 + zeta_60^10 - 2*zeta_60^8 - 2*zeta_60^6 - 4*zeta_60^4 + zeta_60^2 + 2) ; 1//5*(-3*zeta_60^14 + zeta_60^12 + zeta_60^10 + 3*zeta_60^8 + 3*zeta_60^6 + zeta_60^4 - 4*zeta_60^2 - 3) 1//5*(zeta_60^14 + 3*zeta_60^12 + 3*zeta_60^10 - zeta_60^8 - zeta_60^6 - 2*zeta_60^4 - 2*zeta_60^2 + 1)])
                    s3 = matspace(K[1 0 ; 0 zeta_60^12])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end
                
            # G20
            elseif t == 20

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_15^5 0 ; -2*zeta_15^5 + zeta_15^4 + zeta_15 - 2 1])
                    s2 = matspace(K[1 zeta_15^5 + 1 ; 0 zeta_15^5])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 ; 0 zeta_15^5])
                    s2 = matspace(K[1//3*(-2*zeta_15^7 + 2*zeta_15^5 - zeta_15^4 + 2*zeta_15^3 - 2*zeta_15^2 - zeta_15 + 3) 1//3*(4*zeta_15^7 - zeta_15^5 + 2*zeta_15^4 - 4*zeta_15^3 + 4*zeta_15^2 + 2*zeta_15 - 3) ; 1//15*(4*zeta_15^7 - zeta_15^5 + 2*zeta_15^4 - 4*zeta_15^3 + 4*zeta_15^2 + 2*zeta_15 - 3) 1//3*(2*zeta_15^7 + zeta_15^5 + zeta_15^4 - 2*zeta_15^3 + 2*zeta_15^2 + zeta_15)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end

            # G21
            elseif t == 21

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_60 = cyclotomic_field(60)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_60^11 - 2*zeta_60^10 + zeta_60^9 - zeta_60 + 1 1])
                    s2 = matspace(K[1 zeta_60^10 ; 0 zeta_60^10 - 1])

                    push!(gens, s1)
                    push!(gens, s2)

                elseif model == :CHEVIE
                    K,zeta_60 = cyclotomic_field(60)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//3*(-2*zeta_60^15 - 2*zeta_60^13 + zeta_60^11 + zeta_60^9 + 2*zeta_60^7 + 2*zeta_60^5 - zeta_60) 1//3*(zeta_60^15 - 2*zeta_60^13 + zeta_60^11 + zeta_60^9 + 2*zeta_60^7 - 4*zeta_60^5 - zeta_60) ; 1//15*(zeta_60^15 - 2*zeta_60^13 + zeta_60^11 + zeta_60^9 + 2*zeta_60^7 - 4*zeta_60^5 - zeta_60) 1//3*(2*zeta_60^15 + 2*zeta_60^13 - zeta_60^11 - zeta_60^9 - 2*zeta_60^7 - 2*zeta_60^5 + zeta_60)])
                    s2 = matspace(K[1 0 ; 0 zeta_60^10 - 1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)

                end
            
            # G22
            elseif t == 22

                if model == :LT
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

                elseif model == :Magma
                    K,zeta_20 = cyclotomic_field(20)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 ; zeta_20^6 + zeta_20^5 - zeta_20^4 1])
                    s2 = matspace(K[1 zeta_20^6 - zeta_20^5 - zeta_20^4 ; 0 -1])
                    s3 = matspace(K[zeta_20^7 - zeta_20^6 - zeta_20^5 + zeta_20^4 + zeta_20^3 - 1 -2*zeta_20^6 + 2*zeta_20^4 ; -zeta_20^7 + 2*zeta_20^5 - zeta_20^3 + 1 -zeta_20^7 + zeta_20^6 + zeta_20^5 - zeta_20^4 - zeta_20^3 + 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_20 = cyclotomic_field(20)
                    matspace = matrix_space(K, 2, 2)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//5*(zeta_20^7 - 2*zeta_20^5 + 3*zeta_20^3 - 4*zeta_20) 1//5*(3*zeta_20^7 - zeta_20^5 - zeta_20^3 - 2*zeta_20) ; 1//5*(3*zeta_20^7 - zeta_20^5 - zeta_20^3 - 2*zeta_20) 1//5*(-zeta_20^7 + 2*zeta_20^5 - 3*zeta_20^3 + 4*zeta_20)])
                    s2 = matspace(K[0 -zeta_20 ; zeta_20^7 - zeta_20^5 + zeta_20^3 - zeta_20 0])
                    s3 = matspace(K[0 zeta_20^7 - zeta_20^5 + zeta_20^3 - zeta_20 ; -zeta_20 0])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end

            # G23
            elseif t == 23
                
                if model == :LT
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


                elseif model == :Magma
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 ; -zeta_5^3 - zeta_5^2 1 0 ; 0 0 1])
                    s2 = matspace(K[1 -zeta_5^3 - zeta_5^2 0 ; 0 -1 0 ; 0 1 1])
                    s3 = matspace(K[1 0 0 ; 0 1 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 -zeta_5^3 - zeta_5^2 0 ; 0 1 0 ; 0 0 1])
                    s2 = matspace(K[1 0 0 ; -zeta_5^3 - zeta_5^2 -1 1 ; 0 0 1])
                    s3 = matspace(K[1 0 0 ; 0 1 0 ; 0 1 -1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end
                
           
            # G24
            elseif t == 24
                
                if model == :LT
                    # See Lehrer & Taylor (2009), page 108
                    # This group is realized by the line system J_3^{(4)}
                    R,x = polynomial_ring(QQ)
                    K,i = number_field(x^2+1, "i")
                    R,x = polynomial_ring(K)
                    K,sqrt7 = number_field(x^2-7, "√7")
                    i = K(i)
                    V = vector_space(K,3)
                    matspace = matrix_space(K, 3, 3)

                    # The the line system J_3^{(4)}
                    lambda = -1//2*(1 + i*sqrt7)

                    lines = [1//2*V([lambda^2,lambda^2,0]), 1//2*V([lambda^2,-lambda^2,0])]
                    lines = vcat(lines, [1//2*V([lambda^2,0,lambda^2]), 1//2*V([lambda^2,0,-lambda^2])])
                    lines = vcat(lines, [1//2*V([0,lambda^2,lambda^2]), 1//2*V([0,lambda^2,-lambda^2])])
                    lines = vcat(lines, [V([lambda,0,0]), V([0,lambda,0]), V([0,0,lambda])])

                    lines = vcat(lines, [1//2*V([lambda,lambda,K(2)]), 1//2*V([lambda,-lambda,K(2)]), 1//2*V([-lambda,lambda,K(2)]), 1//2*V([-lambda,-lambda,K(2)])])
                    lines = vcat(lines, [1//2*V([lambda,K(2),lambda]), 1//2*V([lambda,K(2),-lambda]), 1//2*V([-lambda,K(2),lambda]), 1//2*V([-lambda,K(2),-lambda])]) 
                    lines = vcat(lines, [1//2*V([K(2),lambda,lambda]), 1//2*V([K(2),lambda,-lambda]), 1//2*V([K(2),-lambda,lambda]), 1//2*V([K(2),-lambda,- lambda])])                 

                    refls = Set([unitary_reflection(l) for l in lines])

                    # I guessed the roots of generating reflections from the proof of 
                    # Lemma 7.34 (page 129)

                    x = 1//2*V([lamba,lamba,K(2)])  

                    # Now, create the corresponding generators
                    #gens = elem_type(matspace)[]

                    

                elseif model == :Magma
                    K,zeta_7 = cyclotomic_field(7)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 ; -1 1 0 ; -zeta_7^4 - zeta_7^2 - zeta_7 - 1 0 1])
                    s2 = matspace(K[1 -1 0 ; 0 -1 0 ; 0 -zeta_7^4 - zeta_7^2 - zeta_7 1])
                    s3 = matspace(K[1 0 zeta_7^4 + zeta_7^2 + zeta_7 ; 0 1 zeta_7^4 + zeta_7^2 + zeta_7 + 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_7 = cyclotomic_field(7)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//2 1//14*(6*zeta_7^4 + 6*zeta_7^2 + 6*zeta_7 + 3) 0 ; 1//2*(-2*zeta_7^4 - 2*zeta_7^2 - 2*zeta_7 - 1) -1//2 0 ; 0 0 1])
                    s2 = matspace(K[1//2 1//14*(-6*zeta_7^4 - 6*zeta_7^2 - 6*zeta_7 - 3) 0 ; 1//2*(2*zeta_7^4 + 2*zeta_7^2 + 2*zeta_7 + 1) -1//2 0 ; 0 0 1])
                    s3 = matspace(K[0 1//7*(-zeta_7^4 - zeta_7^2 - zeta_7 - 4) 1//4*(-zeta_7^4 - zeta_7^2 - zeta_7 - 1) ; 1//3*(zeta_7^4 + zeta_7^2 + zeta_7 - 3) 1//3 1//12*(-3*zeta_7^4 - 3*zeta_7^2 - 3*zeta_7 - 5) ; 1//3*(2*zeta_7^4 + 2*zeta_7^2 + 2*zeta_7) 1//21*(6*zeta_7^4 + 6*zeta_7^2 + 6*zeta_7 - 4) 2//3])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end

            # G25
            elseif t == 25
                
                if model == :LT
                    nothing

                elseif model == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 0 ; zeta_3 + 1 1 0 ; 0 0 1])
                    s2 = matspace(K[1 -zeta_3 - 1 0 ; 0 zeta_3 0 ; 0 -zeta_3 - 1 1])
                    s3 = matspace(K[1 0 0 ; 0 1 zeta_3 + 1 ; 0 0 zeta_3])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 0 ; 0 1 0 ; 0 0 zeta_3])
                    s2 = matspace(K[1//3*(zeta_3 + 2) 1//3*(zeta_3 - 1) 1//3*(zeta_3 - 1) ; 1//3*(zeta_3 - 1) 1//3*(zeta_3 + 2) 1//3*(zeta_3 - 1) ; 1//3*(zeta_3 - 1) 1//3*(zeta_3 - 1) 1//3*(zeta_3 + 2)])
                    s3 = matspace(K[1 0 0 ; 0 zeta_3 0 ; 0 0 1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end

            # G26
            elseif t == 26
                
                if model == :LT
                    nothing

                elseif model == :Magma
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[zeta_3 0 0 ; zeta_3 + 1 1 0 ; 0 0 1])
                    s2 = matspace(K[1 -zeta_3 - 1 0 ; 0 zeta_3 0 ; 0 -zeta_3 + 1 1])
                    s3 = matspace(K[1 0 0 ; 0 1 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 0 ; 0 0 1 ; 0 1 0])
                    s2 = matspace(K[1 0 0 ; 0 1 0 ; 0 0 zeta_3])
                    s3 = matspace(K[1//3*(zeta_3 + 2) 1//3*(zeta_3 - 1) 1//3*(zeta_3 - 1) ; 1//3*(zeta_3 - 1) 1//3*(zeta_3 + 2) 1//3*(zeta_3 - 1) ; 1//3*(zeta_3 - 1) 1//3*(zeta_3 - 1) 1//3*(zeta_3 + 2)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end

            # G27
            elseif t == 27
                
                if model == :LT
                    nothing

                elseif model == :Magma
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 ; -1 1 0 ; zeta_15^7 + zeta_15^4 - zeta_15^3 + zeta_15^2 + zeta_15 - 1 0 1])
                    s2 = matspace(K[1 -1 0 ; 0 -1 0 ; 0 zeta_15^7 + zeta_15^4 - zeta_15^3 + zeta_15^2 + zeta_15 1])
                    s3 = matspace(K[1 0 -zeta_15^4 - zeta_15 ; 0 1 -zeta_15^4 - zeta_15 + 1 ; 0 0 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                elseif model == :CHEVIE
                    K,zeta_15 = cyclotomic_field(15)
                    matspace = matrix_space(K, 3, 3)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1//2 1//10*(-2*zeta_15^7 + 2*zeta_15^6 + 2*zeta_15^5 - 6*zeta_15^4 + 4*zeta_15^3 + 2*zeta_15^2 + 2*zeta_15 + 3) 1//10*(-7*zeta_15^7 + 2*zeta_15^6 + 2*zeta_15^5 - zeta_15^4 - zeta_15^3 - 3*zeta_15^2 - 3*zeta_15 + 3) ; 1//10*(-8*zeta_15^7 + 6*zeta_15^6 + 2*zeta_15^5 - 6*zeta_15^4 + 2*zeta_15^3 - 4*zeta_15^2 + 2*zeta_15 + 5) 1//10*(2*zeta_15^7 - 2*zeta_15^3 + 2*zeta_15^2 - 1) 1//10*(7*zeta_15^7 + 5*zeta_15^4 - 7*zeta_15^3 + 7*zeta_15^2 + 5*zeta_15 - 1) ; 1//15*(2*zeta_15^7 - 4*zeta_15^6 + 2*zeta_15^5 - zeta_15^4 + 2*zeta_15^3 + 6*zeta_15^2 - 3*zeta_15) 1//15*(2*zeta_15^7 - 5*zeta_15^4 - 2*zeta_15^3 + 2*zeta_15^2 - 5*zeta_15 + 4) 1//5*(-zeta_15^7 + zeta_15^3 - zeta_15^2 + 3)])
                    s2 = matspace(K[0 1//5*(2*zeta_15^7 - zeta_15^5 + 3*zeta_15^4 - zeta_15 - 1) 1//10*(-11*zeta_15^7 + 10*zeta_15^6 + 3*zeta_15^5 - 4*zeta_15^4 + 5*zeta_15^3 - 5*zeta_15^2 - 2*zeta_15 + 8) ; 1//5*(3*zeta_15^7 - 4*zeta_15^6 - zeta_15^5 + 3*zeta_15^4 - 3*zeta_15^3 + zeta_15^2 - zeta_15 - 3) 1//5*(zeta_15^7 - zeta_15^3 + zeta_15^2 + 2) 1//10*(-3*zeta_15^7 + 5*zeta_15^5 + 3*zeta_15^3 - 3*zeta_15^2 + 4) ; 1//15*(zeta_15^7 - 8*zeta_15^6 + 3*zeta_15^5 - 4*zeta_15^4 - zeta_15^3 + 7*zeta_15^2 - 2*zeta_15 - 1) 1//15*(-3*zeta_15^7 - 5*zeta_15^5 + 3*zeta_15^3 - 3*zeta_15^2 - 1) 1//5*(-zeta_15^7 + zeta_15^3 - zeta_15^2 + 3)])
                    s3 = matspace(K[-1 0 0 ; 0 1 0 ; 0 0 1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)

                end

            # G28
            elseif t == 28
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K = QQ
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 1 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 ; 1 -1 1 0 ; 0 0 1 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 2 -1 1 ; 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 1 -1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)

                end

            # G29
            elseif t == 29
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K,zeta_4 = cyclotomic_field(4)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 -1])
                    s2 = matspace(K[1//2 1//2 1//2*zeta_4 1//2*zeta_4 ; 1//2 1//2 -1//2*zeta_4 -1//2*zeta_4 ; -1//2*zeta_4 1//2*zeta_4 1//2 -1//2 ; -1//2*zeta_4 1//2*zeta_4 -1//2 1//2])
                    s3 = matspace(K[0 1 0 0 ; 1 0 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)

                end

            # G30
            elseif t == 30
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K,zeta_5 = cyclotomic_field(5)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 -zeta_5^3 - zeta_5^2 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 ; -zeta_5^3 - zeta_5^2 -1 1 0 ; 0 0 1 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 1 -1 1 ; 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 1 -1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)

                end

            # G31
            elseif t == 31
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K,zeta_4 = cyclotomic_field(4)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s2 = matspace(K[0 zeta_4 0 0 ; -zeta_4 0 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s3 = matspace(K[0 1 0 0 ; 1 0 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s4 = matspace(K[1//2 -1//2 -1//2 -1//2 ; -1//2 1//2 -1//2 -1//2 ; -1//2 -1//2 1//2 -1//2 ; -1//2 -1//2 -1//2 1//2])
                    s5 = matspace(K[1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)

                end

            # G32
            elseif t == 32
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 4, 4)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 zeta_3 0 ; 0 0 0 1])
                    s2 = matspace(K[1//3*(zeta_3 + 2) 1//3*(zeta_3 - 1) 1//3*(zeta_3 - 1) 0 ; 1//3*(zeta_3 - 1) 1//3*(zeta_3 + 2) 1//3*(zeta_3 - 1) 0 ; 1//3*(zeta_3 - 1) 1//3*(zeta_3 - 1) 1//3*(zeta_3 + 2) 0 ; 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 ; 0 zeta_3 0 0 ; 0 0 1 0 ; 0 0 0 1])
                    s4 = matspace(K[1//3*(zeta_3 + 2) 1//3*(-zeta_3 + 1) 0 1//3*(-zeta_3 + 1) ; 1//3*(-zeta_3 + 1) 1//3*(zeta_3 + 2) 0 1//3*(zeta_3 - 1) ; 0 0 1 0 ; 1//3*(-zeta_3 + 1) 1//3*(zeta_3 - 1) 0 1//3*(zeta_3 + 2)])
                                        
                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)

                end

            # G3
            elseif t == 33
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 5, 5)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[0 1 0 0 0 ; 1 0 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ; 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 ; 0 0 0 0 1 ; 0 0 1 0 0 ; 0 0 0 1 0 ; 0 1 0 0 0])
                    s3 = matspace(K[1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 0 1 ; 0 0 0 1 0])
                    s4 = matspace(K[1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 0 -zeta_3 - 1 ; 0 0 0 zeta_3 0])
                    s5 = matspace(K[2//3 -1//3 1//3 -1//3 -1//3 ; -1//3 2//3 1//3 -1//3 -1//3 ; 2//3 2//3 1//3 2//3 2//3 ; -1//3 -1//3 1//3 2//3 -1//3 ; -1//3 -1//3 1//3 -1//3 2//3])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)

                end

            # G34
            elseif t == 34
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K,zeta_3 = cyclotomic_field(3)
                    matspace = matrix_space(K, 6, 6)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[0 1 0 0 0 0 ; 1 0 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 0 ; 0 0 0 0 0 1 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 1 0 0 0 0])
                    s3 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 0 1 ; 0 0 0 0 1 0])
                    s4 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 0 -zeta_3 - 1 ; 0 0 0 0 zeta_3 0])
                    s5 = matspace(K[2//3 -1//3 1//3*(zeta_3 + 1) -1//3*zeta_3 -1//3 -1//3 ; -1//3 2//3 1//3*(zeta_3 + 1) -1//3*zeta_3 -1//3 -1//3 ; -1//3*zeta_3 -1//3*zeta_3 2//3 1//3*(zeta_3 + 1) -1//3*zeta_3 -1//3*zeta_3 ; 1//3*(zeta_3 + 1) 1//3*(zeta_3 + 1) -1//3*zeta_3 2//3 1//3*(zeta_3 + 1) 1//3*(zeta_3 + 1) ; -1//3 -1//3 1//3*(zeta_3 + 1) -1//3*zeta_3 2//3 -1//3 ; -1//3 -1//3 1//3*(zeta_3 + 1) -1//3*zeta_3 -1//3 2//3])
                    s6 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 0 zeta_3 0 0 ; 0 0 -zeta_3 - 1 0 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)
                    push!(gens, s6)

                end

            # G35
            elseif t == 35
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K = QQ
                    matspace = matrix_space(K, 6, 6)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 1 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 0 ; 0 -1 0 1 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 1 0 -1 1 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 1 1 -1 1 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
                    s5 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 1 -1 1 ; 0 0 0 0 0 1])
                    s6 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 1 -1])

                    push!(gens, s1)
                    push!(gens, s2)
                    push!(gens, s3)
                    push!(gens, s4)
                    push!(gens, s5)
                    push!(gens, s6)

                end

            # G36
            elseif t == 36
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K = QQ
                    matspace = matrix_space(K, 7, 7)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 1 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 0 0 ; 0 -1 0 1 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 1 0 -1 1 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 1 1 -1 1 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s5 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 1 -1 1 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1])
                    s6 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 1 -1 1 ; 0 0 0 0 0 0 1])
                    s7 = matspace(K[1 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 0 1 0 0 0 ; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 1 -1])

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
                
                if model == :LT
                    nothing

                elseif model == :Magma
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

                elseif model == :CHEVIE
                    K = QQ
                    matspace = matrix_space(K, 8, 8)
                    gens = elem_type(matspace)[]

                    s1 = matspace(K[-1 0 1 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s2 = matspace(K[1 0 0 0 0 0 0 0 ; 0 -1 0 1 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s3 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 1 0 -1 1 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s4 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 1 1 -1 1 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s5 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 1 -1 1 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s6 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 1 -1 1 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 1])
                    s7 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 -1 1 ; 0 0 0 0 0 0 0 1])
                    s8 = matspace(K[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 1 -1])

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

                # The generators correspond to the transpositions r_i=(i,i+1)
                # but note that we consider the irreducible representation of 
                # dimension n-1, so it's not simply the transposition matrix.

                # In Magma, hephardTodd(1,1,n) gives the permutation
                # representation, which is not irreducible...So we will deviate here from
                # Magma.

                # In LT there are no explicit matrices given.

                if model === :CHEVIE || model === :LT || model === :Magma
                    
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

                if model === :LT
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

                elseif model === :Magma
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
                    t[n,n] = z #LT takes t[1,1]=z here

                    #the matrix s
                    if n > 1
                        s = t^-1*transp[n-1]*t #LT takes transp[1] here
                    end

                    #now, create the list of generators
                    #LT orders the transpositions to the end
                    if m == p
                        gens = vcat(transp, [s])
                    elseif p == 1
                        gens = vcat(transp, [t])
                    else
                        gens = vcat(transp, [s], [t^p])
                    end

                elseif model === :CHEVIE
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
                        gens = vcat([s1prime], transp)
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
        set_attribute!(matgrp, :complex_reflection_group_model, [model])
        set_attribute!(matgrp, :is_irreducible, true)

        if refls != nothing
            set_attribute!(matgrp, :complex_reflections, refls)
        end

        # add to list
        push!(component_groups, matgrp)
        push!(modellist, model)
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

        return matgrp
    end

end

# Convenience constructors
complex_reflection_group(i::Int, model::Symbol=:Magma) = complex_reflection_group(ComplexReflectionGroupType(i), model)

complex_reflection_group(m::Int, p::Int, n::Int, model::Symbol=:Magma) = complex_reflection_group(ComplexReflectionGroupType(m,p,n), model)

complex_reflection_group(X::Vector, model::Symbol=:Magma) = complex_reflection_group(ComplexReflectionGroupType(X), model)

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

function complex_reflection_group_cartan_matrix(W::MatrixGroup)

    if !is_complex_reflection_group(W)
        throw(ArgumentError("Group is not a complex reflection group"))
    end

    # We collect roots and coroots of the generators of W
    gen_roots = []
    gen_coroots = []

    # If reflections are assigned already, we take these roots and coroots
    if has_attribute(W, :complex_reflections)
        refls = collect(get_attribute(W, :complex_reflections))
        for g in gens(W)
            i = findfirst(w->matrix(w)==matrix(g), collect(refls))
            push!(gen_roots, root(refls[i]))
            push!(gen_coroots, coroot(refls[i]))
        end
    else
    # Otherwise, we compute roots and coroots
        for g in gens(W)
            b,g_data = is_complex_reflection_with_data(g)
            push!(gen_roots, root(g_data))
            push!(gen_coroots, coroot(g_data))
        end
    end

    K = base_ring(W)
    n = length(gen_roots)
    C = matrix(K,n,n,[ canonical_pairing(gen_coroots[j], gen_roots[i]) for i=1:n for j=1:n ]) 

    return C
end