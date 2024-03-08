
###########################################################################################
# Models from Magma
###########################################################################################
function complex_reflection_group_Magma(n::Int)

  if n == 4
    K,z_3 = cyclotomic_field(3)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_3 0 ; -z_3 - 1 1])
    s2 = matspace(K[1 z_3 + 1 ; 0 z_3])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 5
    K,z_3 = cyclotomic_field(3)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_3 0 ; -2*z_3 - 2 1])
    s2 = matspace(K[1 z_3 + 1 ; 0 z_3])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 6
    K,z_12 = cyclotomic_field(12)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; -z_12^3 - 2*z_12^2 + 1 1])
    s2 = matspace(K[1 z_12^2 ; 0 z_12^2 - 1])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 7
    K,z_12 = cyclotomic_field(12)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-z_12^2 - z_12 z_12^2 ; -2*z_12^2 - 2*z_12 z_12^2 + z_12])
    s2 = matspace(K[1 0 ; -z_12^3 - z_12^2 + z_12 + 2 z_12^2 - 1])
    s3 = matspace(K[z_12^2 - 1 0 ; z_12^3 + z_12^2 - z_12 - 2 1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 8
    K,z_4 = cyclotomic_field(4)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_4 0 ; -z_4 1])
    s2 = matspace(K[1 1 ; 0 z_4])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 9
    K,z_8 = cyclotomic_field(8)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; z_8^3 + z_8^2 + z_8 1])
    s2 = matspace(K[1 -z_8^3 ; 0 -z_8^2])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 10
    K,z_12 = cyclotomic_field(12)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-z_12^3 0 ; -z_12^2 + z_12 1])
    s2 = matspace(K[1 z_12^2 ; 0 z_12^2 - 1])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 11
    K,z_24 = cyclotomic_field(24)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-z_24^6 0 ; -z_24^6 - z_24^3 - 1 1])
    s2 = matspace(K[-z_24^6 - z_24^3 z_24^6 ; -2*z_24^6 - 2*z_24^3 - 1 z_24^6 + z_24^3])
    s3 = matspace(K[-z_24^7 z_24^7 ; -z_24^7 - z_24^4 - z_24 z_24^7 + z_24^4])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    
    W = matrix_group(gens)

  elseif n == 12
    K,z_8 = cyclotomic_field(8)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; -z_8^3 - z_8 + 1 1])
    s2 = matspace(K[1 z_8^3 + z_8 + 1 ; 0 -1])
    s3 = matspace(K[z_8^3 + z_8 - 1 -2 ; -z_8^3 - z_8 - 1 -z_8^3 - z_8 + 1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 13
    K,z_8 = cyclotomic_field(8)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; z_8^3 - z_8 - 2 1])
    s2 = matspace(K[z_8^3 + z_8^2 -z_8^2 ; 2*z_8^3 + 2*z_8^2 - 1 -z_8^3 - z_8^2])
    s3 = matspace(K[-z_8^3 + z_8 + 1 -1 ; -2*z_8^3 + 2*z_8 + 2 z_8^3 - z_8 - 1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 14
    K,z_24 = cyclotomic_field(24)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_24^4 - 1 0 ; -z_24^4 1])
    s2 = matspace(K[1 -z_24^5 + 2*z_24^4 - z_24^3 + z_24 - 1 ; 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    
    W = matrix_group(gens)

  elseif n == 15
    K,z_24 = cyclotomic_field(24)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; z_24^5 - z_24^3 - z_24 - 2 1])
    s2 = matspace(K[-z_24^6 + z_24^4 + z_24^2 + z_24 -z_24 ; -2*z_24^6 + z_24^4 + 2*z_24^2 + 2*z_24 z_24^6 - z_24^2 - z_24])
    s3 = matspace(K[-z_24^5 + z_24^3 + z_24 + 1 -1 ; -2*z_24^5 + 2*z_24^3 + 2*z_24 + 2 z_24^5 - z_24^3 - z_24 - 1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    
    W = matrix_group(gens)

  elseif n == 16
    K,z_5 = cyclotomic_field(5)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_5 0 ; 1 1])
    s2 = matspace(K[1 -z_5 ; 0 z_5])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 17
    K,z_20 = cyclotomic_field(20)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; -z_20^6 + z_20^5 + z_20^4 - 2*z_20^2 + 1 1])
    s2 = matspace(K[1 z_20^2 ; 0 z_20^4])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 18
    K,z_15 = cyclotomic_field(15)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_15^5 0 ; z_15^7 - 2*z_15^6 - z_15^3 + z_15^2 - z_15 - 1 1])
    s2 = matspace(K[1 -z_15^7 + z_15^6 + z_15^3 - z_15^2 + 1 ; 0 z_15^3])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 19
    K,z_60 = cyclotomic_field(60)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_60^12 0 ; z_60^12 + z_60^10 - 1 1])
    s2 = matspace(K[1 -1 ; 0 z_60^10 - 1])
    s3 = matspace(K[-z_60^13 - z_60^11 + z_60^3 + z_60 z_60^13 - z_60^3 ; -z_60^15 - 2*z_60^13 - 2*z_60^11 + z_60^5 + z_60^3 + z_60 z_60^13 + z_60^11 - z_60^3 - z_60])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    
    W = matrix_group(gens)

  elseif n == 20
    K,z_15 = cyclotomic_field(15)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_15^5 0 ; -2*z_15^5 + z_15^4 + z_15 - 2 1])
    s2 = matspace(K[1 z_15^5 + 1 ; 0 z_15^5])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 21
    K,z_60 = cyclotomic_field(60)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; z_60^11 - 2*z_60^10 + z_60^9 - z_60 + 1 1])
    s2 = matspace(K[1 z_60^10 ; 0 z_60^10 - 1])

    push!(gens, s1)
    push!(gens, s2)

    W = matrix_group(gens)

  elseif n == 22
    K,z_20 = cyclotomic_field(20)
    matspace = matrix_space(K, 2, 2)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 ; z_20^6 + z_20^5 - z_20^4 1])
    s2 = matspace(K[1 z_20^6 - z_20^5 - z_20^4 ; 0 -1])
    s3 = matspace(K[z_20^7 - z_20^6 - z_20^5 + z_20^4 + z_20^3 - 1 -2*z_20^6 + 2*z_20^4 ; -z_20^7 + 2*z_20^5 - z_20^3 + 1 -z_20^7 + z_20^6 + z_20^5 - z_20^4 - z_20^3 + 1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 23
    K,z_5 = cyclotomic_field(5)
    matspace = matrix_space(K, 3, 3)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 ; -z_5^3 - z_5^2 1 0 ; 0 0 1])
    s2 = matspace(K[1 -z_5^3 - z_5^2 0 ; 0 -1 0 ; 0 1 1])
    s3 = matspace(K[1 0 0 ; 0 1 1 ; 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 24
    K,z_7 = cyclotomic_field(7)
    matspace = matrix_space(K, 3, 3)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 ; -1 1 0 ; -z_7^4 - z_7^2 - z_7 - 1 0 1])
    s2 = matspace(K[1 -1 0 ; 0 -1 0 ; 0 -z_7^4 - z_7^2 - z_7 1])
    s3 = matspace(K[1 0 z_7^4 + z_7^2 + z_7 ; 0 1 z_7^4 + z_7^2 + z_7 + 1 ; 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 25
    K,z_3 = cyclotomic_field(3)
    matspace = matrix_space(K, 3, 3)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_3 0 0 ; z_3 + 1 1 0 ; 0 0 1])
    s2 = matspace(K[1 -z_3 - 1 0 ; 0 z_3 0 ; 0 -z_3 - 1 1])
    s3 = matspace(K[1 0 0 ; 0 1 z_3 + 1 ; 0 0 z_3])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 26
    K,z_3 = cyclotomic_field(3)
    matspace = matrix_space(K, 3, 3)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_3 0 0 ; z_3 + 1 1 0 ; 0 0 1])
    s2 = matspace(K[1 -z_3 - 1 0 ; 0 z_3 0 ; 0 -z_3 + 1 1])
    s3 = matspace(K[1 0 0 ; 0 1 1 ; 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)

    W = matrix_group(gens)

  elseif n == 27
    K,z_15 = cyclotomic_field(15)
    matspace = matrix_space(K, 3, 3)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 ; -1 1 0 ; z_15^7 + z_15^4 - z_15^3 + z_15^2 + z_15 - 1 0 1])
    s2 = matspace(K[1 -1 0 ; 0 -1 0 ; 0 z_15^7 + z_15^4 - z_15^3 + z_15^2 + z_15 1])
    s3 = matspace(K[1 0 -z_15^4 - z_15 ; 0 1 -z_15^4 - z_15 + 1 ; 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    
    W = matrix_group(gens)

  elseif n == 28
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

    W = matrix_group(gens)

  elseif n == 29
    K,z_4 = cyclotomic_field(4)
    matspace = matrix_space(K, 4, 4)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 0 ; 1 1 0 0 ; z_4 0 1 0 ; 0 0 0 1])
    s2 = matspace(K[1 1 0 0 ; 0 -1 0 0 ; 0 -z_4 + 1 1 0 ; 0 0 0 1])
    s3 = matspace(K[1 0 -z_4 0 ; 0 1 z_4 + 1 0 ; 0 0 -1 0 ; 0 0 1 1])
    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 1 ; 0 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    push!(gens, s4)

    W = matrix_group(gens)

  elseif n == 30
    K,z_5 = cyclotomic_field(5)
    matspace = matrix_space(K, 4, 4)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 0 ; -z_5^3 - z_5^2 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
    s2 = matspace(K[1 -z_5^3 - z_5^2 0 0 ; 0 -1 0 0 ; 0 1 1 0 ; 0 0 0 1])
    s3 = matspace(K[1 0 0 0 ; 0 1 1 0 ; 0 0 -1 0 ; 0 0 1 1])
    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 1 ; 0 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    push!(gens, s4)

    W = matrix_group(gens)

  elseif n == 31
    K,z_4 = cyclotomic_field(4)
    matspace = matrix_space(K, 4, 4)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 0 ; 1 1 0 0 ; z_4 - 1 0 1 0 ; 0 0 0 1])
    s2 = matspace(K[1 1 0 0 ; 0 -1 0 0 ; 0 -z_4 1 0 ; 0 0 0 1])
    s3 = matspace(K[1 0 -z_4 - 1 0 ; 0 1 z_4 0 ; 0 0 -1 0 ; 0 0 1 1])
    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 1 ; 0 0 0 -1])
    s5 = matspace(K[-z_4 0 z_4 - 1 0 ; 0 1 0 0 ; -z_4 - 1 0 z_4 0 ; z_4 0 1 1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    push!(gens, s4)
    push!(gens, s5)

    W = matrix_group(gens)

  elseif n == 32
    K,z_3 = cyclotomic_field(3)
    matspace = matrix_space(K, 4, 4)
    gens = elem_type(matspace)[]

    s1 = matspace(K[z_3 0 0 0 ; -z_3 - 1 1 0 0 ; 0 0 1 0 ; 0 0 0 1])
    s2 = matspace(K[1 z_3 + 1 0 0 ; 0 z_3 0 0 ; 0 z_3 + 1 1 0 ; 0 0 0 1])
    s3 = matspace(K[1 0 0 0 ; 0 1 -z_3 - 1 0 ; 0 0 z_3 0 ; 0 0 -z_3 - 1 1])
    s4 = matspace(K[1 0 0 0 ; 0 1 0 0 ; 0 0 1 z_3 + 1 ; 0 0 0 z_3])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    push!(gens, s4)

    W = matrix_group(gens)

  elseif n == 33
    K,z_3 = cyclotomic_field(3)
    matspace = matrix_space(K, 5, 5)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 0 0 ; 1 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ; 0 0 0 0 1])
    s2 = matspace(K[1 1 0 0 0 ; 0 -1 0 0 0 ; 0 1 1 0 0 ; 0 -z_3 - 1 0 1 0 ; 0 0 0 0 1])
    s3 = matspace(K[1 0 0 0 0 ; 0 1 1 0 0 ; 0 0 -1 0 0 ; 0 0 1 1 0 ; 0 0 0 0 1])
    s4 = matspace(K[1 0 0 0 0 ; 0 1 0 z_3 0 ; 0 0 1 1 0 ; 0 0 0 -1 0 ; 0 0 0 1 1])
    s5 = matspace(K[1 0 0 0 0 ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 1 ; 0 0 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    push!(gens, s4)
    push!(gens, s5)

    W = matrix_group(gens)

  elseif n == 34
    K,z_3 = cyclotomic_field(3)
    matspace = matrix_space(K, 6, 6)
    gens = elem_type(matspace)[]

    s1 = matspace(K[-1 0 0 0 0 0 ; 1 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
    s2 = matspace(K[1 1 0 0 0 0 ; 0 -1 0 0 0 0 ; 0 1 1 0 0 0 ; 0 -z_3 - 1 0 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
    s3 = matspace(K[1 0 0 0 0 0 ; 0 1 1 0 0 0 ; 0 0 -1 0 0 0 ; 0 0 1 1 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 1])
    s4 = matspace(K[1 0 0 0 0 0 ; 0 1 0 z_3 0 0 ; 0 0 1 1 0 0 ; 0 0 0 -1 0 0 ; 0 0 0 1 1 0 ; 0 0 0 0 0 1])
    s5 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 1 0 ; 0 0 0 0 -1 0 ; 0 0 0 0 1 1])
    s6 = matspace(K[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 1 ; 0 0 0 0 0 -1])

    push!(gens, s1)
    push!(gens, s2)
    push!(gens, s3)
    push!(gens, s4)
    push!(gens, s5)
    push!(gens, s6)

    W = matrix_group(gens)

  elseif n == 35
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

    W = matrix_group(gens)

  elseif n == 36
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

    W = matrix_group(gens)

  elseif n == 37
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

    W = matrix_group(gens)

  end

  return W

end

function complex_reflection_group_Magma(t::Tuple)

  (m,p,n) = t

  if m == 1 && p == 1
     # In Magma, ShephardTodd(1,1,n) gives the permutation representation, 
     # which is not irreducible...
     error("No model in Magma for the symmetric group case")
  end

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

  W = matrix_group(gens)

  return W

end
