test_sample_weights = Dict(
    1 => [
        [QQ(  8,10)],
        [QQ(  7, 2)],
        [QQ(-10, 4)],
        [QQ(-10, 8)],
        [QQ( -6, 1)],
    ],

    2 => [
        [QQ(  3, 9), QQ(  2, 4)],
        [QQ( -1, 2), QQ( -8, 8)],
        [QQ( 10,10), QQ( -2,10)],
        [QQ( -2, 3), QQ( -4,10)],
        [QQ(  2, 1), QQ( -4, 3)],
    ],

    3 => [
        [QQ(  8, 9), QQ(  7, 1), QQ( 10, 8)],
        [QQ( 10, 3), QQ( -1, 6), QQ( -4, 7)],
        [QQ( -5, 3), QQ(  9, 9), QQ( -2, 9)],
        [QQ(  9, 5), QQ( -6,10), QQ(  0, 2)],
        [QQ( -5,10), QQ(  3, 2), QQ( -2,10)],
    ],

    4 => [
        [QQ(  7,10), QQ( -7, 9), QQ(  4, 9), QQ( -2, 5)],
        [QQ( -9, 5), QQ( -1, 9), QQ(  7, 6), QQ( -9, 5)],
        [QQ( -4, 8), QQ( -1, 7), QQ( -2, 3), QQ(  3, 6)],
        [QQ(  5, 4), QQ( -5, 9), QQ(  9, 4), QQ(  0, 6)],
        [QQ( -7, 6), QQ( -5, 2), QQ( -6, 3), QQ(  0, 7)],
    ],

    5 => [
        [QQ( -5, 5), QQ(  9, 2), QQ(  2, 3), QQ(  3, 1), QQ( -4, 8)],
        [QQ(  7, 2), QQ(  1, 8), QQ( -2, 7), QQ(  5, 9), QQ( -2, 3)],
        [QQ(  8, 7), QQ(  6, 8), QQ( -9, 2), QQ(  1, 7), QQ( -5,10)],
        [QQ(  7, 3), QQ(  0, 7), QQ( -4,10), QQ(-10, 2), QQ( -4, 4)],
        [QQ( -1, 3), QQ( -7, 3), QQ(-10, 9), QQ( -6, 2), QQ(  3, 8)],
    ],

    6 => [
        [QQ(  1, 3), QQ(-10, 3), QQ(  4, 2), QQ( 10, 5), QQ( -8, 7), QQ( -2, 6)],
        [QQ( -9, 4), QQ(  1, 7), QQ( -7, 9), QQ(  8, 1), QQ( -1, 1), QQ(  8, 2)],
        [QQ( -4, 3), QQ( -3, 8), QQ( -7, 9), QQ( -9, 4), QQ( -1, 2), QQ( -4, 3)],
        [QQ(  8, 4), QQ(  9, 9), QQ(-10, 4), QQ(  2, 1), QQ(  7,10), QQ(  8, 4)],
        [QQ(  3, 4), QQ(  4, 2), QQ( -5, 9), QQ( -7, 2), QQ( -7, 7), QQ( -6, 6)],
    ],

    7 => [
        [QQ(  6, 7), QQ(  9, 3), QQ(-10, 1), QQ(  5,10), QQ(  8, 4), QQ( -6, 8), QQ(  4, 2)],
        [QQ(  1, 1), QQ( -3,10), QQ(  5, 2), QQ( -6, 4), QQ(  5, 3), QQ(  3, 1), QQ( -8, 6)],
        [QQ(  9, 8), QQ(  3, 1), QQ( -1, 3), QQ( -1, 2), QQ( -2, 7), QQ(  5, 1), QQ(  8, 8)],
        [QQ( -2, 5), QQ( -6, 8), QQ( -4, 8), QQ(  3, 1), QQ(-10, 3), QQ(  4, 8), QQ(  6, 8)],
        [QQ(  5, 6), QQ(  5, 1), QQ(  4, 1), QQ(  4, 4), QQ(  4, 1), QQ(  4, 3), QQ(  4, 1)],
    ],

    8 => [
        [QQ(  6, 6), QQ( -4, 6), QQ( -8,10), QQ(  8, 2), QQ(  7, 4), QQ( -4, 3), QQ(  9, 3), QQ(  3, 7)],
        [QQ(  4, 8), QQ( -2, 1), QQ( -2, 6), QQ(  6, 3), QQ(  3, 7), QQ(  7, 5), QQ( -8, 3), QQ( -5, 7)],
        [QQ(  6, 6), QQ(  4, 6), QQ( -5, 8), QQ(  4,10), QQ(  1, 4), QQ(  0, 7), QQ( -1, 4), QQ(  2, 3)],
        [QQ(  4, 7), QQ(-10, 6), QQ(-10, 6), QQ( 10, 2), QQ(  7, 9), QQ(  6,10), QQ( -5, 7), QQ(  2,10)],
        [QQ(  0, 2), QQ(  7, 7), QQ( -9, 6), QQ(  9, 4), QQ(  1, 3), QQ(  9, 4), QQ(  0, 5), QQ(  6, 8)],
    ],

    9 => [
        [QQ(  4, 8), QQ(  6, 4), QQ(  8, 1), QQ(  6, 6), QQ(  1, 2), QQ( -6,10), QQ(  1, 6), QQ(  9, 3), QQ(-10,10)],
        [QQ(  4, 9), QQ(  9, 5), QQ(  8, 7), QQ( -3, 2), QQ(  5, 6), QQ( -8, 6), QQ( 10, 9), QQ( -4, 6), QQ( 10, 6)],
        [QQ( -7,10), QQ(  8, 6), QQ( -8, 8), QQ( -5, 2), QQ(  2, 4), QQ( -6, 4), QQ( -5, 7), QQ( -4, 2), QQ( -5, 6)],
        [QQ(  0, 2), QQ( -1, 2), QQ(  8, 2), QQ(  7, 7), QQ(  2, 2), QQ(  5, 5), QQ( -5, 7), QQ(  3, 5), QQ( -3, 9)],
        [QQ(  1, 7), QQ( -5, 4), QQ( -6, 9), QQ( -4, 6), QQ(  8, 7), QQ(  0, 8), QQ( -3,10), QQ( -5, 8), QQ( 10, 5)],
    ],

    10 => [
        [QQ(  5, 7), QQ(  8, 5), QQ(  1, 1), QQ(-10, 4), QQ(  4, 6), QQ(-10, 4), QQ( -7, 3), QQ( -6, 5), QQ( -4, 1), QQ( -7, 4)],
        [QQ( -2,10), QQ(  1, 4), QQ( -5, 6), QQ( -1, 9), QQ( -2, 1), QQ(  5, 9), QQ(  2, 6), QQ( -5, 6), QQ( -2, 1), QQ(  8, 9)],
        [QQ(  8,10), QQ( -9, 7), QQ(  0, 6), QQ(  5, 5), QQ( -8, 2), QQ(  5, 4), QQ(  8, 2), QQ(  8, 7), QQ( -4, 9), QQ( -1, 1)],
        [QQ( -5, 9), QQ( -6, 8), QQ( -8, 8), QQ( -7, 7), QQ(  8, 3), QQ(  7, 3), QQ( -9, 5), QQ(-10, 2), QQ( -7, 6), QQ( -2, 3)],
        [QQ(  8, 1), QQ( -9,10), QQ( -9, 8), QQ( -5, 4), QQ(  7,10), QQ(  9, 3), QQ(  3, 5), QQ(  0, 3), QQ(  4, 9), QQ( -6, 8)],
    ],

)

# Code to generate fixed test-examples, that can be copy-pasted into console or files:

# using Random

# function generate_random_QQ()
#     numerator = rand(-10:10)
#     denominator = rand(1:10)
#     return ("QQ($numerator,$denominator)", numerator, denominator)
# end

# function pad_QQs(qq_list, max_len_numerator, max_len_denominator)
#     vec = [qq[1] for qq in qq_list]
#     for i in 1:length(vec)
#         parts = split(vec[i], ',')
#         num = parts[1][4:end]
#         denom = parts[2][1:end-1]
#         vec[i] = "QQ(" * " "^(max_len_numerator - length(num)) * num * "," * " "^(max_len_denominator - length(denom)) * denom * ")"
#     end
#     return vec
# end

# function generate_examples()
#    all_qqs = []
    
    # First generate all QQs and record their lengths
#     for i in 1:10
#         for _ in 1:5
#             push!(all_qqs, [generate_random_QQ() for _ in 1:i])
#         end
#     end

#     max_len_numerator = maximum(length(string(qq[2])) for qq_list in all_qqs for qq in qq_list)
#     max_len_denominator = maximum(length(string(qq[3])) for qq_list in all_qqs for qq in qq_list)
    
#     result = "test_sample_weights = Dict(\n"
#     idx = 1
#     for i in 1:10
#         result *= "    $i => [\n"
        
#         lines = []
#         for _ in 1:5
#            vec = pad_QQs(all_qqs[idx], max_len_numerator, max_len_denominator)
#             push!(lines, "        [" * join(vec, ", ") * "]")
#             idx += 1
#         end

        # Add each line to the result
#         for line in lines
#             result *= "$line,\n"
#         end

#         result *= "    ],\n\n"
#     end
#     result *= ")"
#     return result
# end

# println(generate_examples())
