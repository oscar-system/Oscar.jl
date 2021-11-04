import msolve_jll: libmsolve

export msolve

@doc Markdown.doc"""
    function get_rational_parametrization(nr::Int32, lens::Array{Int32,1}, cfs::Ptr{BigInt})

Construct the rational parametrization of the solution set computed via msolve.

**Note**: This is an internal function and should only be used inside `msolve()`.
"""
function get_rational_parametrization(
        nr::Int32,
        lens::Vector{Int32},
        cfs::Ptr{BigInt}
    )
    C, x  = PolynomialRing(QQ,"x")
    ctr   = 0

    elim  = C([unsafe_load(cfs, i) for i in 1:lens[1]])
    ctr += lens[1]

    denom = 0*x
    for i in 1:lens[2]
        denom +=  BigInt(unsafe_load(cfs, i+ctr))*x^(i-1)
    end
    ctr +=  lens[2]

    size  = nr-2
    p = Array{PolyElem,1}(undef, size)
    c = Array{BigInt,1}(undef, size)
    k = 1
    for i in 3:nr
        p[k]  = 0*x
        for j in 1:lens[i]-1
            p[k]  +=  BigInt(unsafe_load(cfs, j+ctr))*x^(j-1)
        end
        c[k]  =   (-1) * BigInt(unsafe_load(cfs, lens[i]+ctr))
        ctr   +=  lens[i]
        k     +=  1
    end

    return elim, denom, p, c
end

@doc Markdown.doc"""
    msolve(I, <keyword arguments>)

Compute the real roots of the ideal `I` with a given precision (default 32 bits) if the solution set over the
complex numbers is finite.
See [BES21](@cite) for more information.

**Note**: At the moment only QQ is supported as ground field. If the dimension of `I`
is greater then zero an empty array is returned.

# Arguments
- `Ì::MPolyIdeal`: input ideal.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).
- `precision::Int=32`: bit precision for the computed solutions.
- `get_param::Bool=false`: return rational parametrization of the solution set over the complex numbers.

# Examples
```jldoctest
julia> R,(x1,x2,x3) = PolynomialRing(QQ, ["x1","x2","x3"])
(Multivariate Polynomial Ring in x1, x2, x3 over Rational Field, fmpq_mpoly[x1, x2, x3])

julia> I = ideal(R, [x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
ideal(x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2, 2*x1*x2 + 2*x2*x3 - x2)

julia> msolve(I)
4-element Vector{Vector{fmpq}}:
 [3197531698215911246794018079661//5070602400912917605986812821504, 1598765849107955623397009039829//5070602400912917605986812821504, -662230497759452443800611668909//5070602400912917605986812821504]
 [1, 0, 0]
 [143587366392252266220763399489//633825300114114700748351602688, 71793683196126133110381699745//633825300114114700748351602688, 173325283664805084153412401855//633825300114114700748351602688]
 [52818775009509558395695966891//158456325028528675187087900672, 3//2535301200456458802993406410752, 845100400152152934331135470251//2535301200456458802993406410752]

julia> msolve(I, get_param=true)
((84*x^4 - 40*x^3 + x^2 + x, 336*x^3 - 120*x^2 + 2*x + 1, PolyElem[-184*x^3 + 80*x^2 - 4*x - 1, -36*x^3 + 18*x^2 - 2*x], BigInt[-1, -1]), Vector{fmpq}[[3197531698215911246794018079661//5070602400912917605986812821504, 1598765849107955623397009039829//5070602400912917605986812821504, -662230497759452443800611668909//5070602400912917605986812821504], [1, 0, 0], [143587366392252266220763399489//633825300114114700748351602688, 71793683196126133110381699745//633825300114114700748351602688, 173325283664805084153412401855//633825300114114700748351602688], [52818775009509558395695966891//158456325028528675187087900672, 3//2535301200456458802993406410752, 845100400152152934331135470251//2535301200456458802993406410752]])
```
"""
function msolve(
        I::MPolyIdeal;                        # input generators
        initial_hts::Int=17,                  # hash table size, default 2^17
        nr_thrds::Int=1,                      # number of threads
        max_nr_pairs::Int=0,                  # number of pairs maximally chosen
                                              # in symbolic preprocessing
        la_option::Int=2,                     # linear algebra option
        info_level::Int=0,                    # info level for print outs
        precision::Int=32,                    # precision of the solution set
        get_param::Bool=false                 # return rational parametrization of
                                              # solution set
        )
    singular_assure(I)
    SI    = I.gens.S
    R     = I.gens.Sx
    # skip zero generators in ideal
    ptr = Singular.libSingular.id_Copy(SI.ptr, R.ptr)
    J   = Singular.Ideal(R, ptr)
    Singular.libSingular.idSkipZeroes(J.ptr)
    # get number of variables
    nr_vars = Singular.nvars(R)
    nr_gens = Singular.ngens(J)
    vars    = Singular.gens(R)

    variable_names = map(string, Singular.symbols(R))

    field_char  = Singular.characteristic(R)

    #= add new variables and linear forms, =#
    genericity_handling = 2

    reset_ht  = 0
    print_gb  = 0

    #= monomial order defaults to zero =#
    mon_order = 0

    if field_char != 0
        error("At the moment msolve only supports the rationals as ground field.")
    end
    # convert Singular ideal to flattened arrays of ints
    lens, cfs, exps   = convert_singular_ideal_to_array(J)

    res_ld    = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    res_dim   = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    res_dquot = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    res_len   = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    res_cf    = ccall(:malloc, Ptr{Ptr{Cvoid}}, (Csize_t, ), sizeof(Ptr{Cvoid}))
    nb_sols   = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    sols_num  = ccall(:malloc, Ptr{Ptr{Cvoid}}, (Csize_t, ), sizeof(Ptr{Cvoid}))
    sols_den  = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    ccall((:msolve_julia, libmsolve), Cvoid,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Cvoid}, Ptr{Cint},
         Ptr{Cvoid}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid},
         Ptr{Ptr{Cchar}}, Ptr{Cchar}, Int, Int, Int, Int,
         Int, Int, Int, Int, Int, Int, Int, Int, Int, Int),
        res_ld, res_dim, res_dquot, res_len, res_cf, nb_sols, sols_num, sols_den, lens,
        exps, cfs, variable_names, "/dev/null", field_char, mon_order, nr_vars,
        nr_gens, initial_hts, nr_thrds, max_nr_pairs, reset_ht, la_option,
        print_gb, get_param, genericity_handling, precision, info_level)
    # convert to julia array, also give memory management to julia
    jl_ld       = unsafe_load(res_ld)

    jl_dim      = unsafe_load(res_dim)
    jl_dquot    = unsafe_load(res_dquot)

    jl_nb_sols  = unsafe_load(nb_sols)
    jl_len      = Base.unsafe_wrap(Array, unsafe_load(res_len), jl_ld)

    nterms  = 0
    
    # set dimension
    I.dim = jl_dim
    if jl_dim > 0
        @info "Dimension is greater than zero, no solutions provided."
        return Vector{fmpq}[]
    end
    if jl_nb_sols == 0
        @info "The system has no solution."
        return Vector{fmpq}[]
    end

    [nterms += jl_len[i] for i=1:jl_ld]
    # if 0 == field_char
        res_cf_conv   = unsafe_load(res_cf)
        jl_cf         = reinterpret(Ptr{BigInt}, res_cf_conv)
        sols_num_conv = unsafe_load(sols_num)
        jl_sols_num   = reinterpret(Ptr{BigInt}, sols_num_conv)
        sols_den_conv = unsafe_load(sols_den)
        jl_sols_den   = reinterpret(Ptr{Int32}, sols_den_conv)
    # elseif isprime(field_char)
    #     res_cf_conv = unsafe_load(res_cf)
    #     jl_cf       = reinterpret(Ptr{Int32}, res_cf_conv)
    #     sols_num_conv = unsafe_load(sols_num)
    #     jl_sols_num   = reinterpret(Ptr{Int32}, sols_num_conv)
    # end

    for i in 1:jl_nb_sols*nr_vars
        tmpcf = BigInt(unsafe_load(jl_sols_num, i))
    end
    solutions = Array{Array{fmpq, 1}, 1}(undef, jl_nb_sols)
    for i in 1:jl_nb_sols
        tmp = Array{Nemo.fmpq, 1}(undef, nr_vars)
        for j in 1:nr_vars
            tmp[j]  = fmpq(unsafe_load(jl_sols_num, (i-1)*nr_vars+j)) >> unsafe_load(jl_sols_den, (i-1)*nr_vars+j)
        end
        solutions[i]  = tmp
    end
    if get_param 
        rat_param = get_rational_parametrization(jl_ld, jl_len, jl_cf)
        return rat_param, solutions
    else
        return solutions
    end
end
