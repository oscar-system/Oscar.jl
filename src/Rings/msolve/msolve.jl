import msolve_jll: libmsolve

export msolve

@doc Markdown.doc"""
    get_rational_parametrization(nr::Int32, lens::Vector{Int32}, cfs::Ptr{BigInt})

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

    denom = C([unsafe_load(cfs, i+ctr) for i in 1:lens[2]])
    ctr +=  lens[2]

    size  = nr-2
    p = Vector{PolyElem}(undef, size)
    c = Vector{fmpz}(undef, size)
    k = 1
    for i in 3:nr
        p[k]  = C([unsafe_load(cfs, j+ctr) for j in 1:lens[i]-1])
        c[k]  =   (-1) * fmpz(unsafe_load(cfs, lens[i]+ctr))
        ctr   +=  lens[i]
        k     +=  1
    end

    return elim, denom, p, c
end

@doc Markdown.doc"""
    msolve(I::MPolyIdeal, <keyword arguments>)

Given an ideal `I` with a finite solution set over the complex numbers, return a pair `p,r` where `p` is the rational parametrization of the solution set and `r` represents the real roots of `Ì`  with a given precision (default 32 bits).
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

# Examples
```jldoctest
julia> R,(x1,x2,x3) = PolynomialRing(QQ, ["x1","x2","x3"])
(Multivariate Polynomial Ring in x1, x2, x3 over Rational Field, fmpq_mpoly[x1, x2, x3])

julia> I = ideal(R, [x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
ideal(x1 + 2*x2 + 2*x3 - 1, x1^2 - x1 + 2*x2^2 + 2*x3^2, 2*x1*x2 + 2*x2*x3 - x2)

julia> msolve(I)
((84*x^4 - 40*x^3 + x^2 + x, 336*x^3 - 120*x^2 + 2*x + 1, PolyElem[-184*x^3 + 80*x^2 - 4*x - 1, -36*x^3 + 18*x^2 - 2*x], fmpz[-1, -1]), Vector{fmpq}[[744483363399261433351//1180591620717411303424, 372241681699630716673//1180591620717411303424, -154187553040555781639//1180591620717411303424], [1, 0, 0], [71793683196126133110381699745//316912650057057350374175801344, 71793683196126133110381699745//633825300114114700748351602688, 173325283664805084153412401855//633825300114114700748351602688], [196765270119568550571//590295810358705651712, 1//590295810358705651712, 196765270119568550571//590295810358705651712]])
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
        precision::Int=32                     # precision of the solution set
        )
    R = base_ring(I)
    nr_vars     = nvars(R)
    nr_gens     = ngens(I)
    field_char  = Int(characteristic(R))

    variable_names = map(string, Singular.symbols(R))

    field_char  = Singular.characteristic(R)

    #= add new variables and linear forms, =#
    genericity_handling = 2

    reset_ht  = 0
    print_gb  = 0
    get_param = 1

    mon_order       = 0
    elim_block_size = 0
    use_signatures  = 0

    if field_char != 0
        error("At the moment msolve only supports the rationals as ground field.")
    end
    # convert Singular ideal to flattened arrays of ints
    lens, cfs, exps   = convert_oscar_ideal_to_array(I)

    res_ld      = Ref(Cint(0))
    res_nr_vars = Ref(Cint(0))
    res_dim     = Ref(Cint(0))
    res_dquot   = Ref(Cint(0))
    nb_sols     = Ref(Cint(0))
    res_len     = Ref(Ptr{Cint}(0))
    res_vnames  = Ref(Ptr{Ptr{Cchar}}(0))
    res_cf_lf   = Ref(Ptr{Cvoid}(0))
    res_cf      = Ref(Ptr{Cvoid}(0))
    sols_num    = Ref(Ptr{Cvoid}(0))
    sols_den    = Ref(Ptr{Cint}(0))

    ccall((:msolve_julia, libmsolve), Cvoid,
        (Ptr{Nothing}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}},
        Ptr{Ptr{Ptr{Cchar}}}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cint}, Ptr{Cvoid}, Ptr{Ptr{Cint}}, Ptr{Cint},
        Ptr{Cint}, Ptr{Cvoid}, Ptr{Ptr{Cchar}}, Ptr{Cchar}, Cint, Cint,
        Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint,
        Cint, Cint, Cint),
        cglobal(:jl_malloc), res_ld, res_nr_vars, res_dim, res_dquot, res_len, res_vnames,
        res_cf_lf, res_cf,
        nb_sols, sols_num, sols_den, lens, exps, cfs, variable_names,
        "/dev/null", field_char, mon_order, elim_block_size, nr_vars,
        nr_gens, initial_hts, nr_thrds, max_nr_pairs, reset_ht, la_option,
        use_signatures, print_gb, get_param, genericity_handling, precision,
        info_level)
    # convert to julia array, also give memory management to julia
    jl_ld       = res_ld[]

    jl_dim      = res_dim[]
    jl_dquot    = res_dquot[]

    jl_nb_sols	= nb_sols[]

    # set dimension
    I.dim = jl_dim
    if jl_dim > 0
        error("Dimension of ideal is greater than zero, no solutions provided.")
    end

    # get rational parametrization
    jl_len  = Base.unsafe_wrap(Array, res_len[], jl_ld)
    nterms  = 0

    [nterms += jl_len[i] for i=1:jl_ld]
    # if 0 == field_char
    jl_cf       = reinterpret(Ptr{BigInt}, res_cf[])

    rat_param = get_rational_parametrization(jl_ld, jl_len, jl_cf)

    if jl_nb_sols == 0
        return rat_param, Vector{fmpq}[]
    end

    # get solutions
    jl_sols_num = reinterpret(Ptr{BigInt}, sols_num[])
    jl_sols_den = reinterpret(Ptr{Int32}, sols_den[])
    # elseif is_prime(field_char)
    #     jl_cf       = reinterpret(Ptr{Int32}, res_cf[])
    #     jl_sols_num   = reinterpret(Ptr{Int32}, sols_num[])
    # end

    #= solutions are returned as intervals, i.e. a minimum and a maximum entry for
     = the numerator and denominator; thus we sum up and divide by  =#
    solutions = Vector{Vector{fmpq}}(undef, jl_nb_sols)

    len = 2*jl_nb_sols*nr_vars
    i = 1
    k = 1
    while i <= len
        j = 1
        tmp = Vector{Nemo.fmpq}(undef, nr_vars)
        while j <= nr_vars
            tmp[j]  = fmpq(unsafe_load(jl_sols_num, i)) >> Int64(unsafe_load(jl_sols_den, i))
            tmp[j] += fmpq(unsafe_load(jl_sols_num, i+1)) >> Int64(unsafe_load(jl_sols_den, i+1))
            tmp[j] = tmp[j] >> 1
            i += 2
            j += 1
        end
        solutions[k] = tmp
        k += 1
    end

    ccall((:free_msolve_julia_result_data, libmsolve), Nothing ,
        (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
        Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cint}}, Cint, Cint, Cint),
        cglobal(:jl_free), res_len, res_cf, sols_num, sols_den,
        jl_ld, jl_nb_sols, field_char)

    return rat_param, solutions
end
