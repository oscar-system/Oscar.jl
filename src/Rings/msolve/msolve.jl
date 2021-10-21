import msolve_jll: libmsolve
import Libdl: dlopen, dlsym, dlclose

export msolve



function get_rational_parametrization(
        nr::Int32,
        lens::Array{Int32,1},
        cfs::Ptr{BigInt}
    )
    C, x  = Nemo.PolynomialRing(Nemo.FlintQQ,"x")
    ctr   = 0

    elim  = 0*x
    for i in 1:lens[1]
        elim  +=  BigInt(unsafe_load(cfs, i))*x^(i-1)
    end
    ctr += lens[1]

    denom = 0*x
    for i in 1:lens[2]
        denom +=  BigInt(unsafe_load(cfs, i+ctr))*x^(i-1)
    end
    ctr +=  lens[2]

    size  = nr-2
    p = Array{Nemo.PolyElem,1}(undef, size)
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

"""
    msolve(I[, initial_hts::Int=17, nr_thrds::Int=1, max_nr_pairs::Int=0,
            la_option::Int=1, infolevel::Int=0, input_file::String="/tmp/in.ms",
            output_file="/tmp/out.ms", precision::Int=67, get_param::Bool=false])

Compute the solution set of the given ideal `I` if the ideal is zero dimensional. At the moment only QQ is supported as ground field.
Sets the dimension of the ideal. If the dimension of `I` is greater then zero the emoty array is returned.

# Arguments
* `I::ideal`: ideal to compute solutions for.
* `initial_hts::Int=17`: hash table size log_2; default is 17, i.e. 2^17 as initial hash
                table size.
* `nr_thrds::Int=1`:  number of threads; default is 1. (not completely supported yet)
* `max_nr_pairs::Int=0`:  maximal number of pairs selected for one F4 matrix; default is
                      0, i.e. no restriction. If matrices get too big or consume
                      too much memory this is a good parameter to play with.
* `la_option::Int=2`: option for linear algebra to be used in F4. there are different linear algebra routines implemented:
    -  `1`: exact sparse-dense computation,
    -  `2`: exact sparse computation, (default)
    - `42`: probabilistic sparse-dense computation,
    - `43`: exact sparse then probabilistic dense computation,
    - `44`: probabilistic sparse computation.
* `info_level::Int=0`: info level for printout:
    - `0`: no printout (default),
    - `1`:  a summary of the computational data is printed at the beginning and the end of the computation,
    - `2`: also dynamical information for each round resp. matrix is printed.
* `input_file::String="/tmp/in.ms"`: input file name for msolve binary; default: /tmp/in.ms.
* `output_file::String="/tmp/in.ms"`: output file name for msolve binary; default: /tmp/out.ms.
* `precision::Int=67`: precision for computing solutions; default is 32.
* `get_param::Bool=false`: get rational parametrization of solution set; default is false.
"""
function msolve(
        I::MPolyIdeal;                        # input generators
        initial_hts::Int=17,                  # hash table size, default 2^17
        nr_thrds::Int=1,                      # number of threads
        max_nr_pairs::Int=0,                  # number of pairs maximally chosen
                                              # in symbolic preprocessing
        la_option::Int=2,                     # linear algebra option
        info_level::Int=0,                    # info level for print outs
        output_file::String="/dev/null",      # msolve output file
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

    variable_names  = Array{String, 1}(undef, nr_vars)
    for i in 1:nr_vars
        variable_names[i] = string(Singular.gens(R)[i])
    end

    field_char  = Singular.characteristic(R)

    #= add new variables and linear forms, =#
    genericity_handling = 2

    reset_ht  = 0
    print_gb  = 0

    #= monomial order defaults to zero =#
    mon_order = 0

    # convert Singular ideal to flattened arrays of ints
    if 0 == field_char
      lens, cfs, exps   = convert_qq_singular_ideal_to_array(J, nr_vars, nr_gens)
    # elseif isprime(field_char)
    #   lens, cfs, exps   = convert_ff_singular_ideal_to_array(J, nr_vars, nr_gens)
    else
        # error("At the moment GroebnerBasis only supports finite fields and the rationals.")
        @error "At the moment msolve only supports the rationala as ground field."
    end
    # lib = Libdl.dlopen("/home/ederc/repos/master-msolve/src/msolve/.libs/libmsolve-0.1.2.so")
    lib = Libdl.dlopen(libmsolve)
    sym = Libdl.dlsym(lib, :msolve_julia)

    res_ld    = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    res_dim   = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    res_dquot = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    res_len   = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    res_cf    = ccall(:malloc, Ptr{Ptr{Cvoid}}, (Csize_t, ), sizeof(Ptr{Cvoid}))
    nb_sols   = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    sols_num  = ccall(:malloc, Ptr{Ptr{Cvoid}}, (Csize_t, ), sizeof(Ptr{Cvoid}))
    sols_den  = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    ccall(sym, Cvoid,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Cvoid}, Ptr{Cint},
         Ptr{Cvoid}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid},
         Ptr{Ptr{Cchar}}, Ptr{Cchar}, Int, Int, Int, Int,
         Int, Int, Int, Int, Int, Int, Int, Int, Int, Int),
        res_ld, res_dim, res_dquot, res_len, res_cf, nb_sols, sols_num, sols_den, lens,
        exps, cfs, variable_names, output_file, field_char, mon_order, nr_vars,
        nr_gens, initial_hts, nr_thrds, max_nr_pairs, reset_ht, la_option,
        print_gb, get_param, genericity_handling, precision, info_level)
    Libdl.dlclose(lib)
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
        return []
    end
    if jl_nb_sols == 0
        @info "The system has no solution."
        return []
    end

    [nterms += jl_len[i] for i=1:jl_ld]
    if 0 == field_char
        res_cf_conv   = unsafe_load(res_cf)
        jl_cf         = reinterpret(Ptr{BigInt}, res_cf_conv)
        sols_num_conv = unsafe_load(sols_num)
        jl_sols_num   = reinterpret(Ptr{BigInt}, sols_num_conv)
        sols_den_conv = unsafe_load(sols_den)
        jl_sols_den   = reinterpret(Ptr{Int32}, sols_den_conv)
    elseif Nemo.isprime(Nemo.FlintZZ(field_char))
        res_cf_conv = unsafe_load(res_cf)
        jl_cf       = reinterpret(Ptr{Int32}, res_cf_conv)
        sols_num_conv = unsafe_load(sols_num)
        jl_sols_num   = reinterpret(Ptr{Int32}, sols_num_conv)
    end

    for i in 1:jl_nb_sols*nr_vars
        tmpcf = BigInt(unsafe_load(jl_sols_num, i))
    end
    solutions = Array{Array{Nemo.fmpq, 1}, 1}(undef, jl_nb_sols)
    for i in 1:jl_nb_sols
        tmp = Array{Nemo.fmpq, 1}(undef, nr_vars)
        for j in 1:nr_vars
            tmp[j]  = BigInt(unsafe_load(jl_sols_num, (i-1)*nr_vars+j)) // BigInt(2)^unsafe_load(jl_sols_den, (i-1)*nr_vars+j)
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
