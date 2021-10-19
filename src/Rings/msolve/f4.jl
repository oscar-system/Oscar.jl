import msolve_jll: libneogb
import Libdl: dlopen, dlsym, dlclose

"""
    f4(I[, hts::Int=17, nthrds::Int=1, maxpairs::Int=0, resetht::Int=0,
            laopt::Int=1, reducegb::Int=0, pbmfiles::Int=0,
            infolevel::Int=0, monorder::Symbol=:degrevlex])

Compute a Groebner basis of the given ideal I w.r.t. to the given monomial
order using Faugere's F4 algorithm. The function takes a Singular ideal as
input and returns a Singular ideal. At the moment only finite fields up to 31-bit and the degree reverse lexicographical ordering is supported.

# Arguments
* `I::Singular.sideal`: ideal to compute a Groebner basis for.
* `initial_hts::Int=17`: hash table size log_2; default is 17, i.e. 2^17 as initial hash
                table size.
* `nr_thrds::Int=1`:  number of threads; default is 1.
* `max_nr_pairs::Int=0`:  maximal number of pairs selected for one matrix; default is
                      0, i.e. no restriction. If matrices get too big or consume
                      too much memory this is a good parameter to play with.
* `la_option::Int=1`: option for linear algebra to be used. there are different
                  linear algebra routines implemented:
    -  `1`: exact sparse-dense computation (default),
    -  `2`: exact sparse computation,
    - `42`: probabilistic sparse-dense computation,
    - `43`: exact sparse then probabilistic dense computation,
    - `44`: probabilistic sparse computation.
* `reduce_gb::Int=1`:  reduce final basis; default is 1.
* `info_level::Int=0`: info level for printout:
    - `0`: no printout (default),
    - `1`:  a summary of the computational data is printed at the beginning and
    the end of the computation,
    - `2`: also dynamical information for each round resp. matrix is printed.
"""
function f4(
        I::MPolyIdeal;                # input generators
        initial_hts::Int=17,                  # hash table size, default 2^17
        nr_thrds::Int=1,                # number of threads
        max_nr_pairs::Int=0,              # number of pairs maximally chosen
                                      # in symbolic preprocessing
        la_option::Int=2,                 # linear algebra option
        reduce_gb::Int=1,             # reduce final basis
        info_level::Int=0              # info level for print outs
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

    mon_order   = 0
    field_char  = Singular.characteristic(R)

    # convert Singular ideal to flattened arrays of ints
    if isprime(field_char)
      lens, cfs, exps = convert_ff_singular_ideal_to_array(J, nr_vars, nr_gens)
    else
        error("At the moment f4 only supports finite fields.")
    end
    #= dir = joinpath(dirname(pathof(GroebnerBasis)),"../deps")
     = lib = Libdl.dlopen("$dir/libmsolve.so.0.2.0") =#
    lib = Libdl.dlopen(libneogb)
    sym = Libdl.dlsym(lib, :f4_julia)

    gb_ld   = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    gb_len  = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    gb_exp  = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    gb_cf   = ccall(:malloc, Ptr{Ptr{Cvoid}}, (Csize_t, ), sizeof(Ptr{Cvoid}))
    nr_terms  = ccall(sym, Int,
        (Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
          Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Int, Int, Int, Int, Int,
          Int, Int, Int, Int, Int, Int, Int),
        gb_ld, gb_len, gb_exp, gb_cf, lens, exps, cfs, field_char, mon_order, nr_vars,
        nr_gens, initial_hts, nr_thrds, max_nr_pairs, 0, la_option, reduce_gb, 0, info_level)

    if nr_terms == 0
        error("Something went wrong in the C code of F4.")
    end
    # convert to julia array, also give memory management to julia
    jl_ld   = unsafe_load(gb_ld)
    jl_len  = Base.unsafe_wrap(Array, unsafe_load(gb_len), jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, unsafe_load(gb_exp), nr_terms*nr_vars)
    # if 0 == char
    #     gb_cf_conv  = unsafe_load(gb_cf)
    #     jl_cf       = reinterpret(Ptr{BigInt}, gb_cf_conv)
    # elseif Nemo.isprime(Nemo.FlintZZ(char))
      gb_cf_conv  = Ptr{Ptr{Int32}}(gb_cf)
      jl_cf       = Base.unsafe_wrap(Array, unsafe_load(gb_cf_conv), nr_terms)
    # end

    # construct Singular ideal
    # if 0 == char
    #     basis = convert_qq_gb_array_to_singular_ideal(
    #       jl_ld, jl_len, jl_exp, jl_cf, R)
    # else
        basis = convert_ff_gb_array_to_singular_ideal(
          jl_ld, jl_len, jl_exp, jl_cf, R)
    # end
    sym = Libdl.dlsym(lib, :free_julia_data)
    ccall(sym, Nothing , (Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
                Int, Int), gb_len, gb_exp, gb_cf, jl_ld, field_char)
    # free data
    ccall(:free, Nothing , (Ptr{Cint}, ), gb_ld)

    # for letting the garbage collector free memory
    jl_len      = Nothing
    jl_exp      = Nothing
    gb_cf_conv  = Nothing
    jl_cf       = Nothing

    I.gb = BiPolyArray(base_ring(I), basis)
    
    return I.gb
end
