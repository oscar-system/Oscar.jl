import msolve_jll: libneogb
import Libdl: dlopen, dlsym, dlclose

"""
    f4(I[, hts::Int=17, nthrds::Int=1, maxpairs::Int=0, resetht::Int=0,
            laopt::Int=1, reducegb::Int=0, pbmfiles::Int=0,
            infolevel::Int=0, monorder::Symbol=:degrevlex])

Compute a Groebner basis of the given ideal I w.r.t. to the given monomial
order using Faugere's F4 algorithm. The function takes a Singular ideal as
input and returns a Singular ideal. At the moment only finite fields up to
31-bit and the rationals are supported as ground fields.

# Arguments
* `I::Singular.sideal`: ideal to compute a Groebner basis for.
* `hts::Int=17`: hash table size log_2; default is 17, i.e. 2^17 as initial hash
                table size.
* `nthrds::Int=1`:  number of threads; default is 1.
* `maxpairs::Int=0`:  maximal number of pairs selected for one matrix; default is
                      0, i.e. no restriction. If matrices get too big or consume
                      too much memory this is a good parameter to play with.
* `resetht::Int=0`: Resets the hash table after `resetht` steps in the algorthm;
                    default is 0, i.e. no resetting at all. Since we add
                    monomials to the matrices which are only used for reduction
                    purposes, but have no further meaning in the basis, this
                    parameter might also help when memory get a problem.
* `laopt::Int=1`: option for linear algebra to be used. there are different
                  linear algebra routines implemented:
    -  `1`: exact sparse-dense computation (default),
    -  `2`: exact sparse computation,
    - `42`: probabilistic sparse-dense computation,
    - `43`: exact sparse then probabilistic dense computation,
    - `44`: probabilistic sparse computation.
* `reducegb::Int=0`:  reduce final basis; default is 0. Note that for
                      computations over Q we do not normalize the polynomials,
                      the basis is only minimal and tailreduced. Normalize by
                      hand if necessary.
* `pbmfiles::Int=0`: option for generating pbm files of matrices:
    - `0`: off (default),
    - `1`:  on.
* `infolevel::Int=0`: info level for printout:
    - `0`: no printout (default),
    - `1`:  a summary of the computational data is printed at the beginning and
    the end of the computation,
    - `2`: also dynamical information for each round resp. matrix is printed.
* `monorder::Symbol=:degrevlex`: monomial order w.r.t. which the computation is
                                done;
    - `degrevlex`: the degree-reverse-lexicographical (DRL) order (default),
    - `lex`: the lexicographical order (LEX).
"""
function f4(
        I::Singular.sideal;           # input generators
        hts::Int=17,                  # hash table size, default 2^17
        nthrds::Int=1,                # number of threads
        maxpairs::Int=0,              # number of pairs maximally chosen
                                      # in symbolic preprocessing
        resetht::Int=0,               # resetting global hash table
        laopt::Int=2,                 # linear algebra option
        reducegb::Int=0,              # reduce final basis
        pbmfiles::Int=0,              # generation of pbm files
        infolevel::Int=0,             # info level for print outs
        monorder::Symbol=:dregrevlex  # monomial order
        )
    R     = I.base_ring
    # skip zero generators in ideal
    ptr = Singular.libSingular.id_Copy(I.ptr, R.ptr)
    J   = Singular.Ideal(R, ptr)
    Singular.libSingular.idSkipZeroes(J.ptr)
    # get number of variables
    nvars   = Singular.nvars(R)
    ngens   = Singular.ngens(J)

    ord = 0
    if monorder == :degrevlex
        ord = 0
    end
    if monorder == :lex
        ord = 1
    end
    char  = Singular.characteristic(R)

    # convert Singular ideal to flattened arrays of ints
    if 0 == char
      lens, cfs, exps   = convert_qq_singular_ideal_to_array(J, nvars, ngens)
    elseif Nemo.isprime(Nemo.FlintZZ(char))
      lens, cfs, exps   = convert_ff_singular_ideal_to_array(J, nvars, ngens)
    else
        error("At the moment GroebnerBasis only supports finite fields and the rationals.")
    end
    #= dir = joinpath(dirname(pathof(GroebnerBasis)),"../deps")
     = lib = Libdl.dlopen("$dir/libmsolve.so.0.2.0") =#
    lib = Libdl.dlopen("/home/ederc/repos/master-msolve/src/neogb/.libs/libneogb-0.1.2.so")
    # lib = Libdl.dlopen(libneogb)
    sym = Libdl.dlsym(lib, :f4_julia)

    gb_ld   = ccall(:malloc, Ptr{Cint}, (Csize_t, ), sizeof(Cint))
    gb_len  = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    gb_exp  = ccall(:malloc, Ptr{Ptr{Cint}}, (Csize_t, ), sizeof(Ptr{Cint}))
    gb_cf   = ccall(:malloc, Ptr{Ptr{Cvoid}}, (Csize_t, ), sizeof(Ptr{Cvoid}))
    nterms  = ccall(sym, Int,
        (Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
          Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Int, Int, Int, Int, Int,
          Int, Int, Int, Int, Int, Int, Int),
        gb_ld, gb_len, gb_exp, gb_cf, lens, exps, cfs, char, ord, nvars,
        ngens, hts, nthrds, maxpairs, resetht, laopt, reducegb, pbmfiles, infolevel)

    if nterms == 0
        error("Something went wrong in the C code of F4.")
    end
    # convert to julia array, also give memory management to julia
    jl_ld   = unsafe_load(gb_ld)
    jl_len  = Base.unsafe_wrap(Array, unsafe_load(gb_len), jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, unsafe_load(gb_exp), nterms*nvars)
    if 0 == char
        gb_cf_conv  = unsafe_load(gb_cf)
        jl_cf       = reinterpret(Ptr{BigInt}, gb_cf_conv)
    elseif Nemo.isprime(Nemo.FlintZZ(char))
      gb_cf_conv  = Ptr{Ptr{Int32}}(gb_cf)
      jl_cf       = Base.unsafe_wrap(Array, unsafe_load(gb_cf_conv), nterms)
    end

    # construct Singular ideal
    if 0 == char
        basis = convert_qq_gb_array_to_singular_ideal(
          jl_ld, jl_len, jl_exp, jl_cf, R)
    elseif Nemo.isprime(Nemo.FlintZZ(char))
        basis = convert_ff_gb_array_to_singular_ideal(
          jl_ld, jl_len, jl_exp, jl_cf, R)
    end
    sym = Libdl.dlsym(lib, :free_julia_data)
    ccall(sym, Nothing , (Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
                Int, Int), gb_len, gb_exp, gb_cf, jl_ld, char)
    # free data
    ccall(:free, Nothing , (Ptr{Cint}, ), gb_ld)

    # for letting the garbage collector free memory
    jl_len      = Nothing
    jl_exp      = Nothing
    gb_cf_conv  = Nothing
    jl_cf       = Nothing

    basis.isGB  = true;
    

    return basis
end
