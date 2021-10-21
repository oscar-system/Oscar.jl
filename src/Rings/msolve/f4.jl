import msolve_jll: libneogb
import Libdl: dlopen, dlsym, dlclose

export f4

@doc Markdown.doc"""
    f4(I, <keyword arguments>)

Compute a Groebner basis of the given ideal I w.r.t. to the degree reverse lexicographical monomial ordering using Faugère's F4 algorithm, see [Fau99](@cite) for more information.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.

# Arguments
- `Ì::MPolyIdeal`: input ideal.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `reduce_gb::Int=1`: compute a reduced Gröbner basis for `I`
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> R,(x,y,z) = PolynomialRing(FiniteField(101), ["x","y","z"])
(Multivariate Polynomial Ring in x, y, z over Galois field with characteristic 101, gfp_mpoly[x, y, z])

julia> I = ideal(R, [x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
ideal(x + 2*y + 2*z + 100, x^2 + 100*x + 2*y^2 + 2*z^2, 2*x*y + 2*y*z + 100*y)

julia> f4(I)
Oscar.BiPolyArray{gfp_mpoly}(gfp_mpoly[#undef, #undef, #undef, #undef], Singular Ideal over Singular Polynomial Ring (ZZ/101),(x,y,z),(dp(3),C) with generators (x + 2*y + 2*z - 1, y*z - 19*z^2 + 10*y + 40*z, y^2 - 41*z^2 + 20*y - 20*z, z^3 + 28*z^2 - 37*y + 13*z), Multivariate Polynomial Ring in x, y, z over Galois field with characteristic 101, Singular Polynomial Ring (ZZ/101),(x,y,z),(dp(3),C), true, #undef)
```
"""
function f4(
        I::MPolyIdeal;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        reduce_gb::Int=1,
        info_level::Int=0
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

    I.gb      = BiPolyArray(base_ring(I), basis)
    I.gb.isGB = true
    
    return I.gb
end
