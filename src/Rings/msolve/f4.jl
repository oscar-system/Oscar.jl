import msolve_jll: libneogb

export f4

@doc Markdown.doc"""
    f4(I::MPolyIdeal, <keyword arguments>)

Compute a Groebner basis of the given ideal `I` w.r.t. to the degree reverse lexicographical monomial ordering using Faugère's F4 algorithm.
See [Fau99](@cite) for more information.

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
julia> R,(x,y,z) = PolynomialRing(GF(101), ["x","y","z"])
(Multivariate Polynomial Ring in x, y, z over Galois field with characteristic 101, gfp_mpoly[x, y, z])

julia> I = ideal(R, [x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
ideal(x + 2*y + 2*z + 100, x^2 + 100*x + 2*y^2 + 2*z^2, 2*x*y + 2*y*z + 100*y)

julia> f4(I)
Oscar.BiPolyArray{gfp_mpoly}(gfp_mpoly[#undef, #undef, #undef, #undef], Singular ideal over Singular Polynomial Ring (ZZ/101),(x,y,z),(dp(3),C) with generators (x + 2*y + 2*z - 1, y*z - 19*z^2 + 10*y + 40*z, y^2 - 41*z^2 + 20*y - 20*z, z^3 + 28*z^2 - 37*y + 13*z), Multivariate Polynomial Ring in x, y, z over Galois field with characteristic 101, Singular Polynomial Ring (ZZ/101),(x,y,z),(dp(3),C), true, #undef, true)
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

    mon_order       = 0
    elim_block_size = 0
    field_char  = Singular.characteristic(R)

    # convert Singular ideal to flattened arrays of ints
    if !(isprime(field_char))
        error("At the moment f4 only supports finite fields.")
    end

    lens, cfs, exps = convert_singular_ideal_to_array(J)

    gb_ld  = Ref(Cint(0))
    gb_len = Ref(Ptr{Cint}(0))
    gb_exp = Ref(Ptr{Cint}(0))
    gb_cf  = Ref(Ptr{Cvoid}(0))

    nr_terms  = ccall((:f4_julia, libneogb), Int,
        (Ptr{Nothing}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cvoid}},
        Ptr{Cint}, Ptr{Cint}, Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cint,
        Cint, Cint, Cint, Cint, Cint, Cint, Cint),
        cglobal(:jl_malloc), gb_ld, gb_len, gb_exp, gb_cf, lens, exps, cfs,
        field_char, mon_order, elim_block_size, nr_vars, nr_gens, initial_hts,
        nr_thrds, max_nr_pairs, 0, la_option, reduce_gb, 0, info_level)

    if nr_terms == 0
        error("Something went wrong in the C code of F4.")
    end
    # convert to julia array, also give memory management to julia
    jl_ld   = gb_ld[]
    jl_len  = Base.unsafe_wrap(Array, gb_len[], jl_ld)
    jl_exp  = Base.unsafe_wrap(Array, gb_exp[], nr_terms*nr_vars)
    if 0 == field_char
        ptr   = reinterpret(Ptr{BigInt}, gb_cf[])
        jl_cf = Base.unsafe_wrap(Array, ptr, nr_terms)
    elseif isprime(field_char)
        ptr   = reinterpret(Ptr{Int32}, gb_cf[])
        jl_cf = Base.unsafe_wrap(Array, ptr, nr_terms)
    end

    # construct Singular ideal
    # if 0 == field_char
    #     basis = convert_qq_gb_array_to_singular_ideal(
    #       jl_ld, jl_len, jl_cf, jl_exp, R)
    # else
        basis = convert_ff_gb_array_to_singular_ideal(
          jl_ld, jl_len, jl_cf, jl_exp, R)
    # end
    ccall((:free_f4_julia_result_data, libneogb), Nothing ,
          (Ptr{Nothing}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
           Ptr{Ptr{Cvoid}}, Int, Int),
          cglobal(:jl_free), gb_len, gb_exp, gb_cf, jl_ld, field_char)

    I.gb        = BiPolyArray(base_ring(I), basis)
    I.gb.isGB   = true
    I.gb.S.isGB = true
    
    return I.gb
end
