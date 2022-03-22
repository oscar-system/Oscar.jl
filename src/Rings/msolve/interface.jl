# we take a Singular ideal and extract the following data:
# an int array lengths storing the lengths of each generator
# an int array cfs storing the coefficients of each generator
# an int array exps storing the exponent vectors of each generator
function convert_oscar_ideal_to_array(
        I::MPolyIdeal)

    R = base_ring(I)
    nr_vars     = nvars(R)
    nr_gens     = ngens(I)
    field_char  = characteristic(R)
    lens        = Int32[length(I[i]) for i in 1:nr_gens]
    nterms      = sum(lens)

    if field_char > 2^31
        error("At the moment f4 only supports finite fields up to prime characteristic < 2^31.")
    end
    # get coefficients
    if field_char == 0
        cfs = BigInt[]
    else
        cfs = Int32[]
    end
    if field_char == 0
        for i in 1:nr_gens
            for cf in coefficients(I[i])
                push!(cfs, BigInt(numerator(cf)))
                push!(cfs, BigInt(denominator(cf)))
            end
        end
    else
        for i in 1:nr_gens
            for cf in coefficients(I[i])
                push!(cfs, Int32(cf.data))
            end
        end
    end

    # get exponent vectors
    exps  = Int32[]
    for i in 1:nr_gens
        for ev in exponent_vectors(I[i])
            append!(exps, convert(Vector{Int32}, ev))
        end
    end

    return lens, cfs, exps
end

function convert_ff_gb_array_to_oscar_array(
        bld::Int32,
        blen::Vector{Int32},
        bcf::Vector{Int32},
        bexp::Vector{Int32},
        R::GFPMPolyRing
        )

    nr_gens = bld
    nr_vars = nvars(R)
    CR      = coefficient_ring(R)

    basis = Vector{gfp_mpoly}(undef, nr_gens)

    len   = 0

    for i in 1:nr_gens
        g  = MPolyBuildCtx(R)
        for j in 1:blen[i]
            push_term!(g, CR(bcf[len+j]),
                       convert(Vector{Int}, bexp[(len+j-1)*nr_vars+1:(len+j)*nr_vars]))
        end
        len +=  blen[i]
        basis[i]  = g.poly
    end

    return basis
end

function convert_singular_ideal_to_array(id::Singular.sideal)
    ngens   = Singular.ngens(id)
    nvars   = Singular.nvars(base_ring(id))
    char    = Singular.characteristic(base_ring(id))
    lens    = Int32[Singular.length(id[i]) for i in 1:ngens]
    nterms  = sum(lens)
    exps    = Int32[]
    if char > 2^31
        error("At the moment f4 only supports finite fields up to prime characteristic < 2^31.")
    end
    if char ==  0
        cfs = BigInt[]
    else
        cfs = Int32[]
    end
    for i = 1:Singular.ngens(id)
        if char == 0
            for c in Singular.coefficients(id[i])
                push!(cfs, Singular.libSingular.n_GetMPZ(numerator(c).ptr, Singular.ZZ.ptr))
                push!(cfs, Singular.libSingular.n_GetMPZ(denominator(c).ptr, Singular.ZZ.ptr))
            end
        else
            for c in Singular.coefficients(id[i])
                push!(cfs, Base.Int(c))
            end
        end
        for e in Singular.exponent_vectors(id[i])
            for j = 1:nvars
                push!(exps, Base.Int(e[j]))
            end
        end
    end
    return lens, cfs, exps
end

# we know that the terms are already sorted and they are all different
# coming from GroebnerBasis' F4 computation, so we do not need p_Add_q for
# the terms, but we can directly set the next pointers of the polynomials
function convert_ff_gb_array_to_singular_ideal(
        bld::Int32,
        blen::Vector{Int32},
        bcf::Vector{Int32},
        bexp::Vector{Int32},
        R::Singular.PolyRing
        )
    ngens = bld

    nvars = Singular.nvars(R)
    basis = Singular.Ideal(R, ) # empty ideal
    # first entry in exponent vector is module component => nvars+1
    exp   = zeros(Cint, nvars+1)

    list  = Singular.elem_type(R)[]
    # we generate the singular polynomials low level in order
    # to avoid overhead due to many exponent operations etc.
    j   = ngens + 1 + 1
    len = 1
    for i = 1:ngens
        # do the first term
        p = Singular.libSingular.p_Init(R.ptr)
        Singular.libSingular.pSetCoeff0(p, Clong(bcf[len]), R.ptr)
        for k = 1:nvars
            exp[k+1]  = bexp[(len-1) * nvars + k]
        end
        Singular.libSingular.p_SetExpV(p, exp, R.ptr)
        lp  = p
        for j = 2:blen[i]
          pterm = Singular.libSingular.p_Init(R.ptr)
          Singular.libSingular.pSetCoeff0(pterm, Clong(bcf[len+j-1]), R.ptr)
          for k = 1:nvars
              exp[k+1]  = bexp[(len+j-1-1) * nvars + k]
          end
          Singular.libSingular.p_SetExpV(pterm, exp, R.ptr)
          Singular.libSingular.SetpNext(lp, pterm)
          lp  = pterm
        end
        push!(list, R(p))
        len += blen[i]
    end
    return Singular.Ideal(R, list)
end

# function convert_qq_gb_array_to_singular_ideal(
#         bld::Int32,
#         blen::Vector{Int32},
#         bcf::Ptr{T} where {T <: Signed},
#         bexp::Vector{Int32},
#         R::Singular.PolyRing
#         )
#     ngens = bld
#
#     nvars = Singular.nvars(R)
#     basis = Singular.Ideal(R, ) # empty ideal
#     # first entry in exponent vector is module component => nvars+1
#     exp   = zeros(Cint, nvars+1)
#
#     list  = Singular.elem_type(R)[]
#     # we generate the singular polynomials low level in order
#     # to avoid overhead due to many exponent operations etc.
#     j   = ngens + 1 + 1
#     len = 1
#     for i = 1:ngens
#         # do the first term
#         p = Singular.libSingular.p_Init(R.ptr)
#         Singular.libSingular.p_SetCoeff0(p,
#                 Singular.libSingular.n_InitMPZ(BigInt(unsafe_load(bcf, len)),
#                     Singular.QQ.ptr), R.ptr)
#         for k = 1:nvars
#             exp[k+1]  = bexp[(len-1) * nvars + k]
#         end
#         Singular.libSingular.p_SetExpV(p, exp, R.ptr)
#         lp  = p
#         for j = 2:blen[i]
#           pterm = Singular.libSingular.p_Init(R.ptr)
#         Singular.libSingular.p_SetCoeff0(pterm,
#                 Singular.libSingular.n_InitMPZ(BigInt(unsafe_load(bcf, len+j-1)),
#                     Singular.QQ.ptr), R.ptr)
#           for k = 1:nvars
#               exp[k+1]  = bexp[(len+j-1-1) * nvars + k]
#           end
#           Singular.libSingular.p_SetExpV(pterm, exp, R.ptr)
#           Singular.libSingular.SetpNext(lp, pterm)
#           lp  = pterm
#         end
#         push!(list, R(p))
#         len += blen[i]
#     end
#     return Singular.Ideal(R, list)
# end
