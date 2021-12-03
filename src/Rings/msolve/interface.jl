# we take a Singular ideal and extract the following data:
# an int array lengths storing the lengths of each generator
# an int array cfs storing the coefficients of each generator
# an int array exps storing the exponent vectors of each generator
function convert_singular_ideal_to_array(id::Singular.sideal)
    ngens   = Singular.ngens(id)
    nvars   = Singular.nvars(base_ring(id))
    char    = Singular.characteristic(base_ring(id))
    lens    = Int32[Singular.length(id[i]) for i in 1:ngens]
    nterms  = sum(lens)
    exps    = Int32[]
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
        blen::Array{Int32,1},
        bcf::Array{Int32,1},
        bexp::Array{Int32,1},
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
#         blen::Array{Int32,1},
#         bcf::Ptr{T} where {T <: Signed},
#         bexp::Array{Int32,1},
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
