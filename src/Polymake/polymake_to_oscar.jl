## conversions of Polymake objects to Oscar objects

function convert(p::Polymake.Polynomial{Polymake.Rational, Polymake.to_cxx_type(Int64)};
                          parent::Union{MPolyRing, Nothing} = nothing)
    coeff_vec = convert(Vector{fmpq}, Polymake.coefficients_as_vector(p))
    monomials = Matrix{Int}(Polymake.monomials_as_matrix(p))
    n_vars = length(monomials[:, 1])
    # not sure if the numbering is the best choice but it matches Polymake
    if isnothing(parent)
        parent, _ = PolynomialRing(QQ, "x" => 0:n_vars - 1, cached=false)
    end

    return parent(coeff_vec, [monomials[i, :] for i in 1:nrows(monomials)])
end

function convert(O::Polymake.BigObject)
    big_object_name = Polymake.type_name(O)

    if "Ideal" == big_object_name
        n_vars = O.N_VARIABLES
        R, _ = PolynomialRing(QQ, "x" => 0:n_vars - 1, cached=false)
        converted_generators = map(p -> convert(p, parent=R), O.GENERATORS)
        
        return ideal(R, converted_generators)
    else
        throw(MethodError)
    end
end

