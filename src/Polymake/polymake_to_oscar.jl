## conversions of Polymake objects to Oscar objects

function convert(p::Polymake.Polynomial{Polymake.Rational, Polymake.to_cxx_type(Int64)};
                          parent::Union{MPolyRing, Nothing} = nothing)
    coeff_vec = convert(Vector{QQFieldElem}, Polymake.coefficients_as_vector(p))
    monomials = Matrix{Int}(Polymake.monomials_as_matrix(p))
    n_vars = ncols(monomials)

    # variable labeling from polymake is shifted
    if isnothing(parent)
        parent, _ = polynomial_ring(QQ, :x => 1:n_vars, cached=false)
    end

    return parent(coeff_vec, [monomials[i, :] for i in 1:nrows(monomials)]) 
end

function convert(::Type{MPolyIdeal{QQMPolyRingElem}}, O::Polymake.BigObject)
    n_vars = O.N_VARIABLES
    R, _ = polynomial_ring(QQ, :x => 1:n_vars, cached=false)
    converted_generators = map(p -> convert(p, parent=R), O.GENERATORS)
    
    return ideal(R, converted_generators)
end

_pmdata_for_oscar(::Nothing, coeff::Field) = nothing
_pmdata_for_oscar(v::Union{Bool,Int64,Float64,String}, coeff::Field) = v
if Polymake.CxxWrap.CxxLong != Int64
  _pmdata_for_oscar(i::Polymake.CxxWrap.CxxLong, coeff::Field) = Int64(i)
end

_pmdata_for_oscar(im::IncidenceMatrix, coeff::Field) = im

_pmdata_for_oscar(g::Polymake.Graph{T}, coeff::Field) where T = Graph{T}(g)

_pmdata_for_oscar(m::Polymake.Matrix, coeff::Field) = matrix(coeff, m)
_pmdata_for_oscar(m::Polymake.Matrix{Polymake.to_cxx_type(Int)}, coeff::Field) = Matrix{Int}(m)
_pmdata_for_oscar(m::Polymake.Matrix{<:Polymake.Integer}, coeff::Field) = matrix(ZZ, m)
_pmdata_for_oscar(m::Polymake.Matrix{<:Polymake.Rational}, coeff::Field) = matrix(QQ, m)

_pmdata_for_oscar(m::Polymake.SparseMatrix, coeff::Field) = _pmdata_for_oscar(Polymake.common.dense(m), coeff)

_pmdata_for_oscar(v::Polymake.Vector, coeff::Field) = collect(elem_type(coeff), map(coeff, v))
_pmdata_for_oscar(v::Polymake.Vector{Polymake.to_cxx_type(Int)}, coeff::Field) = Vector{Int}(v)
_pmdata_for_oscar(v::Polymake.Vector{<:Polymake.Integer}, coeff::Field) = collect(ZZRingElem, map(ZZ, v))
_pmdata_for_oscar(v::Polymake.Vector{<:Polymake.Rational}, coeff::Field) = collect(QQFieldElem, map(QQ, v))

_pmdata_for_oscar(v::Polymake.SparseVector, coeff::Field) = _pmdata_for_oscar(Polymake.common.dense(v), coeff)

_pmdata_for_oscar(nm::Polymake.NodeMap, coeff::Field) = _pmdata_for_oscar(Polymake.Array(nm), coeff)
_pmdata_for_oscar(bd::Polymake.BasicDecoration, coeff::Field) = (_pmdata_for_oscar(Polymake.decoration_face(bd), coeff), _pmdata_for_oscar(Polymake.decoration_rank(bd), coeff))

_pmdata_for_oscar(s::Polymake.Integer, coeff::Field) = ZZ(s)
_pmdata_for_oscar(s::Polymake.Rational, coeff::Field) = QQ(s)
_pmdata_for_oscar(s::Polymake.OscarNumber, coeff::Field) = coeff(s)

_convert_pm_minormax(::Type{Polymake.Max}) = typeof(max)
_convert_pm_minormax(::Type{Polymake.Min}) = typeof(min)
_pmdata_for_oscar(s::Polymake.TropicalNumber{A}, coeff::Field) where A = tropical_semiring(_convert_pm_minormax(A))(s)

_pmdata_for_oscar(s::Polymake.CxxWrap.StdString, coeff::Field) = String(s)

_pmdata_for_oscar(a::Polymake.Array, coeff::Field) = [_pmdata_for_oscar(e, coeff) for e in a]
_pmdata_for_oscar(a::Polymake.Array{T}, coeff::Field) where T <: Union{Polymake.Matrix, Polymake.Vector} = Tuple(_pmdata_for_oscar.(a, Ref(coeff)))

_pmdata_for_oscar(s::Polymake.Set, coeff::Field) = Set(_pmdata_for_oscar(e, coeff) for e in s)
