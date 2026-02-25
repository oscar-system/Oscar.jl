# Use real_solutions to compute roots of a univariate polynomial
function _real_roots(f::QQMPolyRingElem; selected_precision::Int)
  (result, _) = Oscar.real_solutions(
    Vector{QQFieldElem}, ideal(f); precision=selected_precision
  )
  result = [r[1] for r in result]
  sort!(result)
  return result
end

function _analyse_singularity_msolve(
  f::QQMPolyRingElem,
  dfx::QQMPolyRingElem,
  projy::Map,
  xpt::Vector{QQFieldElem},
  ypt::Vector{QQFieldElem};
  selected_precision::Int,
)
  Rxy = parent(f)
  y = gens(Rxy)[2]
  x = gens(Rxy)[1]
  type = :singularity
  if (xpt[1] == xpt[2])
    # Exact case
    fpy = derivative(f, y)
    fpx = derivative(f, x)
    issingular =
      is_zero(Oscar.evaluate(fpy, [xpt[1], ypt[1]])) &&
      is_zero(Oscar.evaluate(fpx, [xpt[1], ypt[1]]))
    if !issingular
      type = :ytangent
    end
    pts = _real_roots(projy(Oscar.evaluate(f, [Rxy(xpt[1]), y])); selected_precision)
    singindex = findall(p -> p == ypt[1], pts)
    @assert length(singindex) == 1 "Something went wrong for exact solution"
    return (pts, singindex[1], type)
  else
    issingular =
      iszero(
        (
          1.0 +
          Float64(constant_coefficient(Oscar.evaluate(dfx, [Rxy(xpt[1]), Rxy(ypt[1])])))
        ) - 1.0,
      ) &&
      iszero(
        (
          1.0 +
          Float64(constant_coefficient(Oscar.evaluate(dfx, [Rxy(xpt[2]), Rxy(ypt[1])])))
        ) - 1.0,
      ) &&
      iszero(
        (
          1.0 +
          Float64(constant_coefficient(Oscar.evaluate(dfx, [Rxy(xpt[1]), Rxy(ypt[2])])))
        ) - 1.0,
      ) &&
      iszero(
        (
          1.0 +
          Float64(constant_coefficient(Oscar.evaluate(dfx, [Rxy(xpt[2]), Rxy(ypt[2])])))
        ) - 1.0,
      )
    if !issingular
      type = :ytangent
    end
    xbefore = xpt[1]
    ptsbefore = _real_roots(projy(Oscar.evaluate(f, [Rxy(xbefore), y])); selected_precision)
    xafter = xpt[2]
    ptsafter = _real_roots(projy(Oscar.evaluate(f, [Rxy(xafter), y])); selected_precision)

    # Just guessing some precision, not optimal...
    yinterval = [ypt[1] - QQ(1,2)^64, ypt[2] + QQ(1,2)^64]
    result = ptsbefore
    diff = length(ptsbefore) - length(ptsafter)
    if diff < 0
      # Choose longer vector
      result = ptsafter
    end

    # Just verify that everything fits with float
    singy = Float64(result[findfirst(r -> yinterval[1] <= r && r <= yinterval[2], result)])
    singindices = Int[]
    outr = typeof(result)()
    singfound = -1
    for i in 1:length(result)
      pt = result[i]
      if yinterval[1] <= pt && pt <= yinterval[2]
        if singfound == -1
          push!(outr, pt)
          singfound = i
        end
        if Float64(pt) != singy
          type = :broken
        end
        push!(singindices, i)
      else
        push!(outr, pt)
      end
    end

    if type == :ytangent
      ls = length(singindices)
      if ls != 2
        type = :broken
      end
    end

    return (outr, singfound, type)
  end
end

function msolve_sings(I; info_level::Int=0, precision::Int=128)
  # Just to remember how to get more debug info
  # (sings, rs2) = real_solutions(ideal(f) + ideal(derivative(f, y)); info_level=2, precision=selected_precision, interval=true)
  (sings, _) = Oscar.real_solutions(Vector{Vector{QQFieldElem}}, I; info_level, precision)
  sort!(sings; by=first)
  return sings
end

function isotopy_graph_from_msolve(
  IG::_IsotopyGraph, f_in::QQMPolyRingElem, random_transform::QQMatrix,
  selected_precision::Int,
)::Bool
  Rxy = parent(f_in)
  @assert nvars(Rxy) == 2 "Need curve in affine plane"

  # This matrix is a generic rotation by definition
  f = f_in((random_transform * gens(Rxy))...)

  (x, y) = gens(Rxy)
  dfx = derivative(f, x)
  K = base_ring(Rxy)
  Ry, t = polynomial_ring(K, [:y])
  projy = hom(Rxy, Ry, [0, t[1]])
  sings = msolve_sings(ideal([f, derivative(f, y)]); precision=selected_precision)
  @vprintln :DrawingCurves 2  "There are $(length(sings)) points to investigate"
  if(length(sings) == 0)
    throw(NotImplementedError(:isotopy_graph, "Curve has no critical points, it may be definite, this situation is not implemented yet."))
  end
  svecs = []
  sindices = []
  stypes = Symbol[]

  # Compute points above and below the singularities.
  i = 0
  for (xpt, ypt) in sings
    i += 1
    (svec, si, type) = _analyse_singularity_msolve(
      f, dfx, projy, xpt, ypt; selected_precision
    )
    @vprintln :DrawingCurves 2 "$i $(Float64(xpt[1]))"
    @vprintln :DrawingCurves 2 "$i Was this singular? $type"
    if type == :broken
      return false
    end
    push!(svecs, svec)
    push!(sindices, si)
    push!(stypes, type)
  end

  singcoords = [[sings[i][1][1], sings[i][2][1]] for i in 1:length(svecs)]
  _assemble_isotopy_graph!(
    IG,
    f,
    singcoords,
    stypes,
    svecs,
    sindices,
    random_transform,
    selected_precision,
    _real_roots,
  )

  return true
end
