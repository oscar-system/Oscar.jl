# Use real_solutions to compute roots of a univariate polynomial
function _real_roots(f::QQMPolyRingElem; solver_precision::Int)
  (result, _) = real_solutions(
    Vector{QQFieldElem}, ideal(f); precision=solver_precision
  )
  result = [r[1] for r in result]
  sort!(result)
  return result
end

function _analyse_singularity_msolve(
  f::QQMPolyRingElem,
  projy::Map,
  xpt::Vector{QQFieldElem},
  ypt::Vector{QQFieldElem},
  type::Symbol,
  solver_precision::Int,
  ybox_tolerance::Int,
)
  Rxy = parent(f)
  x, y = gens(Rxy)
  if (xpt[1] == xpt[2])
    pts = _real_roots(projy(evaluate(f, [Rxy(xpt[1]), y])); solver_precision)
    singindex = findall(p -> ypt[1] <= p && p <= ypt[2], pts)
    @req length(singindex) == 1 "Something went wrong for exact solution"
    return (pts, singindex[1], :success)
  else
    xbefore = xpt[1]
    ptsbefore = _real_roots(projy(evaluate(f, [Rxy(xbefore), y])); solver_precision)
    xafter = xpt[2]
    ptsafter = _real_roots(projy(evaluate(f, [Rxy(xafter), y])); solver_precision)

    # Just guessing some precision, not optimal...
    ytol = QQ(1, 2)^ybox_tolerance
    yinterval = [ypt[1] - ytol, ypt[2] + ytol]
    result = ptsbefore
    diff = length(ptsbefore) - length(ptsafter)
    if diff < 0
      # Choose longer vector
      result = ptsafter
    end
    if type == :ytangent
      @req abs(diff) == 2 "Could not verify ytangent status"
    end

    # Just verify that everything fits with float
    pti = findfirst(r -> yinterval[1] <= r && r <= yinterval[2], result)

    if isnothing(pti)
      # We assume this is isolated
      if length(ptsbefore) != length(ptsafter)
        # Something is wrong if points before and after do not agree.
        return (result, 0, :fail)
      end
      coords = QQ(1, 2) * (ypt[1] + ypt[2])
      a = filter(r -> r < coords, result)
      b = filter(r -> r > coords, result)
      result = vcat(a, [coords], b)
      return (result, length(a) + 1, :success)
    end
    singy = Float64(result[pti])
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
        push!(singindices, i)
      else
        push!(outr, pt)
      end
    end
    if type == :ytangent
      if length(singindices) != 2
        return (result, 0, :fail)
      end
    end

    return (outr, singfound, :success)
  end
end

function msolve_wrapper(I; info_level::Int=0, precision::Int=128)
  # Just to remember how to get more debug info
  # (sings, rs2) = real_solutions(ideal(f) + ideal(derivative(f, y)); info_level=2, precision=solver_precision, interval=true)
  (sings, _) = real_solutions(Vector{Vector{QQFieldElem}}, I; info_level, precision)
  sort!(sings; by=first)
  return sings
end

_interval_contains(i1::QQFieldElem, i2::QQFieldElem, x::QQFieldElem) = i1 <= x && x <= i2

function _boxes_intersect(B1, B2)
  noxintersection = B1[1][1] > B2[1][2] || B2[1][1] > B1[1][2]
  noyintersection = B1[2][1] > B2[2][2] || B2[2][1] > B1[2][2]
  return !(noxintersection || noyintersection)
end

function isotopy_graph_from_msolve(
  IG::_IsotopyGraph, f_in::QQMPolyRingElem, random_transform::QQMatrix,
  solver_precision::Int, ybox_tolerance::Int,
)::Bool
  Rxy = parent(f_in)
  @req nvars(Rxy) == 2 "Need curve in affine plane"

  # This matrix is a generic rotation by definition
  f = f_in((random_transform * gens(Rxy))...)

  (x, y) = gens(Rxy)
  dfx = derivative(f, x)
  K = base_ring(Rxy)
  Ry, t = polynomial_ring(K, [:y]; cached=false)
  projy = hom(Rxy, Ry, [0, t[1]])
  critpt = []
  try
    critpt = msolve_wrapper(ideal([f, derivative(f, y)]); precision=solver_precision)
  catch Exception
    @vprintln :DrawingCurves 2 "Found no critical points."
    return false
  end
  @vprintln :DrawingCurves 2 "There are $(length(critpt)) critical points to investigate"
  for (a, b) in subsets(critpt, 2)
    # Check for genericity
    a1 = a[1][1]
    a2 = a[1][2]
    b1 = b[1][1]
    b2 = b[1][2]
    if !(b2 < a1 || a2 < b1)
      @vprintln :DrawingCurves 2 "Critical points were not in general position."
      return false
    end
  end
  # Check that x-axis intervals are disjoint
  for i in 1:(length(critpt) - 1)
    @req critpt[i][1][2] < critpt[i + 1][1][1] "Intervals must be disjoint. Try increasing precision."
  end

  if (length(critpt) == 0)
    throw(
      NotImplementedError(
        :isotopy_graph,
        "Curve has no critical points, it may be definite, this situation is not implemented yet.",
      ),
    )
  end
  singpts = msolve_wrapper(ideal([f]) + jacobian_ideal(f); precision=solver_precision)
  @vprintln :DrawingCurves 2 "There are $(length(singpts)) singular points."
  nsings = 0
  stypes = Symbol[]

  # Decide which points are ytangent and which are singular
  for cpt in critpt
    sct = 0
    for spt in singpts
      check = _boxes_intersect(spt, cpt)
      if check
        sct += 1
      end
    end
    @req sct <= 1 "Cannot uniquely verify singular pt"
    if sct == 1
      nsings += 1
      push!(stypes, :singularity)
    else
      push!(stypes, :ytangent)
    end
  end
  @req nsings == length(singpts) "Did not find all singularities"

  svecs = []
  sindices = []

  # Compute points above and below the singularities.
  i = 0
  for (xpt, ypt) in critpt
    i += 1
    type = stypes[i]
    (svec, si, state) = _analyse_singularity_msolve(
      f, projy, xpt, ypt, type, solver_precision, ybox_tolerance
    )
    if state == :fail
      @vprintln :DrawingCurves 2 "Failed analysing critical point of type $type."
      return false
    end
    @vprintln :DrawingCurves 2 "$i $(Float64(xpt[1]))"
    @vprintln :DrawingCurves 2 "$i Was this singular? $type"
    push!(svecs, svec)
    push!(sindices, si)
  end

  singcoords = [
    [(critpt[i][1][1] + critpt[i][1][2]) / 2, (critpt[i][2][1] + critpt[i][2][2]) / 2] for
    i in 1:length(svecs)
  ]
  _assemble_isotopy_graph!(
    IG,
    f,
    singcoords,
    stypes,
    svecs,
    sindices,
    random_transform,
    solver_precision,
    _real_roots,
  )

  return true
end
