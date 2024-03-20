
########################################################################
# Required methods for AbsCoveredScheme                                #
########################################################################
### Forwarding of the non-documented getters
Base.in(U::AbsAffineScheme, X::AbsCoveredScheme) = (U in underlying_scheme(X))::Bool
getindex(X::AbsCoveredScheme, i::Int) = coverings(underlying_scheme(X))[i]::Covering
getindex(X::AbsCoveredScheme, C::Covering, D::Covering) = getindex(underlying_scheme(X), C, D)::CoveringMorphism
#setindex(X::AbsCoveredScheme, f::CoveringMorphism, C::Covering, D::Covering) = setindex(underlying_scheme(X), f, C, D)::CoveringMorphism
gluings(X::AbsCoveredScheme) = gluings(underlying_scheme(X))::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:AbsGluing}

########################################################################
# Finding coverings, patches, and gluings                             #
########################################################################
getindex(X::CoveredScheme, C::CoveringType, D::CoveringType) where {CoveringType<:Covering} = X.refinements[(C, D)]
setindex(X::CoveredScheme, f::CoveringMorphismType, C::CoveringType, D::CoveringType) where {CoveringMorphismType<:CoveringMorphism, CoveringType<:Covering} = X.refinements[(C, D)]
getindex(X::CoveredScheme, i::Int) = coverings(X)[i]


function getindex(X::AbsCoveredScheme, C::Covering)
  for i in 1:length(coverings(X))
    C == coverings(X)[i] && return i
  end
  error("covering not listed")
end

getindex(X::AbsCoveredScheme, i::Int, j::Int) = X[X[i], X[j]]

Base.in(U::AbsAffineScheme, X::CoveredScheme) = any(C->(U in C), coverings(X))

########################################################################
# Printing                                                             #
########################################################################

# We use several trickes here to ensure a nice alignment of the information.
# First, we agreed on printing the detailed on the covering of the scheme - one
# may always be interested in the local information about the (covered) scheme.
#
# We label each patch so that in the case of many (>= 5) patches, it is easier
# to see to what patch we refer (and the ordering match with the ordering of
# `default_covering(X)` here).
#
# We first then create a list of strings for all the coordinates, for each
# chart, and look for the length of the longest string: this allows us to
# estimate an offset to add everything we print the coordinates. Hence, in front
# of each coordinate system, one prints the description of the charts and with a
# good handling of the offsets, the description should all be aligned on the
# lefts, 3 spaces after the longest set of coordinates.
#
# If we have more than 10 patches, we need to take care of a left offset, so
# that the label of the patches and aligned on the right: we use here the
# variables `l` and `li` which check the number of digits
function Base.show(io::IO, ::MIME"text/plain", X::AbsCoveredScheme)
  io = pretty(io)
  println(io, "Scheme")
  println(io, Indent(), "over ", Lowercase(), base_ring(X))
  C = default_covering(X)
  if length(C) > 0
    l = ndigits(length(C))
    println(io, Dedent(), "with default covering")
    print(io, Indent(), "described by patches")
    print(io, Indent())
      for i in 1:length(C)
      li = ndigits(i)
      println(io)
      print(io, " "^(l-li)*"$i: ", Lowercase(), C[i])
    end
    println(io)
    print(io, Dedent(), "in the coordinate(s)")
    print(io, Indent())
    for i in 1:length(C)
      li = ndigits(i)
      println(io)
      print(io, " "^(l-li)*"$i: ")
      co = coordinates(C[i])
      str = "["*join(co, ", ")*"]"
      print(io, str)
    end
    print(io, Dedent())
    print(io, Dedent())
  else
    print(io, "with empty default covering")
  end
end

_covering_for_printing(io::IO, X::AbsCoveredScheme) = get(io, :covering, get_attribute(X, :simplified_covering, default_covering(X)))

# The explanation of the printing are the same as the function above. We need
# this `_show_semi_compact` for detailed nested printing. Here, we assume that
# in our nest, we already made precised the base ring of our covered scheme.
# 
# We then only print detailed on the covering. Note that in some cases, we night
# not want to only print the default covering, but a simplified one or a very
# specific one (when one prints ideal sheaves, for instance). In that case, we
# need the flexibility on the covering, using the variable `cov`.
#
# When one has to deal with morphisms of covered schemes, we want to distinguish
# charts of the domain and chart of the codomain. So, in addition to the
# numbering of the labels, we add the possibility to have more using the string
# `n`. Note that, in default use in CoveredSchemeMor printing, `n` is either "a"
# or "b" depending on whether `X` is the domain or the codomain.
function _show_semi_compact(io::IO, X::AbsCoveredScheme, cov::Covering, n::String)
  io = pretty(io)
  show(IOContext(io, :covering => cov, :show_semi_compact => false), X)
  print(io, Indent())
  co_str = String[""]
  for i in 1:length(cov)
    U = cov[i]
    co = coordinates(U)
    str = "["*join(co, ", ")*"]"
    push!(co_str, str)
  end
  k = max(length.(co_str)...)
  for i in 1:length(cov)
    l = ndigits(length(cov))
    li = ndigits(i)
    U = cov[i]
    kc = length(co_str[i+1])
    println(io)
    print(io, " "^(l-li)*"$(i)"*n*": "*co_str[i+1]*" "^(k-kc+3), Lowercase(), U)
  end
  print(io, Dedent())
end

function Base.show(io::IO, X::AbsCoveredScheme)
  cov = _covering_for_printing(io, X)
  if get(io, :show_semi_compact, false)
    l = get(io, :label, "")
    _show_semi_compact(io, X, cov, l)
  else
    io = pretty(io)
    n = n_patches(cov)
    if has_name(X)
      print(io, name(X))
    elseif get(io, :supercompact, false)
      print(io, "Covered scheme")
    else
      if get_attribute(X, :is_empty, false) || n == 0
        print(io, "Empty covered scheme over ")
      else
        print(io, "Scheme over ")
      end
      print(IOContext(io, :supercompact => true), Lowercase(), base_ring(X))
    end
    n > 0 && print(io, " covered with ", ItemQuantity(n, "patch"))
  end
end

########################################################################
# Base change
########################################################################
function base_change(phi::Any, X::AbsCoveredScheme)
  C = default_covering(X)
  CC, f_CC = base_change(phi, C)
  XX = CoveredScheme(CC)
  return XX, CoveredSchemeMorphism(XX, X, f_CC)
end




"""
    normalization(X::AbsCoveredScheme; check::Bool=true) -> Vector{Tuple{CoveredScheme, CoveredSchemeMorphism}}

Return the normalization of the reduced scheme ``X``.

# Input:
- A reduced scheme ``X``
- if `check` is `true` confirm that ``X`` is reduced; this is expensive

# Output:
A list of pairs ``(Y_i, f_i)`` where ``Y_i`` is a normal scheme and
``f_i`` is a morphism from ``Y_i`` to ``X``.
The disjoint union of the ``Y_i`` is the normalization of ``X``
and the ``f_i`` are the restrictions of the normalization morphism to ``Y_i``.
"""
function normalization(X::AbsCoveredScheme; check::Bool=true)
  @check is_reduced(X) "The scheme X=$(X) needs to be reduced."
  irred_comps_sheaf = Oscar.maximal_associated_points(ideal_sheaf(X))
  inc_maps = [Oscar.CoveredClosedEmbedding(scheme(irred_comps_sheaf[i]),irred_comps_sheaf[i])
                for i in 1:length(irred_comps_sheaf)]
  irred_comps = [domain(inc_i) for inc_i in inc_maps]
  norm_pairs = [_normalization_integral(Y_i; check) for Y_i in irred_comps]
  ret_value = [( domain(norm_pairs[i][2]),
                 compose(norm_pairs[i][2],inc_maps[i]) )
                 for i in 1:length(norm_pairs)]
  return ret_value
end

function _normalization_integral(X::AbsCoveredScheme; check::Bool=true)
  @check is_integral(X) "The scheme X=$(X) needs to be integral."
  orig_cov = default_covering(X)
  if has_attribute(X, :simplified_covering)
    orig_cov = simplified_covering(X)
  end
  X_norm_cov, pr_cov = _normalization_integral(orig_cov; check)
  Y = CoveredScheme(X_norm_cov)
  pr = CoveredSchemeMorphism(Y, X, pr_cov)
  return Y, pr
end

function _normalization_integral(C::Covering; check::Bool=true)
  X_i_norm_outputs = [_normalization(U; check, algorithm=:isPrime) for U in patches(C) if !is_empty(U)]
  C_norm = Covering(first.(first.(X_i_norm_outputs)))
  n = n_patches(C_norm)
  for i in 1:n
    for j in i+1:n
      gluing_i_j = _normalization_integral(
        C[i,j],
        X_i_norm_outputs[i][1],
        X_i_norm_outputs[j][1]
      )
      add_gluing!(C_norm, gluing_i_j)
    end
  end
  second(xs) = xs[2]
  morphisms = second.(first.(X_i_norm_outputs))
  covering_morphism_dict = IdDict(zip(domain.(morphisms), morphisms))
  normalization_map = CoveringMorphism(C_norm, C, covering_morphism_dict)
  return C_norm, normalization_map
end

struct NormalizationIntegralGluingData
  original_gluing::AbsGluing
  norm_output1::Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:Map}
  norm_output2::Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:Map}
  check::Bool
end

function _normalization_integral(
    G::AbsGluing,
    X_1_norm_output, # _normalization(X_1)[1]
    X_2_norm_output; # _normalization(X_2)[1]
    check::Bool=true
  )
  data = NormalizationIntegralGluingData(G, X_1_norm_output, X_2_norm_output, check)
  X_1_norm = X_1_norm_output[1]
  X_2_norm = X_2_norm_output[1]
  return LazyGluing(X_1_norm, X_2_norm, _compute_normalization_integral, data)
end

# Warning: assume patches irreducible
function _compute_normalization_integral(
    data::NormalizationIntegralGluingData
  )
  # Initialize the variables
  G = data.original_gluing
  X_1_norm_output = data.norm_output1
  X_2_norm_output = data.norm_output2
  check = data.check

  A, B = patches(G)
  #@check is_integral(A) && is_integral(B) "schemes must be integral"
  (X_1_norm, F_1, hom_to_K_1) = X_1_norm_output
  (X_2_norm, F_2, hom_to_K_2) = X_2_norm_output
  X_1, X_2 = patches(G)
  U_1, U_2 = gluing_domains(G)
  g_1, g_2 = gluing_morphisms(G)
  codomain(F_1) === X_1 || error("Codomain of F_1=$(F_1) should be X_1=$(X_1).")
  codomain(F_2) === X_2 || error("Codomain of F_2=$(F_2) should be X_2=$(X_2).")

  # We assume X_i_norm, F_i, hom_to_K_i are defined as below
  # X_1_norm, F_1, hom_to_K_1 = _normalization_affine(X_1)[1]
  # X_2_norm, F_2, hom_to_K_2 = _normalization_affine(X_2)[1]

  U_1_norm = PrincipalOpenSubset(X_1_norm, pullback(F_1)(complement_equation(U_1)))
  U_2_norm = PrincipalOpenSubset(X_2_norm, pullback(F_2)(complement_equation(U_2)))

  U_1_norm_mor = restrict(F_1, U_1_norm, U_1; check=false)
  U_2_norm_mor = restrict(F_2, U_2_norm, U_2; check=false)

  hom_to_U_2_norm_coordinates = elem_type(OO(U_1_norm))[]
  for a in hom_to_K_2.(coordinates(X_2_norm))
    b_num = pullback(U_1_norm_mor)(pullback(g_1)(OO(U_2)(numerator(a))))
    b_denom = pullback(U_1_norm_mor)(pullback(g_1)(OO(U_2)(denominator(a))))
    # TODO remove check that b_denom is invertible
    push!(hom_to_U_2_norm_coordinates, b_num * inv(b_denom))
  end
  hom_to_U_2_norm = morphism(U_1_norm, U_2_norm, hom_to_U_2_norm_coordinates)

  hom_to_U_1_norm_coordinates = elem_type(OO(U_2_norm))[]
  for a in hom_to_K_1.(coordinates(X_1_norm))
    b_num = pullback(U_2_norm_mor)(pullback(g_2)(OO(U_1)(numerator(a))))
    b_denom = pullback(U_2_norm_mor)(pullback(g_2)(OO(U_1)(denominator(a))))
    # TODO remove check that b_denom is invertible
    push!(hom_to_U_1_norm_coordinates, b_num * inv(b_denom))
  end
  hom_to_U_1_norm = morphism(U_2_norm, U_1_norm, hom_to_U_1_norm_coordinates; check=false)

  return SimpleGluing(
    X_1_norm,
    X_2_norm,
    hom_to_U_2_norm,
    hom_to_U_1_norm;
    check=false
  )
end
