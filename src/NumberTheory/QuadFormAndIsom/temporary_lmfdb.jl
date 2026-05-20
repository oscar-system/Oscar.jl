# There is a compatibility issue with the LMFDB connection, so one has to dev
# it for what follows:
# - dev https://github.com/thofma/LMFDB.jl
# - in line 10 of `src/NumberField.jl`, add the signature
#   `db::LMFDB.LMFDBLite.LMFDBConnection` to the first input of the function
#   `Hecke.number_field(db, label::String).
# - open julia, and run `using Oscar, LMFDB`
# - fetch the database: `db = LMFDB.LMFDBLite.LMFDBConnection()`
# - load some genera: z.B. `l6=LMFDB.genera(db; rank = 6, nplus=6, det = >=(1) & <=(div(300, 8, RoundUp)));`
# - call the functions needed below

# To use this file: include it in the branch of the PR

function _primitive_embeddings_lmfdb(
  G1s::Vector{ZZGenus},
  G2s::Vector{ZZGenus};
  _save_path::Union{Nothing, String}=nothing,
)
  results = Tuple{String, String, ZZMatrix}[]
  for G2 in G2s
    append!(results, _primitive_embeddings_lmfdb(G1s, G2; _save_path))
    GC.gc()
  end
  return results
end

function _primitive_embeddings_lmfdb(
  G1s::Vector{ZZGenus},
  G2::ZZGenus;
  _save_path::Union{Nothing, String}=nothing,
)
  results = Tuple{String, String, ZZMatrix}[]
  G1s = filter(>(rank(G2))∘rank, G1s)
  filter!(G1 -> reduce(&, signature_pair(G1).>=signature_pair(G2)), G1s)
  if !is_even(G2)
    filter!(!is_even, G1s)
  end
  q2 = discriminant_group(G2)
  Ctx = ZZLatGluingCtx()
  push!(Ctx.modules, q2)

  for G1 in G1s
    if is_unimodular(G1)
      append!(results, _primitive_embeddings_in_unimodular_safe_lmfdb(G1, G2; _save_path))
    elseif isone(gcd(numerator(det(G1)), numerator(det(G2))))
      append!(results, _primitive_embeddings_coprime_det_safe_lmfdb(G1, G2; _save_path))
    else
      append!(results, _primitive_embeddings_generic_safe_lmfdb(G1, G2; Ctx, vi=(-1, 1), q2, _save_path))
    end
  end
  return results
end

function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, String, ZZMatrix},
)
  m = parse(Int, first(x[1]))
  _path = joinpath(_save_path, x[2], "prim_emb_$(m).txt")
  isfile(_path) || touch(_path)
  str = x[1]*"_["
  v = x[3]
  for j in 1:ncols(v)
    k = v[1,j]
    str *= "$k"
    if j != ncols(v)
      str *= ","
    end
  end
  str *= "]\n"
  f = open(_path, "a")
  write(f, str)
  close(f)
end

function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, String, Hecke.IntegerUnion, Vector{QQFieldElem}},
)
  _path = joinpath(_save_path, x[1], "overlat.txt")
  isfile(_path) || touch(_path)
  str = x[2]*"_$(x[3])_["
  v = x[4]
  for j in 1:length(v)
    k = v[j]
    str *= "$(numerator(k))"
    if j != length(v)
      str *= ","
    end
  end
  str *= "]\n"
  f = open(_path, "a")
  write(f, str)
  close(f)
end

function save_data_lmfdb(
  _save_path::String,
  x::Tuple{String, Int},
)
  _path = joinpath(_save_path, x[1], "k3_embed.txt")
  isfile(_path) || touch(_path)
  f = open(_path, "w")
  write(f, "$(x[2])")
  close(f)
end

function _primitive_embeddings_in_unimodular_safe_lmfdb(
  G1::ZZGenus,
  G2::ZZGenus;
  _save_path::Union{Nothing, String}=nothing,
)
  results = Tuple{String, String, ZZMatrix}[]
  if is_even(G1)
    parity = :even
    par = :even
  else
    parity = :odd
    par = :both
  end
  sign = signature_pair(G1) .- signature_pair(G2)
  q = rescale(discriminant_group(G2), -1; cached=false)
  GKs = _integer_genera(q, sign, par)
  isempty(GKs) && return results

  Ns = ZZLat[]
  for GK in GKs
    append!(Ns, representatives(GK))
  end
  Ms = get_attribute(G2, :representatives)
  for M in Ms
    label_bottom = get_attribute(M, :lmfdb_label)
    for N in Ns
      _, tmp = unimodular_primitive_extensions(M, N; parity)
      for (T, M2, N2) in tmp
        @assert gram_matrix(M) == gram_matrix(M2)
        genus(T) != G1 && continue #Should never happen but still
        for T2 in get_attribute(G1, :representatives)
          if is_definite(G1)
            ok, f_top = is_isometric_with_isometry(T, T2)
          elseif length(get_attribute(G1, :representatives)) == 1
            ok = true
          else
            ok = is_isometric(T, T2)
          end
          !ok && continue
          label_top = get_attribute(T2, :lmfdb_label)
          if !is_definite(G1)
            v = zero_matrix(ZZ, 0, 0)
          else
            v = map_entries(ZZ, coordinates(basis_matrix(N2), T)*f_top)
          end
          __x = (label_top, label_bottom, v)
          @show __x
          push!(results, __x)
          if !isnothing(_save_path)
            save_data_lmfdb(_save_path, __x)
          end
        end
      end
    end
  end
  return results
end

### Coprime det

function _primitive_embeddings_coprime_det_safe_lmfdb(
  G1::ZZGenus,
  G2::ZZGenus;
  _save_path::Union{Nothing, String}=nothing,
)
  results = Tuple{String, String, ZZMatrix}[]
  if is_even(G1)
    parity = :even
    par = :even
  else
    parity = :odd
    par = :both
  end
  sign = signature_pair(G1) .- signature_pair(G2)
  q1 = discriminant_group(G1)
  q2 = discriminant_group(G2)
  q, _ = direct_sum(q1, rescale(q2, -1; cached=false); cached=false, as_bilinear_module=(par !== :even))
  GKs = _integer_genera(q, sign, par)
  isempty(GKs) && return results

  Ns = ZZLat[]
  for GK in GKs
    append!(Ns, representatives(GK))
  end
  Ms = get_attribute(G2, :representatives)
  for M in Ms
    GM, _ = image_in_Oq(M)
    label_bottom = get_attribute(M, :lmfdb_label)
    for N in Ns
      GN, _ = image_in_Oq(N)
      _, tmp = _primitive_extensions_coprime_left(M, N, parity, GM, GN)
      for (T, M2, N2) in tmp
        @assert gram_matrix(M) == gram_matrix(M2)
        genus(T) != G1 && continue #Should never happen but still
        for T2 in get_attribute(G1, :representatives)
          if is_definite(G1)
            ok, f_top = is_isometric_with_isometry(T, T2)
          elseif length(get_attribute(G1, :representatives)) == 1
            ok = true
          else
            ok = is_isometric(T, T2)
          end
          !ok && continue
          label_top = get_attribute(T2, :lmfdb_label)
          if !is_definite(G1)
            v = zero_matrix(ZZ, 0, 0)
          else
            v = map_entries(ZZ, coordinates(basis_matrix(N2), T)*f_top)
          end
          __x = (label_top, label_bottom, v)
          @show __x
          push!(results, __x)
          if !isnothing(_save_path)
            save_data_lmfdb(_save_path, __x)
          end
        end
      end
    end
  end
  return results
end

function _primitive_embeddings_generic_safe_lmfdb(
  G1::ZZGenus,
  G2::ZZGenus;
  _save_path::Union{Nothing, String}=nothing,
  Ctx=nothing,
  vi::NTuple{2, Int}=(-1, -1),
  q2::TorQuadModule=discriminant_group(G2),
  Ms::Vector{ZZLat}=ZZLat[],
)
  results = Tuple{String, String, ZZMatrix}[]
  R = rescale(representative(G1), -1; cached=false)
  U = hyperbolic_plane_lattice()
  if is_even(G1)
    parity = :even
    T, _ = direct_sum(R, U; cached=false)
  else
    parity = :both
    T, _ = direct_sum(R, U, U; cached=false)
  end
  q1n = discriminant_group(T)
  signK = signature_pair(G1) .- signature_pair(G2)

  Fac, lgm = _local_glue_maps(q2, q1n, parity; Ctx, vi=reverse(vi))
  isempty(lgm) && return results
  DKs = Dict{ZZGenus, Vector{ZZLat}}()

  for x in lgm
    qx = _form_over(x, parity)
    qK = rescale(qx, -1; cached=false)
    GKs = _integer_genera(qK, signK, parity)
    isempty(GKs) && continue
    # At that point, we know that we have an extension
    for GK in GKs
      haskey(DKs, GK) && continue
      DKs[GK] = representatives(GK)
    end
    isempty(Ms) && append!(Ms, get_attribute(G2, :representatives))
    for M in Ms
      label_bottom = get_attribute(M, :lmfdb_label)
      qM = discriminant_group(M)
      D = _gluing_ambient(qM, q1n, parity)
      GM, _ = image_in_Oq(M)
      _ok, phiM = is_isometric_with_isometry(qM, q2)
      @assert _ok
      xM = _pullback_left(x, phiM, parity)
      xMs = _split_orbit_left_group(xM, GM, parity)
      for _x in xMs, y in _all_glue_maps(_x)
        V, M2, T2, GV = _overlattice_with_glue_stabilizer(y, D, parity)
        resV = NTuple{3, ZZLat}[]
        qV = domain(GV)
        for GK in GKs, K in DKs[GK]
          _, pe = unimodular_primitive_extensions(V, K; right_discriminant_action=GV, parity)
          append!(resV, pe)
        end
        for (S, V2, W2) in resV
          T3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(T2), zero_matrix(QQ, rank(T2), degree(W2)-degree(T2))))
          L = lll(orthogonal_submodule(S, T3))
          genus(L) == G1 || continue
          M3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(M2), zero_matrix(QQ, rank(M2), degree(W2)-degree(M2))))
          @assert is_sublattice(L, M3)
          @assert is_primitive(L, M3)
          N = orthogonal_submodule(L, M3)
          bM = coordinates(basis_matrix(M3), L)
          bN = coordinates(basis_matrix(N), L)
          L = integer_lattice(; gram=gram_matrix(L))
          M3 = lattice_in_same_ambient_space(L, bM)
          N = lattice_in_same_ambient_space(L, bN)
          for L2 in get_attribute(G1, :representatives)
            if is_definite(G1)
              ok, f_top = is_isometric_with_isometry(L, L2)
            elseif length(get_attribute(G1, :representatives)) == 1
              ok = true
            else
              ok = is_isometric(L, L2)
            end
            !ok && continue
            label_top = get_attribute(L2, :lmfdb_label)
            if is_definite(G1)
              v = map_entries(ZZ, coordinates(bN*f_top, L2))
            else
              v = zero_matrix(ZZ, 0, 0)
            end
            __x = (label_top, label_bottom, v)
            @show __x
            push!(results, __x)
            if !isnothing(_save_path)
              save_data_lmfdb(_save_path, __x)
            end
          end
        end
      end
    end
  end
  return results
end

function _prime_index_overlattices_lmfdb(
  Gs::Vector{ZZGenus};
  _save_path::Union{Nothing, String}=nothing,
)
  res = Tuple{String, String, Int, Vector{QQFieldElem}}[]
  for G in Gs
    q = discriminant_group(G)
    o = order(q)
    if isone(o)
      continue
    elseif radical(o) == o
      continue
    end
    if is_even(G)
      q = Hecke._as_finite_bilinear_module(q)
    end
    pds = prime_divisors(last(elementary_divisors(q)))
    D = Dict{eltype(pds), Tuple{TorQuadModule, Vector{Tuple{TorQuadModuleMap, GAPGroupHomomorphism}}}}()
    for p in pds
      qp, _ = torsion_subgroup(q, p)
      r = Oscar._stabilizer_isotropic_elementary(qp, 1)
      D[p] = (qp, r)
    end
    for M in get_attribute(G, :representatives)
      GM, _ = image_in_Oq(M)
      label_bottom = get_attribute(M, :lmfdb_label)
      T = domain(GM)
      if is_even(G)
        T = Hecke._as_finite_bilinear_module(T)
        GM = Oscar._orthogonal_group(T, matrix.(gens(GM)))
      end
      for p in pds
        qp, tmp = D[p]
        is_empty(tmp) && continue
        Tp, TpinT = torsion_subgroup(T, p)
        OTp = orthogonal_group_bilinear(Tp)
        iso = isomorphism(PermGroup, OTp)
        Gp, _ = restrict_automorphism_group(GM, TpinT)
        isoGp, _ = iso(Gp)
        ok, phi = is_isometric_with_isometry(Tp, qp)
        @assert ok
        iphi = inv(phi)
        for (Hinqp, j) in tmp
          _H, _HinTp = sub(Tp, elem_type(Tp)[phi\(qp(lift(a))) for a in gens(domain(Hinqp))])
          _, jTp = sub(codomain(iso), elem_type(codomain(iso))[iso(OTp(phi * hom(j(g)) * iphi; check=false)) for g in gens(domain(j))])
          for _h in double_cosets(codomain(iso), domain(jTp), isoGp)
            h = iso\(representative(_h))
            _Hh, _ = sub(Tp, [h(Tp(lift(a))) for a in gens(_H)])
            L = cover(_Hh)
            _v = p*lift(first(gens(_Hh)))
            G2 = genus(L)
            k = findfirst(isequal(G2), Gs)
            isnothing(k) && continue
            if length(get_attribute(Gs[k], :representatives)) == 1
              label_top = get_attribute(only(get_attribute(Gs[k], :representatives)), :lmfdb_label)
            else
              for L2 in get_attribute(Gs[k], :representatives)
                !is_isometric(L, L2) && continue
                label_top = get_attribute(L2, :lmfdb_label)
              end
            end
            __x = (label_bottom, label_top, p, _v)
            @show __x
            push!(res, __x)
            if !isnothing(_save_path)
              save_data_lmfdb(_save_path, __x)
            end
          end
        end
      end
    end
  end
  return res
end

function _embeddings_in_K3_lattice(
  Gs::Vector{ZZGenus};
  _save_path::Union{Nothing, String}=nothing,
)
  L = rescale(k3_lattice(), -1)
  set_attribute!(L, :lmfdb_label, "")
  GK3 = genus(L)
  set_attribute!(GK3, :representatives, ZZLat[L])
  res = Tuple{String, Int}[]
  for G in Gs
    !is_even(G) && continue
    p, n = signature_pair(G)
    n in [1,2] || continue
    r = _primitive_embeddings_in_unimodular_safe_lmfdb(GK3, G)
    u = unique!([rr[2] for rr in r])
    for s in u
      i = count(rr -> rr[2] == s, r)
      push!(res, (s, i))
      @show (s, i)
      if !isnothing(_save_path)
        save_data_lmfdb(_save_path, (s, i))
      end
    end
  end
  return res
end

# - `G1s` consists of all genera of lattices of rank n+1, with given
#   signatures and determinant in absolute value bounded by d/(n+1)
#   for some d.
# - `G2s` consists of all genera of lattice of rank n, with same
#   possible signatures and determinant in absolute value bounded by d/n.
# - `__save_path` is the absolute path to a folder where to save data (then
#   the dispatch is made automatically; see the code below).
function _computations_lmfdb(
  G1s::Vector{ZZGenus},
  G2s::Vector{ZZGenus};
  __save_path::Union{String, Nothing}=nothing,
)
  n = rank(first(G2s))
  if !isnothing(__save_path)
    _save_path = joinpath(__save_path, "rank$n")
    isdir(_save_path) || mkdir(_save_path)
    for G2 in G2s, L in get_attribute(G2, :representatives)
      label = get_attribute(L, :lmfdb_label)
      path_label = joinpath(_save_path, label)
      isdir(path_label) || mkdir(path_label)
    end
  else
    _save_path = nothing
  end

  overlat = _prime_index_overlattices_lmfdb(G2s; _save_path)
  k3_embed = _embeddings_in_K3_lattice(G2s; _save_path)
  prim_emb = _primitive_embeddings_lmfdb(G1s, G2s; _save_path)
  return prim_emb, overlat, k3_embed
end

