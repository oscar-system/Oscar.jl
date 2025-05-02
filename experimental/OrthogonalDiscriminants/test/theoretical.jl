@testset "group order" begin
  d = 8
  q = 3
  order_oplus = order(omega_group(1, d, q))
  order_ominus = order(omega_group(-1, d, q))
  f = Oscar.OrthogonalDiscriminants.order_omega_mod_N
  @test f(d, q, order_oplus) == (true, false)
  @test f(d, q, order_ominus) == (false, true)
  @test f(d, q, order(omega_group(0, d-1, q))) == (true, true)
  @test f(d, q, q*order_oplus) == (false, false)
  @test f(d, q, (q-1)*order_oplus) == (false, false)

  @test_throws ArgumentError f(5, 2, 1)

  for entry in all_od_infos(comment_matches => "order")
    chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
    @test Oscar.OrthogonalDiscriminants.od_from_order(chi) == (true, entry[:valuestring])
  end
  for entry in all_od_infos(identifier => "A8")
    if ! comment_matches(entry, "order")
      chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
      @test Oscar.OrthogonalDiscriminants.od_from_order(chi) == (false, "")
    end
  end
end

@testset "eigenvalues" begin
  for entry in all_od_infos(comment_matches => "ev")
    chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
    @test Oscar.OrthogonalDiscriminants.od_from_eigenvalues(chi) == (true, entry[:valuestring])
  end
  for entry in all_od_infos(identifier => "A8")
    if ! comment_matches(entry, "ev")
      chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
      @test Oscar.OrthogonalDiscriminants.od_from_eigenvalues(chi) == (false, "")
    end
  end
end

@testset "Specht modules" begin
  for entry in all_od_infos(comment_matches => "specht")
    chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
    @test Oscar.OrthogonalDiscriminants.od_for_specht_module(chi) == (true, entry[:valuestring])
  end
  for entry in all_od_infos(identifier => "A8")
    if ! comment_matches(entry, "specht")
      chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
      @test Oscar.OrthogonalDiscriminants.od_for_specht_module(chi) == (false, "")
    end
  end
end

@testset "p-groups" begin
  rx = r"syl\((?<p>\d+)\)"
  info = all_od_infos(characteristic => 0,
                 comment_matches => (x -> any(startswith("syl"), x)))
  for entry in info[1:100]
    chi = Oscar.OrthogonalDiscriminants.character_of_entry(entry)
    c = filter(startswith("syl"), entry[:comment])
    for e in c
      p = parse(Int, match(rx, e)[:p])
      (flag, valstring) = Oscar.OrthogonalDiscriminants.od_from_p_subgroup(chi, p)
      @test flag
      @test valstring == entry[:valuestring]
    end
  end
end
