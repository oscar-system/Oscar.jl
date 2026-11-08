@testset "Printing" begin
  old_flag = Oscar.is_unicode_allowed()

  @test Oscar.is_unicode_allowed() == false
  @test allow_unicode(true; temporary=true) == false
  @test Oscar.is_unicode_allowed() == true

  struct AtoB
  end

  function Base.show(io::IO, ::AtoB)
    if Oscar.is_unicode_allowed()
      print(io, "A→B")
    else
      print(io, "A->B")
    end
  end

  allow_unicode(false; temporary=true)
  @test sprint(show, AtoB()) == "A->B"
  allow_unicode(true; temporary=true)
  @test sprint(show, AtoB()) == "A→B"

  # Restore old flag
  allow_unicode(old_flag; temporary=true)
end

# Test statistics metadata generation with and without Git information.
@testset "Test statistics metadata" begin
  metadata = Oscar._test_stats_metadata(Dict{Symbol,String}())
  @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}$", metadata.timestamp)
  @test metadata.commit == "v$(Oscar.VERSION_NUMBER)"

  git_info = Dict(
    :commit => "0123456789abcdef",
    :date => "2026-07-24 12:17:12 +0200",
  )
  metadata = Oscar._test_stats_metadata(git_info)
  @test metadata == (timestamp="2026-07-24T12-17-12", commit="0123456")
end
