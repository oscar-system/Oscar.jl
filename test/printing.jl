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
