@testset "Printing" begin
  if VERSION <= v"1.8.0"
    old_flag = Oscar.is_unicode_allowed()

    @test Oscar.is_unicode_allowed() == false
    @test allow_unicode(true) == false
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

    allow_unicode(false)
    @test sprint(show, AtoB()) == "A->B"
    allow_unicode(true)
    @test sprint(show, AtoB()) == "A→B"

    # Restore old flag
    allow_unicode(old_flag)
  end
end
