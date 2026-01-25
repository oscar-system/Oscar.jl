@testset "Printings" begin
  function _show_details(io::IO, X::Union{ZZLatWithIsom, QuadSpaceWithIsom})
    return show(io, MIME"text/plain"(), X)
  end
  # Finite order
  L = root_lattice(:A, 2)
  Lf = integer_lattice_with_isometry(L)
  Vf = ambient_space(Lf)
  for X in [Lf, Vf]
    @test sprint(_show_details, X) isa String
    @test sprint(Oscar.to_oscar, X) isa String
    @test sprint(show, X) isa String
    @test sprint(show, X; context=:supercompact => true) isa String
  end
  # Infinite order
  L = integer_lattice(; gram = QQ[1 2; 2 1])
  f = QQ[4 -1; 1 0]
  Lf = integer_lattice_with_isometry(L, f)
  Vf = ambient_space(Lf)
  for X in [Lf, Vf]
    @test sprint(_show_details, X) isa String
    @test sprint(Oscar.to_oscar, X) isa String
    @test sprint(show, X) isa String
    @test sprint(show, X; context=:supercompact => true) isa String
  end
end
