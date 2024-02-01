@testset "the infinite zero double complex" begin
  # Setting up the factory for the entries
  mutable struct ZeroModuleFactory{ChainType} <: Oscar.ChainFactory{ChainType}
    R::MPolyRing

    function ZeroModuleFactory(R::MPolyRing)
      return new{ModuleFP{elem_type(R)}}(R)
    end
  end

  # Overwriting the call syntax as required
  function (fac::ZeroModuleFactory)(D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    return FreeMod(fac.R, 0)
  end

  function Oscar.can_compute(fac::ZeroModuleFactory, D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    return true
  end

  # Setting up the factories for the boundary maps
  mutable struct VerticalZeroMaps{MorphismType} <: Oscar.ChainMorphismFactory{MorphismType}
    R::MPolyRing

    function VerticalZeroMaps(R::MPolyRing)
      return new{ModuleFPHom}(R)
    end
  end

  mutable struct HorizontalZeroMaps{MorphismType} <: Oscar.ChainMorphismFactory{MorphismType}
    R::MPolyRing

    function HorizontalZeroMaps(R::MPolyRing)
      return new{ModuleFPHom}(R)
    end
  end

  # Overwriting their call syntax as required
  function (fac::VerticalZeroMaps)(D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    dom = D[i, j]
    inc = (Oscar.vertical_direction(D) == :chain ? -1 : 1)
    cod = D[i, j + inc]
    return hom(dom, cod, elem_type(cod)[])
  end

  function (fac::HorizontalZeroMaps)(D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    dom = D[i, j]
    inc = (Oscar.horizontal_direction(D) == :chain ? -1 : 1)
    cod = D[i + inc, j]
    return hom(dom, cod, elem_type(cod)[])
  end

  function Oscar.can_compute(fac::VerticalZeroMaps, D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    return true
  end

  function Oscar.can_compute(fac::HorizontalZeroMaps, D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    return true
  end

  # Implement the actual constructor
  function zero_double_complex(R::MPolyRing)
    entry_fac = ZeroModuleFactory(R)
    vert_map_fac = VerticalZeroMaps(R)
    horz_map_fac = HorizontalZeroMaps(R)

    result = Oscar.DoubleComplexOfMorphisms(entry_fac, horz_map_fac, vert_map_fac, horizontal_direction=:chain, vertical_direction=:chain)
    return result
  end

  # Test code for the above implementation
  R, (x, y, z) = QQ[:x, :y, :z]
  Z = zero_double_complex(R)

  @test iszero(Z[10, 20])
  M = Z[0, 0]
  a = zero(M)
  @test parent(a) === Z[0, 0] 
  @test iszero(Oscar.horizontal_map(Z, 20, 25))
  phi = Oscar.horizontal_map(Z, 11, 20)
  @test codomain(phi) === Z[10, 20]
  psi = Oscar.vertical_map(Z, 10, 21)
  @test codomain(psi) === Z[10, 20]
end

@testset "infinite one-line double complexes" begin
  # Set up a factory for the chains.
  struct MyNewChainFactory{ChainType} <: Oscar.ChainFactory{ChainType} 
    original_complex::ComplexOfMorphisms

    function MyNewChainFactory(C::ComplexOfMorphisms{T}) where {T<:ModuleFP}
      return new{ModuleFP}(C)
    end
  end

  function (fac::MyNewChainFactory)(D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    inc = (Oscar.typ(fac.original_complex) == :chain ? -1 : 1)
    lb = (Oscar.typ(fac.original_complex) == :chain ? last(range(fac.original_complex)) : first(range(fac.original_complex)))
    rb = (Oscar.typ(fac.original_complex) == :cochain ? last(range(fac.original_complex)) : first(range(fac.original_complex)))
    R = base_ring(fac.original_complex[lb])

    iszero(j) || error("invalid production request")

    # Extend by zeroes to both directions
    (i < lb || i > rb) && return FreeMod(R, 0)

    # Give back the entry of the original complex wherever applicable
    return fac.original_complex[i]
  end

  function Oscar.can_compute(fac::MyNewChainFactory, D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    return iszero(j)
  end

  # Set up a factory for the morphisms; we only need to worry about the horizontal 
  # ones as the vertical ones will be forbidden to ask for by restrictions on 
  # the bounds.
  struct MyNewMorphismFactory{MorphismType} <: Oscar.ChainMorphismFactory{MorphismType}
    original_complex::ComplexOfMorphisms

    function MyNewMorphismFactory(C::ComplexOfMorphisms{T}) where {T<:ModuleFP}
      return new{ModuleFPHom}(C)
    end
  end

  function (fac::MyNewMorphismFactory)(D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    inc = (Oscar.typ(fac.original_complex) == :chain ? -1 : 1)
    lb = (Oscar.typ(fac.original_complex) == :chain ? last(range(fac.original_complex)) : first(range(fac.original_complex)))
    rb = (Oscar.typ(fac.original_complex) == :cochain ? last(range(fac.original_complex)) : first(range(fac.original_complex)))
    R = base_ring(fac.original_complex[lb])

    iszero(j) || error("invalid production request")

    # Handle the trivial out of bounds case with zero maps
    (i < lb || i > rb) && return hom(D[i, j], D[i + inc, j], elem_type(D[i + inc, j])[zero(D[i + inc, j]) for k in 1:ngens(D[i, j])])

    # Handle the boundary cases
    (i == lb && inc == -1) && return hom(D[i, j], D[i + inc, j], elem_type(D[i + inc, j])[zero(D[i + inc, j]) for k in 1:ngens(D[i, j])])
    (i == rb && inc == 1) && return hom(D[i, j], D[i + inc, j], elem_type(D[i + inc, j])[zero(D[i + inc, j]) for k in 1:ngens(D[i, j])])

    # return the morphisms from the original complex otherwise
    return map(fac.original_complex, i)
  end

  function Oscar.can_compute(fac::MyNewMorphismFactory, D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int)
    return iszero(j)
  end

  struct MyDummyMorphismFactory{MorphismType} <: Oscar.ChainMorphismFactory{MorphismType}
    function MyDummyMorphismFactory(C::ComplexOfMorphisms{T}) where {T<:ModuleFP}
      return new{ModuleFPHom}()
    end
  end

  Oscar.can_compute(fac::MyDummyMorphismFactory, D::Oscar.AbsDoubleComplexOfMorphisms, i::Int, j::Int) = false

  # The actual constructor for the object we want to have.
  function as_infinite_one_line_double_complex(C::ComplexOfMorphisms{T}) where {T<:ModuleFP}
    is_complete(C) || error("implemented only for complete complexes")
    chain_fac = MyNewChainFactory(C)
    mor_fac = MyNewMorphismFactory(C)
    dummy_fac = MyDummyMorphismFactory(C)

    lb = (Oscar.typ(C) == :chain ? last(range(C)) : first(range(C)))
    rb = (Oscar.typ(C) == :cochain ? last(range(C)) : first(range(C)))
    return Oscar.DoubleComplexOfMorphisms(chain_fac, mor_fac, dummy_fac, 
                                          right_bound=rb, left_bound=lb, 
                                          horizontal_direction=Oscar.typ(C), vertical_direction=Oscar.typ(C),
                                          upper_bound=0, lower_bound=0
                                         )
  end

  # Test code
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, [x, y])
  F = FreeMod(R, 1)
  IF, _ = I*F
  M, _ = quo(F, IF)
  res = free_resolution(M).C # for the moment we need to access the field to get a `ComplexOfMorphisms`
  dc = as_infinite_one_line_double_complex(res)

  @test dc[-1, 0] === M # Should reproduce M
  @test dc[0, 0] isa FreeMod && rank(dc[0, 0]) == 1  # A free module over R of rank 1
  @test dc[1, 0] isa FreeMod && rank(dc[1, 0]) == 2  # A free module over R of rank 2
  @test iszero(dc[120, 0]) # A zero module
  @test can_compute_index(dc, 100, 0)
  @test !can_compute_index(dc, 100, 2)
  @test !can_compute_index(dc, 100, -3)
  @test iszero(horizontal_map(dc, 120, 0)) # A zero map
  @test can_compute_horizontal_map(dc, 121, 0)
  @test !can_compute_horizontal_map(dc, 121, 5)
  @test !can_compute_horizontal_map(dc, 121, -5)
  @test !iszero(horizontal_map(dc, 0, 0))   # The augmentation map of the resolution
  @test has_upper_bound(dc) && upper_bound(dc) < 1  # Returns `true`
  @test_throws ErrorException vertical_map(dc, 0, 0)
  @test_throws ErrorException horizontal_map(dc, 0, 2)   # An illegitimate request throwing an error as indicated by the above output
end
