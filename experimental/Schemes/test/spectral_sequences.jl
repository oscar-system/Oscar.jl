@testset "Spectral sequences for toric varieties" begin
  X = projective_space(NormalToricVariety, 2) # IP^2 as a toric variety
  S = cox_ring(X)
  G = grading_group(S)

  F = graded_free_module(S, [zero(G)])
  (x, y, z) = gens(S)
  I, _ = sub(F, [x*y*z*F[1]]) # OO(-3)
  J, _ = sub(F, [x^4*y*z*F[1]]) # OO(-6)
  phi = hom(J, I, [x^3*I[1]]) # a morphism OO(-3) <-- OO(-6)

  # We wrap the above morphism in a 2-term complex. 
  cplx = Oscar.SimpleComplexWrapper(Oscar.ComplexOfMorphisms([phi]))
  res = total_complex(Oscar.CartanEilenbergResolution(cplx))

  # And compute its spectral sequence.
  css = Oscar.CohomologySpectralSequence(X, res);
  # At the moment the code for spectral sequences can not handle 
  # dynamic exponent vectors in the toric context. We are working on 
  # fixing this. In the meantime, the user can set the exponent 
  # vector manually to a global default using the following method.
  Oscar.set_global_exponent_vector!(css, 10)
  p1 = css[1]; # The first page. Printing this still throws an error! sorry...

  # The indexation is a bit awkward (due to internal reasons). 
  # The first index is for the direction of the `phi`-complex.
  # The second index ranges from 0 down to -2 and stands for the Cech-cohomology. 
  # This will give you an impression of the first page:
  @test [ngens(p1[i, j]) for j in 0:-1:-2, i in 0:1] == [0 0; 0 0; 1 10]

  # The map between the two non-trivial entries looks like this:
  @test !is_zero(matrix(map(p1, 1, -2)))

  # Let's look at the second page, then:
  p2 = css[2];
  [ngens(p2[i, j]) for j in 0:-1:-2, i in 0:1]

  # Beware that the entry at (0, -2) is zero!
  @test [!is_zero(p2[i, j]) for j in 0:-1:-2, i in 0:1] == [0 0; 0 0; 0 1]
end

