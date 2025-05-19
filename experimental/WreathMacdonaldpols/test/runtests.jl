@testset "wreath macdonald polynomials" begin
  K,_ = abelian_closure(QQ)
  parent = polynomial_ring(K, [:q,:t], cached=false)
  Q = fraction_field(parent[1])
  (q,t)=parent[2]

  result_p1 = matrix(Q,[t^2 t 1; q t 1; q q^2 1])
  @test result_p1 == wreath_macdonald_polynomials(1,3,(@perm 3 (1,2,3)),[0,1,-1],parent)

  result_p2 =matrix(Q,[[1//q^4, 1//q, (q^3 + 1)//q^3, 1//q^2, q, (q^3 + 1)//q^2, (q^3 + 1)//q, 1, q^3],
              [t//q^3, t^3//q^4, (q*t + t^3)//q^3, t//q,  t^3//q^2, (q^2*t + t^2)//q^3, (q^2*t + t^2)//q^2, 1, t^2//q],
              [1//(q^2*t), 1//(q^2*t), (q + t^2)//(q^2*t), t//q, t//q, (q^2 + t)//(q^2*t), (q*t + 1)//q, 1, 1],
              [1//t^2, q//t^4, (q + t^2)//t^2, t^2, q, (q + t^2)//t^3, (q + t^2)//t, 1, q//t^2],
              [1//t^2, t, t^3 + 1, t^2, t^5, (t^3 + 1)//t, t^4 + t, 1, t^3],
              [1//q, t, q*t + 1, q, q^2*t, (q^2 + t)//q, q^2 + t, 1, q*t],
              [1//t^2, 1//t^2, (q + t^2)//t^2, q, q, 2//t, (q + t^2)//t, 1, 1],
              [q//t, q^4//t, (q^5 + q^2)//t, q^3//t, q^6//t, (q^3 + q*t)//t, (q^4 + q^2*t)//t, 1, q^3],
              [1//t^2, 1//q, (q + t^2)//t^2, q^2//t^2, q, (q + t^2)//(q*t), (q + t^2)//t, 1, t^2//q]])
  @test result_p2 == wreath_macdonald_polynomials(2,3,(@perm 3 ()),[1,-1,0],parent)

  result_p3 = matrix(Q,[q^2   t   q^3 + q*t   q^4   q^2*t   (q^2 + t)//q   q^2 + t   1   t//q^2])
  @test result_p3 == wreath_macdonald_polynomial(multipartition([[1,1],[],[]]),(@perm 3 (1,2,3)), [0,1,-1],parent)
end
