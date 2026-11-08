function check_char(A::Ring, B::Ring)
  @req characteristic(A) == characteristic(B) "wrong characteristic"
end
