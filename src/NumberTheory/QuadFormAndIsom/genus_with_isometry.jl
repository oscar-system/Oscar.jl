function is_in_same_genus(L1::ZZLatWithIsom, L2::ZZLatWithIsom)
  P,x = polynomial_ring(QQ, x; cached=false)  
  charpoly(P,L1) == charpoly(P,L2) || return false
  m = minpoly(P,L1)
  m == minpoly(P,L2) || return false
  for (d,e) in factor(m)
    if e>1
      error("not implemented for non-semisimple isometries")
    end 
    H1 = hermitian_structure(kernel_lattice(L1, d))
    H2 = hermitian_structure(kernel_lattice(L2, d))
    genus(H1) == genus(H2) || return false 
    
  end 
end 
