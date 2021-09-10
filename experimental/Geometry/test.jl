using Oscar
using Oscar.Misc
  IA2 = affine_space( QQ, 2 )
  R = ambient_ring( IA2 )
  vars = gens(R)
  x = vars[1]
  y = vars[2]
  IA2x = localize( IA2, x*(x+y) )
  IA2xy = localize( IA2x, x*y )
  defining_ideal( IA2xy )

# f = x^2+y^2-1
# X = subscheme( IA2, f )
# Y = subscheme( IA2x, f )
#
U = localize( IA2, x )
V = localize( IA2, y )
(W, i, j) = intersect( IA2, V )
pullback(j)
#
# IA2y = localize( IA2, y )
# IA2xy = localize( IA2x, y ) 
# IA2yx = localize( IA2y, x )
#
# g = Glueing(inclusion_in_parent( IA2xy ),
#             inclusion_in_parent( IA2yx ),
#             AffSchMorphism( IA2xy, IA2yx, gens(R) )
#             )
# 
