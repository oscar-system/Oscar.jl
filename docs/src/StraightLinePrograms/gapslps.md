# GAP's SLPs

There are two other available SLP types: `GAPSLProgram` and `AtlasSLProgram`,
and related `GAPSLDecision` and `AtlasSLDecision`, which are constructed
similarly as in GAP:

```julia
julia> prg = GAPSLProgram( [ [1,2,2,3], [3,-1] ], 2 )
# input:
r = [ g1, g2 ]
# program:
r[3] = r[1]^2*r[2]^3
r[4] = r[3]^-1
# return value:
r[4]

julia> SLP.evaluate(prg, [perm1, perm2])
(1,3,4,2)

julia> SLP.evaluate(prg, [x, y])
#1 = ^  x  2  ==>  x^2
#2 = ^  y  3  ==>  y^3
#3 = * #1 #2  ==>  (x^2y^3)
#4 = ^ #3 -1  ==>  (x^2y^3)^-1
return: #4

julia> SLProgram(prg) # direct compilation (with room for optimizations obviously)
#1 =    x     ==>  x
#2 =    y     ==>  y
#3 = ^ #1  2  ==>  x^2
#4 = ^ #2  3  ==>  y^3
#5 = * #3 #4  ==>  (x^2y^3)
#3 =   #5     ==>  (x^2y^3)
keep: #1..#3
#4 = ^ #3 -1  ==>  (x^2y^3)^-1
keep: #1..#4
return: #4

julia> GAPSLProgram( [ [2,3], [ [3,1,1,4], [1,2,3,1] ] ], 2 )
# input:
r = [ g1, g2 ]
# program:
r[3] = r[2]^3
# return values:
[ r[3]*r[1]^4, r[1]^2*r[3] ]

julia> GAPSLDecision([ [ [ 1, 1, 2, 1 ], 3 ], [ "Order", 1, 2 ], [ "Order", 2, 3 ], [ "Order", 3, 5 ] ] )
# input:
r = [ g1, g2 ]
# program:
r[3] = r[1]*r[2]
order( r[1] ) == 2 || return false
order( r[2] ) == 3 || return false
order( r[3] ) == 5 || return false
# return value:
true

julia> SLProgram(ans)
#1 =    x     ==>  x
#2 =    y     ==>  y
#3 = * #1 #2  ==>  (xy)
keep: #1..#3
test: order(#1) == 2 || return false
test: order(#2) == 3 || return false
test: order(#3) == 5 || return false
return: true

julia> d = AtlasSLDecision("inp 2\nchor 1 2\nchor 2 3\nmu 1 2 3\nchor 3 5")
inp 2
chor 1 2
chor 2 3
mu 1 2 3
chor 3 5

376> SLP.evaluate(d, [perm1, perm2])
false

julia> GAPSLDecision(d)
# input:
r = [ g1, g2 ]
# program:
order( r[1] ) == 2 || return false
order( r[2] ) == 3 || return false
r[3] = r[1]*r[2]
order( r[3] ) == 5 || return false
# return value:
true

julia> SLProgram(d)
#1 =    x     ==>  x
#2 =    y     ==>  y
test: order(#1) == 2 || return false
test: order(#2) == 3 || return false
#3 = * #1 #2  ==>  (xy)
keep: #1..#3
test: order(#3) == 5 || return false
return: true
```

