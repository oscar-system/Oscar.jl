```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Vinberg's algorithm

 A Lorentzian lattice $L$ is an integral $\Z$-lattice of signature $(s_+, s_-)$ with $s_+=1$ and $s_->0$. 
 $L$ is called reflective if the set of fundamental roots $\~{R}(L)$ is finite.

 See for example [Tur17](@cite) for the theory of Arithmetic Reflection Groups and Reflective Lorentzian Lattices.

## Description 
 The algorithm constructs a fundamental polyhedron $P$ for a Lorentzian lattice $L$ by computing its fundamental roots $r$.

 Choose $v_0$ in $L$ with $v_0^2 > 0$ as a point that $P$ should contain.

 Let $Q$ be the corresponding gram matrix of $L$. A fundamental root $r$ satisfies
 - the vector $r$ is primitive
 - reflection by $r$ preserves the lattice, i.e. $\frac{2}{r^2}*r*Q$ is an integer matrix.
 - the pair $(r, v_0)$ is positive oriented, i.e. $(r, v_0) > 0$
 - the product $(r, \~{r}) \geq \ 0$ for all roots $\~{r}$ already found
 This implies that $r^2$ divides $2*i$ for $i$ being the level of $Q$, i.e. the last invariant of the smith normal form of $Q$. 

 $P$ can be constructed by solving $(r, v_0) = n$ and $r^2 = k$ by increasing order of the value $\frac{n^2}{k}$ and $r$ satisfying the above conditions.

 Perhaps $v_0$ lies on one or on several roots. Then $P$ is cannot be uniquely determined.
 In that case we need a direction vector $v_1$ that satisfies
 - the pair $(v_0, v_1)$ lies orthogonal on each other, i.e. $(v_0, v_1) = 0$
 - for all possible roots $\~{v}$ with $(v_0, \~{v}) = 0$ it holds $(\~{v}, v_1) \neq 0$ 

 With $v_0$ and $v_1$ fixed $P$ can be uniquely determined for any choice of root lengths and maximal distance $(v_0, r)$.
 We choose the first roots $r$ by increasing order of the value $\frac{(\~{r}, v_1)}{r^2}$ for all possible roots $\~{v}$ with $(v_0, \~{v}) = 0$.
 For any other root length we continue as stated above.
 
 For proofs of the statements above and further explanations see [Vin75](@cite).

 ## Function
 
 ```@docs
 vinberg_algorithm(Q::ZZMatrix, upper_bound::ZZRingElem)
 ```


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en).
* [Stevell Muller](https://www.math.uni-sb.de/ag/brandhorst/index.php?option=com_content&view=article&id=30:muller&catid=10&lang=de&Itemid=104),

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).