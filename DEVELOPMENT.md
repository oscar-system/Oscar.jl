# Oscar - how do we envision development

## What Do We Have:

We have a large codebase for infrastructure in place, comprising at least
 - matrices
 - polynomials (univariate and multivariate)
 - power series
 - number fields
 - (abelian) groups
 - polytopes, cones, linear programs
 - polyhedral fans
 - ... and MUCH more
As well as specialised stuff
 - polymake
 - Singular
 - Gap

So: please use it. It is safe to assume all can be improved, however, if we try to perfect every
single line of code again and again, we we won't get anywhere; there is a balance to be found.
For preference: interoperability is much more
important than speed, readability is much more important than speed. Having
said that: of course, sometimes pure speed matters, but not nearly as often as
people think.

The infrastructure is incomplete, e.g. we do not have
 - varieties, (we've started with curves recently)
 - simplicial complexes, combinatorial manifolds
 - surfaces
 - tropical varieties, in particular hypersurfaces
 - tropical polytopes
 - ....

If you need them, e.g. you want a library/ package for (elliptic) curves, then
step ONE should be to define a curve. It does not need to be perfect, support
all of Hartshorn immediately, but it also should not contradict him either. It
is OK for a function to say: sorry, not (yet) implemented. Then if I need it,
I can add and re-use the infrastructure.

The key idea is to pick a single textbook and follow the conventions in there. This
way we get a consistent interface.

## What Is Missing?

At this point we have, mathematically speaking, almost nothing. This has to
change, quickly.  We work in a team of ~ 10(?) people.
We lack infrastructure as everyone will notice immediately.
As a result, the classical
   I work for 6 month in a separate repo on a branch and then will dazzle you with perfect
   code and cool examples
is not going to work for now. It will result in everyone fixing the same
infrastructure problems over and over again. Please consider to work, e.g. in
a file/directory  in `Oscar/examples` and push on a regular basis, even, or in
particular, incomplete code. Break it down into small pull requests.
 - it is in examples, so no-one is accidentally going to automatically use it
 - it is in the main repo, so everyone can see it and try to use it, producing feedback
 - infrastructure, in Oscar at least, can be added in the same repo, immediately and thus
   is available to others, immediately. (feedback, ...)
 - Talk to us - as early as possible. There is a chance that what you're trying to
   do is already there - in a different package maybe, not linked in, ....
 - Talk to us!
   - on slack: send an email to webmaster-oscar@mathematik.uni-kl.de to be added
   - on oscar-dev@mathematik.uni-kl.de
   - whenever you meet us (at the uni, during conferences, ...)

# Practical Development

For practical development:
 - use [the `Revise` package](https://github.com/timholy/Revise.jl)
 - start a file/ directory in Oscar/examples, call it Blah.jl

Then `Oscar.example("Blah.jl")` will include it, so you don't have to worry about a path;
and `Oscar.revise("Blah.jl")` will put it under revise, i.e. re-load it automatically whenever you make a change.

If Blah.jl starts with

```julia
 module BlahModule
 using Oscar

 ...

 end

 using .BlahModule
``` 

then you can even redefine struct's without quitting Julia (and restarting)


Julia's convention is to use `snake_case`, not `CamelCase` for function names,
so that is the accepted style in Oscar.
Types (and their constructor) tend to be in `CamelCase`. However, for me, as a user,
please ALSO provide the constructor/ a constructor in underscore_case. As a
user I do usually not know if something is a constructor or a function (nor do I
want to).

In Oscar all objects have a default mathematical meaning and we have multiple dispatch. As
a result
 - IsSurjectiveMapForAbelianGroups
 - IsSurjectiveMapForRing
is not necessary
 - is_surjective
will do, the different calling signatures will add the type info. If you have a ring map
and want to test it as a map between additive/ multiplicative groups, then you'll have to 
first convert/coerce/wrap it (in)to a map between groups... (we can always discuss exceptions)
We will have to add "forgetful functor" wrappers.

Whenever you feel the urge to define a new basic type, say a new multivariate
ring, please consider the ramifications:
 - can you do matrices?
 - modules?
 - ideals?
 - graded stuff?
 - "complete" arithmetic?
 - interaction with other types? (map to residue rings, apply automorphisms, ...)
 - in fact, everything the other MPoly type can?

if yes: still reconsider, but go ahead

if no: at least use "our" types to interface your function, better still, use our type and
complain about lack of functionality/ speed/ interface (or provide patches)

This applies to all foundations! They are all incomplete, and they can all
be improved BUT if everyone does their own foundations, we cannot work
together.

The main point about Oscar is for non-experts to allow access, thus it needs to
work as you would expect as the youngest student knowing the objects. Something
like: "but all experts know and expect this - it is always done this way" will
make your function impossible to be used by outsiders. The only vague chance of
having a consistent system is to stick to "global definitions", ie. not to the
"experts" version (for efficiency). Feel free to add the other layer if you
need.

Lastly, so far, unless really really necessary, don't include outside packages
as dependencies. Every outside package means that the testing, releases, ...
are more complicated.  You're at the whim of others.  If you need a highly
non-trivial 10000 line package, well, there is no help and we need to see how
to arrange it. In particular, for the near future, it will be your job to
ensure that this outside package still works with us.  If you need a trivial 10
line function from a 10000 line package, just don't.
