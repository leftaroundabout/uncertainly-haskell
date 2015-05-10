A library for dealing with values that are known not exactly but only to some roughly known precision, because obtaining the exact value is impossible or at least infeasible and perhaps not necessary.
It is nontrivial to see how the errors in such quantities propagate through further computations; this library provides a functor interface that allows "embedding" of any such computations into the uncertainty-range.

**This project is more or less abandoned. I continue work on &ldquo;uncertain values&rdquo; in Haskell, but within the [manifolds](https://github.com/leftaroundabout/manifolds) project, which defines the setting needed to define what &ldquo;deviations&rdquo; are even supposed to be, in anything more general than scalar types / single physical measurements.**
