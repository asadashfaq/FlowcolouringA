Flow Tracing
============

As of July 10 2014 this repository is now depending on *aurespf* as a submodule
rather than old snapshots of the solver.

The *_less.py* files are clones of the old files with similar names. The
difference is updated (less) dependencies.

To run the flow tracing algorithm refer to the wrapper function *calc_usage*
in *usage_old.py*.

For vector flow tracing look at *vector.py*. It contains an implementation
that builds on top of the up/down stream flow tracing algorithm.
