Possible Physics Issues
=======================
* We're summing flowoffs incorrectly, they are always 0 in the root grid if we 
  use more than LOR.

Possible Performance Problems
=============================
* GIDs of engine servers should be stripped everywhere.
* GIDs of octree components should never be passed with credits to the children;
  the lifetime of these components is deterministic and well understood.
* Some of our operations that iterate over our children assume that we have all
  eight children for certain things (like allocating temporary storage), because
  we don't store the number of actual children that we have anywhere.
 
Runtime/INI Configuration
=========================
* Configuration data should be stored in one structure instead of in individual
  global_variable objects.
* global_variable is probably unnecessary.
* Octopus should provide a base configuration class.
** Applications can Clang-style casting dynamic casting of base class pointers
   to the application's derived class.
*** In release these casts should not do any checking to ensure we maintain
    fast access to the parameters
** Octopus should provide two "global" instances of the configuration class:
*** One will be the "running" configuration, which the application is using.
*** One will be the "buffered" configuration.
** Octopus should provide a function that copies the "buffered" configuration
   to the running configuration (similar to Cisco routers).

Science Table/State Class
=========================
* The science table should be deprecated, and replaced with a state base class,
  similar to Dominic's original design.

Pluggable Solver Infrastructure
===============================
* The KT solver should be extracted from the core octree code and encapsulated
  in an implementation that fulfills some sort of Solver concept.

I/O
===
* I/O layer needs to be completely rewritten; it is not versatile enough and it
  is brittle.
** Should be possible to specify multiple output backends.
* Reimagined role of I/O subsystem; should be more general purpose, capable of
  handling auxiliary I/O (e.g. we want to use it for more than just outputting
  the system state).
** Checkpointing subsystem should share common infrastructure with the general
   I/O subsystem.
** Change to Octopus architecture: checkpointing will now sit on top of the
   the I/O subsystem.

Generic Octree Algorithms
=========================
* *Replace generic actions that take hpx::util::function with templated actions*
** octree_server::apply, octree_server::apply_leaf
** octree_server::reduce, octree_server::reduce_zonal

Dead Code/Unnecessary Code
==========================
* Visit code should either be rewritten, removed or moved to an archive branch.
* Performance of indexer2d's index calculations should be investigated.

