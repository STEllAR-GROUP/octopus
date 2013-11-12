Basic Operations
----------------

------------------------------- -------------------------------------------------
Operation                       Behavior 
------------------------------- -------------------------------------------------
set_time                        top-down
clear_refinement_marks          top-down
communicate_ghost_zones         neighbor-to-neighbor *
child_to_parent_state_injection bottom-up
child_to_parent_flux_injection  bottom-up
require_child                   neighbor-to-neighbor, bottom-up **
apply                           top-down
reduce, reduce_zonal            top-down
slice                           top-down

Refinement Operations
---------------------

------------------------------- -------------------------------------------------
Operation                       Behavior 
------------------------------- -------------------------------------------------
mark                            neighbor-to-neighbor, bottom-up, top-down *** ^
populate                        top-down ^ 
link                            neighbor-to-neighbor, top-down **** ^
remark                          neighbor-to-neighbor, bottom-up, top-down ***** ^ 

^ these functions will invoke their kernel before recursing to their children.
FIXME: not sure all of these are labelled.

* communicate_ghost_zones performs communication between a node and all of its
neighbors and with all of the nodes that depend on it for AMR information.

** require_child calls propagate, which may perform a propagation of child
requirements to some of its neighbors. Propagation will occur on all the AMR
boundaries of the parent that the child faces (up to three). propagate and
require_child are sibling-recursive calls that traverse up the tree; their
termination condition is a lack of AMR boundaries (eventually the call chain
will reach a node in which the three relevant faces are either normal neighbors
and/or physical boundaries). TL;DR - complex recursive call structure heading
in an outward and upward direction, dependent on runtime information.

*** mark may invoke require_child on up to three neighboring nodes for each
child that is marked for creation (e.g. theoretical maximum of 8*3 calls). Some
of these calls may actually be redundant.

*** link is applied recursively in the usual top-down manner, invoking
link_kernel on each node. The link_kernel function may communicate directly
with its children, neighbors and indirectly with the children of its neighbors.
The communication will go no further than the aforementioned three groups. Who
link_kernel communicates with is driven by the types of neighbors the node has,
which children the node has marked for creation, which children the node
already has and which children its neighbors have. 
 
**** remark is called on the entire tree in parallel. It is responsible for
calling require_child on each node's "corner" and "edge" relatives (which are
not directly stored by reference in each node). It may call require_child up to
four times (one corner, three edges) for each child marked (e.g. theoretical maximum
of 8*4 calls). Some of these calls may be redundant. 
