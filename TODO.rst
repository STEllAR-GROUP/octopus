Short-term todo list
====================

These are the things necessary to refine the initial timestep (e.g. produce the
equivalent of X.initial.silo.gz). This basically gives us everything else,
though.

* Add state injection from children.
* Set up physical boundaries for root.
* Implement exec_function().
* Implement max_dt().
* Implement enforce_boundaries(). 
* Implement AMR/physical bounds client-side via type-punning.
* Set up AMR/physical bounds in (set|tie).*sibling functions.
* Add interface for defining the problem (e.g. defining the science_table).
* Implement step():
    * Implement check_for_refine()
* Implement sub_step():
    * Implement compute_(x|y|z)_flux() 
    * Implement adjust_(x|y|z)_flux()
    * Implement sum_(x|y|z)_difs()
    * Implement add_difs()
* Add I/O for verification.
 
