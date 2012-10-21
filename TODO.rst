Short-term todo list
====================

These are the things necessary to refine the initial timestep (e.g. produce the
equivalent of X.initial.silo.gz). This basically gives us everything else,
though.

* *DONE* Implement AMR/physical bounds client-side via type-punning.
    * *DONE* Special handling for send_ghost_zone_async() (aka operator(), aka interpolation)
    * *DONE* Special handling for get_sibling()
* *DONE* Set up AMR/physical bounds in (set|tie).*sibling functions.
* *DONE* Set up physical boundaries for root.
* *DONE* Implement enforce_boundaries(). 
* *DONE* Implement exec_function().
* Implement refinment. 
* *DONE* Implement step():
    * Push semantics
* *DONE* Implement sub_step():
    * Implement compute_(x|y|z)_flux() 
    * Implement adjust_(x|y|z)_flux()
    * Implement sum_(x|y|z)_difs()
    * Implement add_difs()
    * Add state injection from children.
    * Push semantics
* *DONE* Add interface for defining the problem (e.g. defining the science_table).
    * *DONE* Add defaults to science_table before passing it to the user.
* Add I/O for verification.
* Implement default_main.

