Short-term todo list
====================

* *DONE* Implement AMR/physical bounds client-side via type-punning.
    * *DONE* Special handling for send_ghost_zone_async() (aka operator(), aka interpolation)
    * *DONE* Special handling for get_sibling()
* *DONE* Set up AMR/physical bounds in (set|tie).*sibling functions.
* *DONE* Set up physical boundaries for root.
* *DONE* Implement enforce_boundaries(). 
* *DONE* Implement exec_function().
* Implement refinement. 
* *DONE* Implement step():
    * Push semantics.
* *DONE* Implement sub_step():
    * *DONE* Implement compute_(x|y|z)_flux() 
    * Implement adjust_(x|y|z)_flux()
    * *DONE* Implement sum_(x|y|z)_difs()
    * *DONE* Implement add_difs()
    * Add state injection from children.
    * Push semantics.
* *DONE* Add interface for defining the problem (e.g. defining the science_table).
    * *DONE* Add defaults to science_table before passing it to the user.
* *DONE* Add I/O for verification.
* Implement default_main.
* *DONE* A lot of recursive functions need to be called from within step/sub-step to break the global barriers.
* An off switch would be nice.
* Visit control component.

Long-term todo list
===================

* Application INI helpers. A paradigm is needed here to make this not suck.
* More descriptive names (this is a problem both internally and at the interface level).
* --help should show INI options (maybe).
* Investigate using DBShowErrors to silence the error dumps that SILO does (their error messages get mangled by something).

