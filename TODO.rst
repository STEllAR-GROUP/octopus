Short-term todo list
====================

* *DONE* Implement AMR/physical bounds client-side via type-punning.
    * *DONE* Special handling for send_ghost_zone_async() (aka operator(), aka interpolation)
    * *DONE* Special handling for get_sibling()
* *DONE* Set up AMR/physical bounds in (set|tie).*sibling functions.
* *DONE* Set up physical boundaries for root.
* *DONE* Implement enforce_boundaries(). 
* *DONE* Implement exec_function().
* Refinement. 
    * *DONE* Implementation.
    * Testing/debugging.
* Regridding.
    * Implementation.
    * Testing/debugging.
* Temporal prediction.
    * Implementation.
    * Testing/debugging.
* *DONE* Implement step():
* Implement sub_step():
    * *DONE* Implement compute_(x|y|z)_flux() 
    * Implement adjust_(x|y|z)_flux()
    * *DONE* Implement sum_(x|y|z)_difs()
    * *DONE* Implement add_difs()
    * *DONE* Add state injection from children.
* *DONE* Add interface for defining the problem (e.g. defining the science_table).
    * *DONE* Add defaults to science_table before passing it to the user.
* *DONE* Add I/O for verification.
* *DONE* A lot of recursive functions need to be called from within step/sub-step to break the global barriers.
* Visit control component.
    * *DONE* Implementation.
    * Testing/debugging.
    * Better cleanup during unexpected shutdown (e.g. signals, crashes).
* User-defined configuration (necessary for dynamic updating of config data).
* Android integration.
 
Long-term todo list
===================

* Application INI helpers. A paradigm is needed here to make this not suck.
* More descriptive names (this is a problem both internally and at the interface level).
* --help should show INI options (maybe).
* Investigate using DBShowErrors to silence the error dumps that SILO does (their error messages get mangled by something).
* Timestep refinement.
* Place grid nodes based on where the corresponding grid node from the previous timestep is located.
* Combine config() and science().
* Rename ghost_zone to ghost_region (to avoid naming clashes with zones, which is the term we use for each independent discrete value).
* Switch to std::array.
* Add grid "age".
* Add refine_or_unrefine to refinement_criteria interface.
* Implement default_main.
* Replace hpx::wait() with continuations and hpx::wait_all().

