## For next version 0.1.1
#### General
- [ ] write unit tests
- [ ] update docs to reflect the major changes made, include changelog
- [ ] try to reconcile a `base` parameters.yaml file with `custom` attributes
- [ ] create a simple pet 2D example (specfem2D acoustic checkerboard from Ryan)
      that can be run serial and can be used for testing of the base classes
- [ ] Bring in all classes from Seisflows and update to Py3 w/ docstrings

#### Preprocess
- [ ] Finish updating to Py3 and writing full docstrings, better integration 
      with ObsPy
- [ ] Include the machinery of Pyatoa as a Preprocess class rather than 
      completely separate as Pyaflowa

#### Plugins
##### Note: I always found it really confusing that all of the optmization and line search machinery is set in plugins, and that the optimization base class just calls a plugin to inherit its machinery. What was the purpose of this and can it be cleaned up a little bit?

- [ ] Clean up the random plugins and maybe organize them better

#### Postprocess
- [ ] Finish writing combine_vol_data_vtk wrapper to generate vtk models

#### Solver
- [ ] Re-introduce Specfem2D and Specfem3D Globe solver classes

#### Workflow
- [ ] Add resume_from capability into Inversion base class
- [ ] Finish updating migration and migration pyatoa, clean up forward

#### Scripts
- [ ] add more descriptive help statements, maybe a step by step way to set up a
      run folder that asks the User to choose the parameters they want and fills
      in the parameter.yaml file and maybe also points them to the directories
      that require custom classes
- [ ] include print statement for submit detailing important parameters such
      as ntasks, walltime, begin, end
- [ ] include print statement for resume detailing parameters such as ntasks,
      walltime, begin, end, resume_from
- [ ] create a status function to check on the status of jobs, maybe by 
      querying the output files. Try not to make it system dependent so that it
      can still be generally applicable to Seisflows
- [ ] Move convert_model() function into the util section of Seisflows


