Preprocessing
=============



Default
-------

Pyaflowa
--------

``Pyaflowa`` is the name of the Preprocessing module in SeisFlows that takes 
advantage of the misfit quantification i
`Pyatoa <https://github.com/adjtomo/pyatoa>`__. Pyatoa provides more control 
over the preprocessing step by providing Users with the ability to window
waveforms, choose from a variety of misfit functions/ adjoint source types,
create useful waveform comparison and map figures, and provide an aggregate 
interface to the misfit information that allows the User to perform statistical
analyses on misfit windows and misfit information.

Re-Running Preprocessing
------------------------

Users may wish to stop a workflow after the generation of forward synthetics, 
and run through preprocessing manually, to make decisions on parameters like
misfit function, waveform bandpass, or windowing parameters. 

seisflows debug
preproces.setup()
# pick source
source_name = solver.source_names[0]

# If you change parameters directly in the paramters.yaml file
# Run misfit quant directly on the login node
preprocess.quantify_misfit(source_name=source_name, iteration=99, step_count=99)

- Look at the waveform figures and edit parameters as you see fit
- Either edit parameters in 'preprocess._config.pyflex_config' or change in the 
  parameter file and then rerun steps above to load in new parameters
