# Run a few SeisFlows commands to generate test data. Useful if the parameter
# file ever changes, at which point we will need to generate new files.
cd ..
rm *yaml

seisflows setup
cp parameters.yaml test_setup_parameters.yaml

seisflows configure -r
cp parameters.yaml test_conf_parameters.yaml

seisflows par -p materials elastic
seisflows par -p density constant
seisflows par -p nt 1000
seisflows par -p dt .01
seisflows par -p f0 .084
seisflows par -p format ascii
seisflows par -p begin 1
seisflows par -p end 1
seisflows par -p case data
seisflows par -p attenuation False
seisflows par -p specfem_bin ./bin
seisflows par -p specfem_data ./DATA
seisflows par -p model_init ./MODEL_INIT
seisflows par -p model_true ./MODEL_TRUE

mv parameters.yaml test_filled_parameters.yaml


