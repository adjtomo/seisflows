#!/usr/bin/env python3
"""
This class provides utilities for the Seisflows solver interactions with
Specfem3D Globe. It is built on top of the Specfem3D solver class.

SPECFEM3D_Globe specific notes:
    - Does not allow SU seismogram outputs, only ASCII, SAC, ASDF, 3D_Array
    - Although SPECFEM3D_GLOBE ASCII synthetics have extension '.ascii'
      the adjoint sources are NOT supposed to have this, i.e., they 
      should have extension '.adj'
        
"""
from seisflows.solver.specfem3d import Specfem3D


class Specfem3DGlobe(Specfem3D):
    """
    Solver SPECFEM3D_GLOBE
    ----------------------
    SPECFEM3D_Globe-specific alterations to the solver.specfem3d module

    Parameters
    ----------
    :type smooth_type: str
    :param smooth_type: choose how smoothing is performed for gradients.
	these are tied to the internal smoothing functions available.
	- 'gaussian': convolve with a 3D gaussian, slow and computationally
	    intensive, but default and matches 2D and 3D_Cartesian smoothing
	- 'laplacian' (default): average points around vertex to smooth. 
	    faster and preferred method for GLOBE code

    Paths
    -----
    ***
    """
    __doc__ = Specfem3D.__doc__ + __doc__

    def __init__(self, smooth_type="laplacian", **kwargs):
        """Instantiate a Specfem3D_Globe solver interface"""
        super().__init__(**kwargs)
	self.smooth_type = smooth_type

	self._acceptable_smooth_types = ["laplacian", "gaussian"]

        if self.materials.upper() == "ACOUSTIC":
            self._parameters += ["vp"]
        elif self.materials.upper() == "ELASTIC":
            self._parameters += ["vp", "vs"]
        elif self.materials.upper() == "ISOTROPIC":
            self._parameters += ["vp", "vs"]
        elif self.materials.upper() == "ANISOTROPIC":
            self._parameters += ["vpv", "vph", "vsv", "vsh", "eta"]

    def check(self):
	"""
	Checks parameter validity for SPECFEM3D_GLOBE parameters
	"""
	super().check()
	
	assert(self.smooth_type in self._acceptable_smooth_types), \
	    f"`smooth_type` must be in {self._acceptable_smooth_types}" 

    def data_wildcard(self, comp="?"):
        """
        Returns a wildcard identifier for synthetic data
        Currently only support for ASCII seismograms

        :rtype: str
        :return: wildcard identifier for channels
        """
        if self.data_format.upper() == "ASCII":
            return f"*.?X{comp}.sem.ascii"

    def forward_simulation(self, executables=None, save_traces=False,
                           export_traces=False, **kwargs):
        """
        Calls SPECFEM3D_GLOBE forward solver, exports solver outputs to traces.
        
        .. note::
            SPECFEM3D_GLOBE does not require 'xgenerate_databases' which is 
            required for Cartesian. This function overwrites that

        :type executables: list or None                                          
        :param executables: list of SPECFEM executables to run, in order, to     
            complete a forward simulation. This can be left None in most cases,  
            which will select default values based on the specific solver        
            being called (2D/3D/3D_GLOBE). It is made an optional parameter      
            to keep the function more general for inheritance purposes.          
        :type save_traces: str                                                   
        :param save_traces: move files from their native SPECFEM output location 
            to another directory. This is used to move output waveforms to       
            'traces/obs' or 'traces/syn' so that SeisFlows knows where to look   
            for them, and so that SPECFEM doesn't overwrite existing files       
            during subsequent forward simulations                                
        :type export_traces: str                                                 
        :param export_traces: export traces from the scratch directory to a more 
            permanent storage location. i.e., copy files from their original     
            location 
        """
        if exeuctables is None:
            executables = ["bin/xspecfem3D"]

	    # Database files only need to be made once, usually at the first
	    # evaluation. Once made, we don't have to run xmeshfem3D anymore.
	    if not glob(os.path.join(self.model_databases, "proc*_Database")):
		executables = ["bin/xmeshfem3D"] + executables

        super().forward_simulation(executables=executables, 
                                   save_traces=save_traces, 
                                   export_traces=export_traces, 
				   **kwargs)

    def adjoint_simulation(self, executables=None, save_kernels=False,
			   export_kernels=False):
	"""
	Supers SPECFEM3D for adjoint solver and removes GLOBE-specific fwd files

	:type executables: list or None                                          
        :param executables: list of SPECFEM executables to run, in order, to     
            complete an adjoint simulation. This can be left None in most cases, 
            which will select default values based on the specific solver        
            being called (2D/3D/3D_GLOBE). It is made an optional parameter      
            to keep the function more general for inheritance purposes.          
        :type save_kernels: str                                                  
        :param save_kernels: move the kernels from their native SPECFEM output   
            location to another path. This is used to move kernels to another    
            SeisFlows scratch directory so that they are discoverable by         
            other modules. The typical location they are moved to is             
            path_eval_grad                                                       
        :type export_kernels: str                                                
        :param export_kernels: export/copy/save kernels from the scratch         
            directory to a more permanent storage location. i.e., copy files     
            from their original location. Note that kernel file sizes are LARGE, 
            so exporting kernels can lead to massive storage requirements.
	"""
	super().adjoint_simulation(executables=executables,                      
			       	   save_kernels=save_kernels,                    
			       	   export_kernels=export_kernels)                
        
	# Working around fact that save_forward arrays have diff naming
        if self.prune_scratch:                                                   
            for glob_key in ["proc??????_reg?_absorb_buffer.bin"]: 
                logger.debug(f"removing '{glob_key}' files from database "       
                             f"directory")                                       
                unix.rm(glob(os.path.join(self.model_databases, glob_key)))

    def smooth(self, input_path, output_path, parameters=None, span_h=None,      
               span_v=None, use_gpu=False):                                      
        """
        Logic function to choose between available smoothing types for GLOBE
                                                                                 
        :type input_path: str                                                    
        :param input_path: path to data                                          
        :type output_path: str                                                   
        :param output_path: path to export the outputs of xcombine_sem           
        :type parameters: list                                                   
        :param parameters: optional list of parameters,                          
            defaults to `self._parameters`                                       
        :type span_h: float                                                      
        :param span_h: horizontal smoothing length in meters                     
        :type span_v: float                                                      
        :param span_v: vertical smoothing length in meters                       
        :type use_gpu: bool                                                      
        :param use_gpu: whether to use GPU acceleration for smoothing. Requires  
            GPU compiled binaries and GPU compute node.                          
        """
        if self.smooth_type == "gaussian":
            super().smooth(input_path=input_path, output_path=output_path, 
                           parameters=parameters, span_h=span_h, span_v=span_v,
                           use_gpu=use_gpu)
        elif self.smooth_type == "laplacian":
            self.smooth_laplacian(
                    input_path=input_path, output_path=output_path, 
                    parameters=parameters, span_h=span_h * 1E-3, 
                    span_v=span_v * 1E-3
                    )
                                  

    def smooth_laplacian(self, input_path, output_path, parameters=None, 
			 span_h=None, span_v=None):
        """                                                                      
        Wrapper for SPECFEM binary: xsmooth_laplacian_sem

        Smooths kernels by with Laplacian smoothing which takes averages of a
	mesh corner with all it's surrounding points.

        .. note::
            Externally this smooth function behaves almost identically to
            the normal gaussian smoothing function
                                                                                 
        .. note::                                                                
            It is ASSUMED that this function is being called by                  
            system.run(single=True) so that we can use the main solver           
            directory to perform the kernel smooth task                          
                                                                                 
        :type input_path: str                                                    
        :param input_path: path to data                                          
        :type output_path: str                                                   
        :param output_path: path to export the outputs of xcombine_sem           
        :type parameters: list                                                   
        :param parameters: optional list of parameters,                          
            defaults to `self._parameters`                                       
        :type span_h: float                                                      
        :param span_h: horizontal smoothing length in km
        :type span_v: float                                                      
        :param span_v: vertical smoothing length in km
        """                                                                      
        unix.cd(self.cwd)                                                        
                                                                                 
        # Assign some default parameters from class attributes if not given      
        if parameters is None:                                                   
            parameters = self._parameters                                        
        if span_h is None:                                                       
            span_h = self.smooth_h 
        if span_v is None:                                                       
            span_v = self.smooth_v
                                                                                 
        logger.debug(f"smoothing {parameters} with laplacian, horizontal span "
                     f"{span_h}m and vertical span {span_v}m")               


        # NOTE: Converting smoothing lengths 'm' -> 'km' as laplacian smoothing
        #   function is epxecting things in 'km' while SeisFlows expects things
        #   in 'm'
        span_h *= 1E-3
        span_v *= 1E-3

        if not os.path.exists(output_path):                                      
            unix.mkdir(output_path)                                              
                                                                                 
        # Ensure trailing '/' character, required by xsmooth_sem                 
        input_path = os.path.join(input_path, "")                                
        output_path = os.path.join(output_path, "")                              

        # mpiexec ./bin/xsmooth_laplacian_sem SIGMA_H SIGMA_V name input output
        for name in parameters:                                                  
            exc = (f"bin/xsmooth_laplacian_sem {str(span_h)} {str(span_v)} "
                   f"{name}_kernel {input_path} {output_path}")
            # e.g., combine_vs.log                                               
            stdout = f"{self._exc2log(exc)}_{name}.log"                          
            self._run_binary(executable=exc, stdout=stdout)                      
                                                                                 
        # Rename output files to remove the '_smooth' suffix which SeisFlows     
        # will not recognize                                                     
        files = glob(os.path.join(output_path, "*"))                             
        unix.rename(old="_smooth", new="", names=files)  
