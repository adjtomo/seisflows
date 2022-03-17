"""
SeisFlows3 messages tool. For providing a uniform look to SeisFlows3 print
and log statements.
"""
from textwrap import wrap


def mjr(val, char="="):
    """
    Message formatter used to block off sections in log files with visually
    distinctive separators. Defined as individual functions to reduce call
    length.

    Major: For important or workflow.main() messages like starting workflow

    .. rubric::
        >>> print(msg.mjr("Important message here"))
        or
        >>> logger.info.(msg.mjr("Important message here"))

    :type val: str
    :param val: formatted message to return
    :type char: str
    :param char: border character to separate the message from remainder of logs
    :rtype: str
    :return: formatted string message to be printed to std out
    """
    return f"\n{char*80}\n{val:^80s}\n{char*80}"


def mnr(val, char="/"):
    """
    Message formatter used to block off sections in log files with visually
    distinctive separators. Defined as individual functions to reduce call
    length.

    Minor: For key messages, describing things like what iteration were at

    .. rubric::
        >>> print(msg.mnr("Semi important message here"))
        OR
        >>> logger.info.(msg.mnr("Semi important message here"))

    :type val: str
    :param val: formatted message to return
    :type char: str
    :param char: border character to separate the message from remainder of logs
    :rtype: str
    :return: formatted string message to be printed to std out
    """
    return f"\n{char * 80}\n{val:^80s}\n{char * 80}"


def sub(val, char="-"):
    """
    Message formatter used to block off sections in log files with visually
    distinctive separators. Defined as individual functions to reduce call
    length.

    Sub: For sub-critical messages, describing things like notes and warnings

    .. rubric::
        >>> print(msg.mnr("Sub-critical message here"))
        OR
        >>> logger.info.(msg.sub("Sub-critical message here"))


    :type val: str
    :param val: formatted message to return
    :type char: str
    :param char: border character to separate the message from remainder of logs
    :rtype: str
    :return: formatted string message to be printed to std out
    """
    return f"\n{val}\n{char*80}"


def cli(text="", items=None, wraplen=80, header=None, border=None, hchar="/"):
    """
    Provide a standardized look to the SeisFlows command line interface messages
    The look we are after is something like:


    $ seisflows cmd

        =======================
                HEADER
                //////

        text

        item1
        item2
        ...
        itemN
        =======================

    $ ls -l

    .. rubric::
        >>> print(msg.cli("stdout text here", items=["a", "b", "c"],\
                          header="warning", border="="))
        ========================================================================
                                        WARNING
                                        ///////
        stdout text here

        a
        b
        c
        ========================================================================

    :type text: str
    :param text: text to format into the cli look
    :type items: list
    :param items: optional list of items that will be displayed on new lines
        after the text. Useful for listing parameters or paths. The items here
        are NOT wrapped.
    :type wraplen: int
    :param wraplen: desired line length to wrap messages.
    :type header: str
    :param header: optional header line that will be centered (wraplen/2) and
        capitalized. Useful for things like 'WARNING' and 'ERROR'
    :type border: str
    :param border: a character to use to block off
    :type hchar: str
    :param hchar: character to underline the header with
    :rtype output_str: str
    :return output_str: formatted string to print out
    """
    # Start with a newline to space from command line arg
    output_str = "\n"
    # Add top border
    if border is not None:
        output_str += f"{border * wraplen}\n"
    # Add header below top border and a line below that
    if header is not None:
        output_str += f"{header.upper():^{wraplen}}\n"
        output_str += f"{hchar * len(header):^{wraplen}}\n"
    # Format the actual input string with a text wrap
    if text:
        output_str += "\n".join(wrap(text, width=wraplen,
                                     break_long_words=False))
    # Add list items in order of list
    if items:
        # Sometimes text is blank so we don't need the double newline
        if text:
            output_str += "\n\n"
        output_str += "\n".join(items)
    # Add bottom border
    if border is not None:
        output_str += f"\n{border * wraplen}"
    # Final newline to space from next cli
    output_str += "\n"
    return output_str


def write_par_file_header(f, paths_or_parameters, name="", tabsize=4,
                          border="=", uline="/"):
    """
    Re-usable function to write docstring comments inside the SeisFlows3
    parameter file. Used by seisflows.SeisFlows.configure()

    Headers look something like this

    # ===========================
    #       MODULE NAME
    #       ///////////
    # PAR (type):
    #     description of par
    # ===========================
    PAR: val

    :type f: _io.TextIO
    :param f: open text file to write to
    :type paths_or_parameters: dict
    :param paths_or_parameters: the paths or parameters that should be written
        to the header
    :type name: str
    :param name: the name of the module that is being written, will be used as
        the header of the docstring
    :type tabsize: int
    :param tabsize: how large to expand tab character '\t' as spaces
    :type border: str
    :param border: character to use as the header and footer border
    :type uline: str
    :param uline: how to underline the header
    """
    # Some aesthetically pleasing dividers to separate sections
    # Length 77 ensure that total line width is no more than 80 characters
    # including the '#' and spaces
    top = (f"\n# {border * 77}"
           f"\n# {name.upper():^77}"
           f"\n# {uline * len(name):^77}"
           f"\n"
           )
    bot = f"# {border * 77}\n"

    # Write top header, all parameters, types and descriptions, and then footer
    f.write(top)
    for key, attrs in paths_or_parameters.items():
        if "type" in attrs:
            f.write(f"# {key} ({attrs['type']}):\n")
        else:
            f.write(f"# {key}:\n")
        docstrs = wrap(attrs["docstr"], width=77 - tabsize,
                       break_long_words=False)
        for line, docstr in enumerate(docstrs):
            f.write(f"#\t{docstr}\n".expandtabs(tabsize=tabsize))
    f.write(bot)


def write_par_file_paths_pars(f, paths_or_parameters, indent=0, tabsize=4):
    """
    Re-usable function to write paths or parameters in yaml format to the
    SeisFlows3 parameter file. Used by seisflows.SeisFlows.configure()

    Parameters are written something like:

    PAR1: val1
    PAR2: val2
    Par3:
        - val3a
        - val3b
        - val3c

    :type f: _io.TextIO
    :param f: open text file to write to
    :type paths_or_parameters: dict
    :param paths_or_parameters: the paths or parameters that should be written
        to the header
    :type indent: int
    :param indent: level of indentation to match yaml style. passed to
        str.expandtabs(tabsize=`indent`)
    :type tabsize: int
    :param tabsize: how large to expand tab character '\t' as spaces
    """
    for key, attrs in paths_or_parameters.items():
        # Lists need to be treated differently in yaml format
        if isinstance(attrs["default"], list):
            f.write(f"{key}:\n")
            for val in attrs["default"]:
                f.write(f"\t- {val}\n".expandtabs(tabsize=tabsize))
        else:
            # Yaml saves NoneType values as 'null' or blank lines
            if attrs["default"] is None:
                f.write(f"\t{key}:\n".expandtabs(tabsize=indent))
            else:
                f.write(
                    f"\t{key}: {attrs['default']}\n".expandtabs(tabsize=indent)
                )


SystemWarning = """

Please double check SYSTEM parameter

    Expected hostname: {}
    Actual hostname: {}

"""

ParameterWarning_SPECFEM = """

PARAMETER WARNING

    There is a conflict between parameters.

    SPECFEM Parameter:  "{}"
    Old Value:  {}
    Overwriting with:  {}

"""

DataFormatWarning = """

DATA FORMAT WARNING

    reader format: {}
    writer format: {}

    Incompatible file formats may result in job failure or other problems.

"""

TaskIDWarning = """

    WARNING: system.taskid() OS environment variable 'SEISFLOWS_TASKID' not 
    found, SeisFlows3 is assuming debug mode and returning taskid=0. 
    If you are not running in debug mode, please check your SYSTEM.run() or
    SYSTEM.run_single() commands, which are responsible for setting 
    the 'SEISFLOWS_TASKID'

"""

FileError = """

FILE NOT FOUND

    {file}

"""

SolverError = """

SOLVER FAILED

    Nonzero exit status returned by the following command:  
    
    {exc}

    Subsequent tasks may fail because expected solver output is not in place.
    Users running on clusters without fault tolerance should consider stopping 
    any remaining workflow tasks to avoid further loss of resources. 

    To troubleshoot solver errors, navigate to ./scratch/solver to browse solver
    output or try running solver manually in the directories set up in
    ./scratch/solver. 

"""


ReceiverError_SPECFEM = """

ERROR READING RECEIVERS

    Error reading receivers.

"""

SourceError_SPECFEM = """

ERROR READING SOURCES

    In DIRECTORY, there must be one or more files matching WILDCARD.

    DIRECTORY:  "{}"
    WILDCARD:  "{}"

"""


ReaderError = """

READER ERROR

   Seismic data reader not found.

   PAR.READER must correspond to an entry in seisflows.plugins.readers

"""

WriterError = """

WRITER ERROR

   Seismic data writer not found.

   PAR.WRITER must correspond to an entry in seisflows.plugins.writers

"""

TaskTimeout = """

TASK TIMED OUT

    Stopping workflow because task time limit exceeded. (To adjust limit,
    add or modify TASKTIME in parameter file.)

        Task name:  {classname}.{method}
        Task id:    {job_id}
        Time limit (minutes): {tasktime}

"""

TaskError_LSF = """

TASK ERROR

    Task failed:  {}.{}

    For more information, see output.lsf/{}

    Stopping workflow...

"""


TaskError_PBS = """

TASK ERROR

    Task failed:  {}.{}

    For more information, see output.pbs/{}

    Stopping workflow...

"""

TaskError_SLURM = """

TASK ERROR

    Task failed:  {classname}.{method}

    For more information, see output.logs/{job_id}

    Stopping workflow...

"""

obspyImportError = """

DEPENDENCY ERROR

    The current data processing workflow requires OBSPY.  Please install it and
    try again.

"""

mpiError1 = """

SYSTEM CONFIGURATION ERROR

    The following system configuration can be used only with single-core
    solvers:

        system.{}

    If your solver requires only a single core, then set NPROC equal to 1.

    If your solver requires multiple cores, then consider using lsf_lg, pbs_lg,
    or slurm_lg system configurations instead.

"""

mpiError2 = """

DEPENDENCY ERROR

    The following system configuration requires MPI4PY:

        system.{}

    Please install MPI4PY and try again, or consider choosing a different system
    configuration.

"""

mpiError3 = """

SYSTEM CONFIGURATION WARNING

    The following system configuration requires 'mpiexec':

        system.{}

    Please make sure than 'mpiexec' is accessible through your shell's PATH
    environment variable. If your executable goes by a different name such as
    'mpirun', consider creating an alias in your shell's configuration file, and
    remember to source the modified configuration file. If MPI is not available
    on your system, consider using the 'multithreaded' system interface instead.

"""

MissingParameter_Workflow = """

Please specify a workflow by adding a line to the parameter file, e.g.

    WORKFLOW='inversion';

for a list of available workflows, see seisflows/workflow in the source code

"""

MissingParameter_System = """

Please specify a system interface by adding a line to the parameter file, e.g.

    SYSTEM='serial';

for a list of available interfaces, see seisflows/system in the source code

"""

CompatibilityError1 = """

Parameter settings have changed.

In your parameter file, please remove
    OPTIMIZE='base'

and add one of the following instead
    OPTIMIZE='LBFGS'
    OPTIMIZE'=NLCG'
    OPTIMIZE='SteepestDescent'

"""

Warning_pbs_sm = """

WARNING:  PBS_SM hasn't been tested for a long while because we don't own a PBS
cluster.  If you have access to one cluster and are willing to help debug, 
please let us know.

"""

Warning_pbs_lg = """

WARNING:  PBS_LG hasn't been tested for a long while because we don't own a PBS
cluster.  If you have access to one cluster and are willing to help debug, 
please let us know.

"""

PoissonsRatioError = """
                
ERROR CHECKING POISSON'S RATIO

    The Poisson's ratio of the given model is out of bounds with respect 
    to the defined range ({min_val}, {max_val}). The model bounds were found 
    to be:

    {pmin:.2f} < PR < {pmax:.2f}

"""

ExportResidualsError = """

PREPROCESSING ERROR

    The Solver function 'export_residuals' expected 'residuals' directories to 
    be created but could not find them and cannot continue the workflow. 
    
    Please check the preprocess.prepare_eval_grad() function which is
    responsible for exporting the 'residuals' directory. Or use 
    'seisflows debug' to investigate the error more closely.


"""


DataFilenamesError = """

    ERROR: The property solver.data_filenames, which is used to search for 
    trace data in ./scratch/solver/*/traces is empty and should not be. 
    Please check solver.data_filenames and solver.data_wildcard against filenames
    in the traces/ directory.
    

"""


def check(cls):
    """
    Standardized log message sent from the check() function that is required
    within each module and submodule. Simply notifies the user that checks
    are being performed within the given module

    :type cls: 'type'
    :param cls: self.__class__ or type(self) which should be passed in from
        INSIDE a function defined INSIDE a class.
        e.g., <class 'seisflows3.preprocess.default.Default'>
    :rtype: str
    :return: formatted check string describing the module and classname
        e.g., preprocess.default
    """
    module, clsname = split_cls(cls)
    return f"check paths/pars module: {module}.{clsname}"


def setup(cls):
    """
    Standardized log message sent from the setup() function that is required
    within each module and submodule. Simply notifies the user that checks
    are being performed within the given module

    :type cls: 'type'
    :param cls: self.__class__ which should be passed in from INSIDE a function
        defined INSIDE a class.
        e.g., <class 'seisflows3.preprocess.default.Default'>
    :rtype: str
    :return: formatted string describing the module and classname
        e.g., preprocess.default

    """
    module, clsname = split_cls(cls)
    return f"setting up module: {module}.{clsname}"


def whoami(cls, append="", prepend=""):
    """
    Standardized log message sent from any function inside a class. Used in logs
    for subclasses to declare who they are before running functions, makes it
    easier to track down what's happening where

    :type cls: 'type'
    :param cls: self.__class__ which should be passed in from INSIDE a function
        defined INSIDE a class.
        e.g., <class 'seisflows3.preprocess.default.Default'>
    :type append: str
    :param append: add any string after the whoami statement
    :type prepend: str
    :param prepend: add any string in front of the whoami statment
    :rtype: str
    :return: formatted string describing the module and classname
        e.g., preprocess.default
    """
    module, clsname = split_cls(cls)
    return f"{prepend}{module}.{clsname}{append}"


def split_cls(cls):
    """
    Repeatedly used function to split the output of type(self) or self.__class__
    into separate strings


    :type cls: 'type'
    :param cls: self.__class__ which should be passed in from INSIDE a function
        defined INSIDE a class.
        e.g., <class 'seisflows3.preprocess.default.Default'>
    :rtype: tuple of strings
    :return: module, classname
    """
    type_, str_, brkt = str(cls).split("'")
    sf3, module, fid, clsname = str_.split(".")
    return module, clsname