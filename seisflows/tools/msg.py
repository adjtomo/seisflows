"""
SeisFlows messages tool. For providing a uniform look to SeisFlows print
and log statements that end up in stdout or in log files.
"""
from textwrap import wrap

# Unicode degree symbol for log statements etc.
DEG = u"\N{DEGREE SIGN}"


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
    output_str = ""
    # Add top border
    if border is not None:
        output_str += f"\n{border * wraplen}\n"
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
    # output_str += "\n"
    return output_str


# Template parameter file used for 'seisflows setup' to initiate a workflow
base_parameter_file = """
# //////////////////////////////////////////////////////////////////////////////
#
#                        SeisFlows YAML Parameter File
#
# //////////////////////////////////////////////////////////////////////////////
#
# Modules correspond to the structure of the source code, and determine
# SeisFlows' behavior at runtime. Each module requires its own sub-parameters.
#
# .. rubric::
#   - Determine available options for modules by running:
#       > seisflows print modules
#   - Auto-fill with docstrings and default values (recommended) by running:
#       > seisflows configure
#   - Swap out module parameters for a configured parameter file by running:
#       > seisflows swap {module} {name} (e.g., seisflows swap solver specfem3d)
#   - To set values as NoneType, use: null
#   - To set values as infinity, use: inf
#
#                                    MODULES
#                                    ///////
# workflow (str):    The types and order of functions for running SeisFlows
# system (str):      Computer architecture of the system being used
# solver (str):      External numerical solver to use for waveform simulations
# preprocess (str):  Preprocessing schema for waveform data
# optimize (str):    Optimization algorithm for the inverse problem
# ==============================================================================
workflow: forward
system: workstation
solver: specfem2d
preprocess: default
optimize: gradient
"""


# SeisFlows 'Globe' logo in ASCII. Used for CLI print statements
ascii_logo = """
                                                                                
                             @@@@@@@@@@@@@@@@@@@@@@@                            
                       @@@@@@@@@.         .(%@&#(     %@@@.                     
                  /@@@@@@@&      *@@@@@@@@@&         ,@@@@  @@(                 
               @@@@@@@@      @@@@@@@@       &@@@@@@@@@@@@@  .@% @@              
            *@@@@@@@      @@@@@@@      (@@@@@@               @@ @@ @*           
          @@@@@@@      @@@@@@@,     /@@@@@                       @ @ @@         
        &@@@@@@      @@@@@@@.     @@@@@@                           @ @ @&       
       @@@@@@@      @@@@@@@      @@@@@@                             @@ @ @      
     %@@@@@@      @@@@@@@@      @@@@@@                               @@ @ @&    
    @@@@@@@      @@@@@@@@      @@@@@@                                 @* @ (@   
   @@@@@@@       @@@@@@@       @@@@@@                                 @@ @@ @@  
   @@@@@@&      @@@@@@@@      ,@@@@@@/                               .@@ .@& @  
  @@@@@@@       @@@@@@@@       @@@@@@@                               @@@  @@ ,@ 
  @@@@@@@       @@@@@@@@       @@@@@@@@                             @@@  (@@  @ 
 @@@@@@@@       @@@@@@@@.       @@@@@@@@@                         @@@@   @@@  @@
 @@@@@@@@       @@@@@@@@@        @@@@@@@@@@@                   @@@@@@   @@@@  @@
 @@@@@@@@       ,@@@@@@@@@         @@@@@@@@@@@@@%         %@@@@@@@@.   @@@@   @@
  @@@@@@@        @@@@@@@@@@*         #@@@@@@@@@@@@@@@@@@@@@@@@@@@     @@@@   @@ 
  @@@@@@@@        @@@@@@@@@@@            @@@@@@@@@@@@@@@@@@@@&      @@@@@    @@ 
   @@@@@@@@        *@@@@@@@@@@@*               .@@@@@@&           @@@@@@    @@  
   &@@@@@@@@         @@@@@@@@@@@@@&                            @@@@@@@(    @@@  
    @@@@@@@@@          @@@@@@@@@@@@@@@@                   %@@@@@@@@@%    #@@@   
     (@@@@@@@@@          *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      @@@#    
       @@@@@@@@@#            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       @@@@      
        (@@@@@@@@@@              @@@@@@@@@@@@@@@@@@@@@@@@@%         @@@@#       
          &@@@@@@@@@@%                   #@@@@@@@@(              @@@@@&         
             @@@@@@@@@@@@                                    %@@@@@@.           
               &@@@@@@@@@@@@@@                          *@@@@@@@@@              
                   @@@@@@@@@@@@@@@@@@@@&/.   ,%@@@@@@@@@@@@@@@                  
                       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                      
                             &@@@@@@@@@@@@@@@@@@@@@&                                
    """

ascii_logo_small = """
                                    @@@@@@@@@@                        
                               .@@@@.    .%&(  %@.          
                            @@@@   @@@@   &@@@@@@ ,%@       
                         @@@@   @@@,  /@@              @    
                        @@@   @@@@   @@@              @     
                      @@@@   @@@@   @@@                @  @ 
                      @@@   @@@@   ,@@@                @ @  
                     @@@@   @@@@    @@@@              @@ @ @
                     @@@@   @@@@@    @@@@@          @@@ @@ @
                     @@@@    @@@@@     @@@@@@@@@@@@@@  @@  @
                      @@@@    @@@@@@        @@@&     @@@  @ 
                      @@@@@     @@@@@@@@         %@@@@#  @@ 
                        @@@@#      @@@@@@@@@@@@@@@@@   @@   
                         &@@@@@          @@@@(       @@&    
                            @@@@@@@             /@@@@       
                                @@@@@@@@@@@@@@@@@
                                    @@@@@@@@@@          
"""
