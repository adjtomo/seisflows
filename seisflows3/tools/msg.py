"""
SeisFlows3 messages tool. For providing a uniform look to SeisFlows3 print
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
            if len(attrs["default"]) == 0:
                f.write(f"{key}: []\n")
            else:
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
