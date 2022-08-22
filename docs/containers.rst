Running SeisFlows in Containers
===============================

As part of the `NSF-funded SCOPED project
<https://www.nsf.gov/awardsearch/showAward?AWD_ID=2104052>`__ SeisFlows has
been containerized using `Docker <https://www.docker.com/>`__, which removes the
oftentimes troubling task of installing and configuring software. These
containers can be used to run SeisFlows on workstations (using Docker) or
on high performance computers using
`Singularity/Apptainer <https://apptainer.org/>`__.

.. note::
    The SeisFlows Docker image can be located here:
    https://github.com/SeisSCOPED/pyatoa

.. note::
    The `Pyatoa` container is shipped with the latest versions of
    `SeisFlows <https://github.com/adjtomo/seisflows>`__,
    `Pyatoa <https://github.com/adjtomo/pyatoa>`__,
    `PySEP <https://github.com/uafgeotools/pysep>`__  and
    their respective dependencies.

The remainder of this documentation page assumes you are familiar with Docker
and containers. To learn more about Docker and containers you can visit:
https://www.docker.com/resources/what-container/

To install Docker on your workstation, visit:
https://docs.docker.com/get-docker/


Running an example with Docker
------------------------------

Here we will step through how running the
`SeisFlows-SPECFEM2D <specfem2d_example.html>`__ example using the available
Docker image.

To get the latest version of the `Pyatoa` Docker image, run:

.. code-block:: bash
    $ docker pull ghcr.io/seisscoped/pyatoa:nightly

Their are two methods for running the SeisFlows example using this Docker image,
either through a `JupyterHub interface <https://jupyter.org/hub>`__, or through
the command line. The former provides a graphical user interface that mimics
a virtual desktop for easier navigation, while the latter provides more
flexibility for scripting and using SeisFlows as a research tool.

From JupyterHub
^^^^^^^^^^^^^^^

We can run SeisFlows through JupyterHub by running through a port:

.. code-block:: bash
    $ docker run -p 8888:8888 ghcr.io/seisscoped/pyatoa:nightly

To open the running JupyterHub instance, open the URL that is pasted to stdout
in your favorite web browser. This will likely look something like
http://127.0.0.1:8888/lab?token=xxx where xxx is your unique token

Once you have opened the JupyterHub instance, you should see a graphical
user interface. The Pyatoa, PySEP and SeisFlows repositories will be downloaded
and navigable in the file system on the left hand side of the window.

From inside the JupyterHub instance, click the 'Terminal' icon to open up a
terminal window, and run the following example command. Note, it's best to run
the example in a new directory to avoid muddling up the home directory.

.. code-block:: bash
    $ mkdir sf_sem2d_example
    $ cd sf_sem2d_example
    $ seisflows examples run 2  # when prompted, simply hit enter

This example will download, configure and compile SPECFEM2D into your
JupyterHub instance, and then run a SeisFlows-Pyatoa-SPECFEM2D inversion
problem.


From the command line
^^^^^^^^^^^^^^^^^^^^^

.. warning::
    Command line implementation does not currently work.

To run the SeisFlows SPECFEM2D example from the command line, we simply need
to point Docker to the image we just downloaded, and call the SeisFlows command
line tool. To run the help message:

.. code-block:: bash
    $ docker run ghcr.io/seisscoped/pyatoa:nightly seisflows -h

Running the example should be as easy as:

.. code-block:: bash
    $ docker run ghcr.io/seisscoped/pyatoa:nightly seisflows examples run 2

In the above example, SeisFlows automatically identifies the SPECFEM2D
installation within the Docker image and uses this as the external numerical
solver. All results will be output to your current working directory.

