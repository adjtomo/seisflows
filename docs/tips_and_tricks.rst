Tips and Tricks
================

Learn some neat tips and tricks for running SeisFlows that may not be evidently
apparent when running examples or looking at the parameter file.


Stopping Mid Workflow
----------------------

Stop a workflow prematurely to look at results or change parameters.
All SeisFlows workflows (except TestFlow) have a parameter called `stop_after` 
which can be used to stop mid workflow.

To check valid options for `stop_after`, run the following command **from 
within a valid working directory**.

.. code:: bash

   seisflows print tasks

To set your `stop_after` parameter you can use the `seisflows par` command. 
For example:

.. code:: bash

   seisflows par stop_after run_adjoint_simulations

To resume a stopped workflow, you only need to re-run `seisflows submit`. The
checkpointing system will ensure that the workflow picks up from where it left
off.

.. code:: bash

   seisflows submit


Checkpointing
-------------

SeisFlows has a checkpointing system which ensures that tasks that have already
been run will not be re-run in the case of job failures and workflow restarts. 

The checkpointing system uses a text file called `sfstate.txt` which simply has 
entries related to tasks in the task list.

Tasks in the task list have three states: 'completed', 'failed' and 'pending'.

- Completed: Task has already been run and will be skipped over if re-run
- Failed: Task has failed and will be re-run
- Pending: Task has not been executed and will be run 

SeisFlows manages the `sfstate.txt` file on its own, however Users can manually
edit the state file if they want certain tasks to be re-run. Simply open
the task file with a text editor and change states.

.. note::

   In the future we hope to improve the checkpointing system with a command 
   line option to edit the file `seisflows state`, and with a more sophisticated
   system that can single out particular job failures to re-run.

Tasktime vs. Walltime
---------------------

Jobs run on Clusters have two time-related parameters `tasktime` and `walltime`.

`Walltime` refers to the submission wall time given to the `main` job, whereas
`tasktime` refers to the submission wall time given to each simulation job.

`Tasktime` is relatively simple to figure out - it should be set to the longest
expected time it takes **one** simulation to finish. If running inversion 
workflows, expect that adjoint simulations will take longer to run w.r.t 
forward simulations. Be sure to add a little buffer time for serial processing 
steps taken before or after simulations.

`Walltime` should represent how long you think an **entire** workflow will take
to run. At an extreme, this can be set to the longest allowable walltime on 
your system (e.g., 24 hours). Or you can try to calculate how long an entire 
workflow will take.

For example, if you are running a 2 iteration inversion where each simulation
(tasktime) takes 10 min, then you may expect 1 forward simulation, 1 
adjoint simulation and 2-3 forward simulations for the line search. Given
open queues (i.e., all array jobs can run at the same time), this will equal
roughly 2 iterations * 5 simulations / iteration * 10 minutes / simulation 
= 100 minutes. 

In the above example, a User might want to add some buffer time for long 
queue times and non-simulation processing steps. An acceptable walltime might 
then be 150-200 minutes.
