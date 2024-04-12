"""
A simple checkpointing system that is used for failure tolerance and seamless
restarts during workflows. Keeps track of both function names and system task
IDs so that functions which partially finish do not rerun already completed
tasks when resuming. Task IDs are tracked with special string notation that
provides some level of compression for the state file.

.. note::

    This checkpointing system is not meant to be used for large scales, it has
    been tested for up to 500 tasks but larger number of tasks may result in
    file I/O issues

.. rubric::

    state = State()
    state("test_function", ntask=60)
"""
import os
import json
import time


class State:
    """
    Simple checkpointing system for completed tasks and sub-tasks
    """
    def __init__(self, fid="sfstate.json", lock_file=".sflock"):
        """
        :type fid: str
        :param fid: full path and filename used to save the checkpoint file.
            File should have the file extension .json. The path to this file 
            will be used to save the lock file
        :type lock_file: str
        :param lock_file: filename used to create a lock file to ensure that
            only one process can write to the state file at a time. This file
            will be created in the same directory as `fid`
        """
        self.fid = os.path.abspath(fid)
        _path = os.path.dirname(self.fid)

        self._lock_file = os.path.join(_path, lock_file)
        # Written a the top of the state file and ignored by rest of class
        self._header = {
            "_header0": "=====================================================",
            "_header1": "               SEISFLOWS STATE FILE",
            "_header2": "Empty strings signify COMPLETED tasks",
            "_header3": "Numerical values signify INCOMPLETE tasks by TASK ID",
            "_header4": "====================================================="
            }

        self.checkpoint = {}
        if os.path.exists(self.fid):
            self.load()
        
    def __str__(self):
        """Return the string representation of the checkpoint"""
        _str_out = ""
        if self.checkpoint:
            _max_key = max([len(key) for key in self.checkpoint.keys()])
            _str_out += f"idx {'state':^{_max_key}}  tasks\n"
            _str_out += f"{'=' * len(_str_out)}\n"
            for i, (key, val) in enumerate(self.checkpoint.items()):
                # e.g., " 1 test_function: 1-2,4-7,9,11-14\n"
                _str_out += \
                    f"{i:<3} {key:<{_max_key}}  {self._arr_to_str(val)}\n"
        return _str_out
    
    def stage(self, name, ntasks):
        """
        Stage a state by inputting its name and the number of tasks that that
        state needs to complete. If the state already exists, then nothing
        will happen. In order to overwrite states, use the `reset` method
        """
        # Ensure we are using the most up to date checkpoint
        if os.path.exists(self.fid):
            self.load()

        if name not in self.checkpoint:
            self.checkpoint[name] = [0] * ntasks  # All tasks incomplete
            self.save()

    def get(self, name):
        """
        Get the status of a specific state by name. If the state does not exist,
        return None
        """
        self.load()
        if name in self.checkpoint:
            return self._arr_to_str(self.checkpoint[name])
        else:
            return None
                    
    def complete(self, name, taskid=None):
        """
        Update the status of a task to complete. The thought is that only 
        complete jobs need to signify status change, else the status stays 0.
        This must be run with lock because the chance is high that multiple 
        sub tasks finish at the same time, so we need to ensure they are
        updating one at a time and not all at once, otherwise only the fastest
        task gets to write
        """
        # If no task id is given, complete ALL tasks for the given state
        if taskid is None:
            self.checkpoint[name] = [1] * len(self.checkpoint[name])
            self.save()
        # If a single taskid is given, complete just that task id. Run with
        # lock file because this may be called by multiple processes at once
        else:
            def _complete_task():
                """Procedure for an individual task to set its status to 
                complete"""
                self.load()
                try:
                    self.checkpoint[name][taskid] = 1
                except IndexError:
                    import pdb;pdb.set_trace()
                self.save()

            self._run_with_lock(_complete_task)
        
    def check_complete(self, name):
        """
        Check if a task or specific task_id is complete. Return True if task
        is complete, False in all other cases
        """
        is_complete = False
        if name in self.checkpoint:
            is_complete = all(i == 1 for i in self.checkpoint[name])

        return is_complete
    
    def reset(self):
        """Reset the state file, usually for the next iteration"""
        self.checkpoint = {}
        self.save()

    def save(self):
        """
        Save the internal checkpoint to JSON file. Overwrite any existing
        Use a lock file to ensure that multiple processes do not try to write
        to this file at the same time, which would result in lost data from one
        of the processes
        """
        # Convert the arrays in the checkpoint to string representation for
        # better readability
        _to_save = {}
        # Add a header to data since JSON cannot have comments
        for key, val in self._header.items():
            _to_save[key] = val
        # Convert all arrays to string representations
        for key, val in self.checkpoint.items():
            _to_save[key] = self._arr_to_str(val)
        
        with open(self.fid, "w") as f:
            json.dump(_to_save, f, indent=4)

    def load(self):
        """Load the internal checkpoint from JSON file"""
        with open(self.fid, "r") as f:
            loaded_checkpoint = json.load(f)

        # Convert the string representations back to arrays
        self.checkpoint = {}
        for key, val in loaded_checkpoint.items():
            if key.startswith("_"):
                continue
            self.checkpoint[key] = self._str_to_arr(val)

    def lock(self):
        """Create a lock file to prevent multiple processes from writing to
        the state file at the same time"""
        os.open(self._lock_file, os.O_CREAT | os.O_EXCL)  # LOCK
    
    def unlock(self):
        """Remove the lock file to allow other processes to write"""
        os.remove(self._lock_file)

    def _run_with_lock(self, func, wait_time_s=0.01, failsafe_s=10, *args, 
                        **kwargs):
        """
        Only attempt to do a task when access to a specific lock file is 
        granted. This prevents multiple parallel processes from trying to
        do the same task at the same time

        .. note::

            https://stackoverflow.com/questions/56178988/\
                    python-script-to-writelock-read-a-file-between-3-processes

        :type func: function
        :param func: function to run when lock is granted, given *args, **kwargs
        :type wait_time_s: float
        :param wait_time_s: time in seconds to wait before trying to access the 
            lock file again if the previous attempt failed
        :type failsafe_s: float
        :param failsafe_s: time in seconds to wait before raising an exception
            if the lock file cannot be accessed
        """
        failsafe = 0
        try:    
            self.lock()
            func(*args, **kwargs)
            self.unlock()
        except FileExistsError:
            time.sleep(wait_time_s)
            failsafe += wait_time_s
            if failsafe > failsafe_s:
                raise Exception("run with lock exceeded failsafe time")
            self._run_with_lock(func, wait_time_s, failsafe_s, *args, **kwargs)


    @staticmethod
    def _arr_to_str(arr, check_val=0):
        """
        Convert a list containing 1's and 0's, signifying complete and 
        incomplete, to a string representation where all 0's are counted and
        all 1's are not. The string is based on SLURM array representation,
        where '-' represent spans, ',' separates individual tasks

        :type arr: list
        :param arr: list of 1's and 0's
        :type check_val: int
        :param check_val: value that signifies a task is incomplete, default 0
        :rtype: str
        :return: string representation of the array `arr`
        """
        # Collect all indices for values of 0 which are incomplete
        zero_indices = [idx for idx, val in enumerate(arr) if val == check_val]

        # If None, that means all jobs are complete, return empty
        if not zero_indices:
            return ""
        
        # First zero value starts us off
        str_out = str(zero_indices[0])
        running = False  # keeps track of runs of 0s
        # Index i tracks the previous index since we start incremeting at 1
        for i, jdx in enumerate(zero_indices[1:]):
            idx = zero_indices[i]
            if jdx - idx == 1:
                running = True
                # Edge case for the end of the array
                if jdx == zero_indices[-1]:
                    str_out += f"-{jdx}"
                continue
            else:
                if running:
                    str_out += f"-{idx},{jdx}"
                else:
                    str_out += f",{jdx}"
                running = False

        return str_out
    
    @staticmethod
    def _str_to_arr(str_, ntasks=None, check_val=0):
        """
        Convert a string representation of job status into a list of 1's and 0's
        Keeping the total number of tasks in mind

        E.g., "1-2,4-7,9,11-14" => [0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1]

        :type str_: str
        :param str_: string representation of array
        :type ntasks: int
        :param ntasks: number of tasks to expect and therefore the size of the
            array returned, if None, will be determined by the max value in the 
            `str_` which may not be correct
        :type check_val: int
        :param check_val: value that signifies a task is incomplete, default 0
        :rtype: list
        :return: list of 1's and 0's based on the string representation `str_`
        """
        # Opposite of check_val
        other_val = 1 if check_val == 0 else 0

        # Get default ntask by the max value, being careful that last val
        # may be a run, and dealing with empty `str_` meaning all tasks complete
        if ntasks is None:
            # Special case for empty `str_` meaning all tasks are complete, but
            # since we don't know the number of tasks, we return an empty list
            if str_ == "":
                return []
            # Get `ntask` by finding the end of the range
            val = str_.split(",")[-1]
            if "-" in val:
                ntasks = int(val.split("-")[-1]) + 1
            else:
                ntasks = int(val) + 1
        else:
            # Special case for empty `str_` meaning all tasks are complete
            if str_ == "":
                return [other_val] * ntasks 

        # Loop through string representation to generate array
        arr = [other_val] * (ntasks)
        for s in str_.split(","):
            if "-" in s:
                start, end = s.split("-")
                for i in range(int(start), int(end) + 1):
                    arr[i] = check_val 
            else:
                arr[int(s)] = check_val
                
        return arr

