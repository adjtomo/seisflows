"""
A simple checkpointing system that is used for failure tolerance and seamless
restarts during workflows. Keeps track of both function names and system task
IDs so that functions which partially finish do not rerun already completed
tasks when resuming. Task IDs are tracked with special string notation that
provides some level of compression for the state file.

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
    def __init__(self, path="./", fid="sfstate.json", lock_file="sflock"):
        """
        :type path: str
        :param path: path to save the checkpoint file and the lock file used
            to make sure multiple processes do not write to the same file
        :type fid: str
        :param fid: file identifier for checkpoint file. filename only
        :type ftype: str
        :param ftype: file type for checkpoint file, either 'npz' for NumPy
            npz file, or 'txt' for text file
        :type ntasks: int
        :param ntasks: number of tasks to expect for each function that is
            checkpointed. should match the `system.ntasks` in SeisFlows 
        """
        self.fid = os.path.join(path, fid)
        self._lock_file = os.path.join(path, lock_file)
        # Written a the top of the state file and ignored by rest of class
        self._header = {
            "_header0": "=====================================================",
            "_header1": "               SEISFLOWS STATE FILE",
            "_header2": "empty strings signify that a task is COMPLETE",
            "_header3": "numerical values signify INCOMPLETE tasks by TASK ID",
            "_header4": "======================================================"
            }

        self.checkpoint = {}
        if os.path.exists(self.fid):
            self.load()

    def __call__(self, name, ntasks):
        """
        Main function that contains all of the logic of the State class so that
        Workflow's only need to run a single function to either checkpoint,
        or retrieve checkpoint data
        """
        # Ensure we are using the most up to date checkpoint
        if os.path.exists(self.fid):
            self.load(ntasks)

        # If tasks is complete, then return True
        if self.check_complete(name):
            status = [1] * ntasks
        # Else, either set up the task or return tasks that need to be rerun
        else:
            # Return array representation of tasks status
            if name in self.checkpoint:
                status = self.checkpoint[name]
            else:
                # Set a task which h
                status = [0] * ntasks
                self.checkpoint[name] = status
                self.save()
                
        return self._arr_to_str(status)
            
    def reset(self):
        """Reset the state file, usually for the next iteration"""
        self.checkpoint = {}
        self.save()
        
    def check_complete(self, name, taskid=None):
        """
        Check if a task or specific task_id is complete. Return True if task
        is complete, False in all other cases
        """
        is_complete = False
        if name in self.checkpoint:
            if taskid:
                is_complete = self.checkpoint[name][taskid] == 1
            else:
                is_complete = all(i == 1 for i in self.checkpoint[name])

        return is_complete

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

    def load(self, ntasks=None):
        """Load the internal checkpoint from JSON file"""
        with open(self.fid, "r") as f:
            loaded_checkpoint = json.load(f)

        # Convert the string representations back to arrays
        self.checkpoint = {}
        for key, val in loaded_checkpoint.items():
            if key.startswith("_"):
                continue
            self.checkpoint[key] = self._str_to_arr(val, ntasks)

    def done(self, name, taskid):
        """
        Update the status of a task to complete. The thought is that only 
        complete jobs need to signify status change, else the status stays 0.
        This must be run with lock because the chance is high that multiple 
        sub tasks finish at the same time, so we need to ensure they are
        updating one at a time and not all at once, otherwise only the fastest
        task gets to write
        """
        def _complete_task():
            """Procedure for an individual task to set its status to complete"""
            self.load()
            self.checkpoint[name][taskid] = 1
            self.save()

        self._run_with_lock(_complete_task)

    def lock(self):
        """Create a lock file to prevent other processes from working"""
        os.open(self._lock_file, os.O_CREAT | os.O_EXCL)

    def unlock(self):
        """Remove the lock file, unlocking the State class to allow new procs"""
        os.remove(self._lock_file)

    def _run_with_lock(self, func, wait_time_s=1, failsafe_s=60, *args, 
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
        while True:
            try:    
                self.lock()
                func(*args, **kwargs)
                self.unlock()
            except FileExistsError:
                time.sleep(wait_time_s)
                failsafe += wait_time_s
                if failsafe > failsafe_s:
                    raise Exception("run with lock exceeded failsafe time")

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
                ntasks = int(val.split("-")[-1])
            else:
                ntasks = int(val)
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
                try:
                    arr[int(s)] = check_val
                except IndexError:
                    import pdb; pdb.set_trace()
                
        return arr

