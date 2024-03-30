"""
A simple checkpointing system that is used for failure tolerance and seamless
restarts during workflows. Keeps track of both function names and system task
IDs so that functions which partially finish do not rerun already completed
tasks when resuming

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
        self.checkpoint = {}
        self._lock_file = os.path.join(path, lock_file)

        if not os.path.exists(self.fid):
            self.save()
        else:
            self.checkpoint = self.load()

    def __call__(self, name, ntasks):
        """
        Main function that contains all of the logic of the State class so that
        Workflow's only need to run a single function to either checkpoint,
        or retrieve checkpoint data
        """
        # If tasks is compelte, then return True
        if self.check_complete(name):
            return True
        # Else, either set up the task or return tasks that need to be rerun
        else:
            if name in self.checkpoint:
                return self.checkpoint[name]
            else:
                self.set_task(name, ntasks=ntasks)
                return False
        
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
        with open(self.fid, "w") as f:
            json.dump(self.checkpoint, f, indent=4)

    def run_with_lock(self, func, wait_time_s=1, failsafe_s=60, *args, 
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
                os.open(self._lockfile, os.O_CREAT | os.O_EXCL)
                func(*args, **kwargs)
                os.remove(self._lock_file)
            except FileExistsError:
                time.sleep(wait_time_s)
                failsafe += wait_time_s
                if failsafe > failsafe_s:
                    raise Exception("run with lock exceeded failsafe time")

    def load(self):
        """Load the internal checkpoint from JSON file"""
        with open(self.fid, "r") as f:
            return json.load(f)

    def set_task(self, name, ntasks=1):
        """
        Set up a new task with all sub-tasks set to incomplete. Tasks are
        set where task ID corresponds to the index in the list, and the 
        corresponding value is either 0 for incomplete/failed/running etc. or 1 
        for complete/successful
        """
        self.checkpoint[name] = [0] * ntasks
        self.save()

    def end_task(self, name, taskid):
        """
        Update the status of a task to complete. The thought is that only 
        complete jobs need to signify status change, else the status stays 0
        """
        self.checkpoint[name][taskid] = 1
        self.save()

    
    